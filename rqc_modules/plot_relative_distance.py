import pandas
import ast
import scipy.stats
import matplotlib.pyplot as plt
import numpy
import pandas
from statsmodels.stats.multitest import multipletests

def plot_relative_distance(args):
    DISTANCE = args.distance
    REFERENCE_LABEL = args.label
    PLOT_COUNTS = False  # TODO: make this an argument
    INPUT = args.input
    OUTPUT = args.output
    FOREGROUND_LABEL = args.foreground
    BACKGROUND_LABEL = args.background

    if OUTPUT:
        OUTPUT_FORMAT = OUTPUT.split(".")[-1] if OUTPUT else "png"

        if OUTPUT_FORMAT not in ["png", "eps", "pdf"]:
            raise ValueError("Output format must be one of: png, eps, pdf")

    print("LOG - plot_relative_distance started")
    print("LOG - reference label: {}".format(REFERENCE_LABEL))
    print("LOG - input file: {}".format(INPUT))
    offsets_file_path = INPUT

    if not offsets_file_path.endswith(".tsv"):
        raise ValueError("Input file must be a TSV file.")
    print("LOG - distance: {}".format(DISTANCE))
    print("LOG - label: {}".format(REFERENCE_LABEL))
    print("LOG - input file: {}".format(INPUT))

    offsets_file_path = INPUT
    df = pandas.read_csv(offsets_file_path, sep='\t')
    df = df.set_index('gene_id')
    keys = [x for x in df.columns if x not in ['gene_id', 'position'] and "count" not in x]
    for k in keys:
        print(df[k])
        def safe_eval(s):
            if pandas.isna(s):  # handles NaN safely
                return []
            return list(ast.literal_eval(s))
        df[k] = df[k].apply(lambda s: safe_eval(s))

    d_coverages = df.to_dict(orient='index')

    # TODO: write coverages as tsv file for plotting without recalculating offsets
    # load file into this variable d_coverages

    d_total_offsets = {}
    d_offset_hists = {}
    d_offset_kdes = {}
    x_ticks = list(range(-DISTANCE, DISTANCE))

    for k in keys:
        d_total_offsets[k] = []
        d_offset_hists[k] = []

    for idx, coverages in d_coverages.items():
        for k in keys:
            if len(d_total_offsets[k]) == 0:
                d_total_offsets[k] = list(coverages[k])
            else:
                d_total_offsets[k] += coverages[k]

    for k, v in d_total_offsets.items():
        d_offset_hists[k] = [v.count(i) for i in range(-DISTANCE, DISTANCE)]

        kernel = scipy.stats.gaussian_kde(v)
        kde = kernel(x_ticks)

        d_offset_kdes[k] = kde

    def gene_level_permutation_test_fast(df, foreground_key, background_key, distance,
                                      n_permutations=100000, seed=42,
                                      base_resolution=True, window_size=10):
        """
        Vectorized gene-level block permutation test for positional enrichment.

        For each gene, keeps the observed number of foreground sites fixed and
        redraws that many positions (without replacement) from the gene's own
        background position list, across ALL permutations at once per gene
        (avoids the n_permutations x n_genes Python-level loop).
        """
        rng = numpy.random.default_rng(seed)

        if base_resolution:
            bin_edges = numpy.arange(-distance, distance + 2)
        else:
            bin_edges = numpy.arange(-distance, distance + 1, window_size)

        n_bins = len(bin_edges) - 1

        # ---- observed histogram ----
        observed_offsets = []
        for offsets in df[foreground_key]:
            observed_offsets.extend(offsets)
        observed_hist, _ = numpy.histogram(observed_offsets, bins=bin_edges)

        # ---- pre-filter genes with no sites (wasted work otherwise) ----
        site_counts = df[foreground_key].apply(len)
        active = df[site_counts > 0]
        active_bg = active[background_key].tolist()
        active_n_sites = site_counts[site_counts > 0].tolist()

        print(f"Active genes: {len(active)} / {len(df)}")

        # ---- accumulate null histograms across all permutations ----
        null_hists = numpy.zeros((n_permutations, n_bins), dtype=numpy.int64)

        for a_positions, n_sites in zip(active_bg, active_n_sites):
            a_arr = numpy.asarray(a_positions)
            n_a = len(a_arr)
            if n_a == 0:
                continue

            n_draw = min(n_sites, n_a)

            if n_draw == n_a:
                # every A gets drawn every permutation - same fixed set, no
                # randomness needed, just bin it once and add to every row
                hist, _ = numpy.histogram(a_arr, bins=bin_edges)
                null_hists += hist  # broadcasts across all permutation rows
                continue

            # vectorized "without replacement" sampling for ALL permutations
            # at once: argsort of random keys gives a random permutation of
            # indices per row, take the first n_draw columns
            random_keys = rng.random((n_permutations, n_a))
            sampled_indices = numpy.argpartition(random_keys, n_draw, axis=1)[:, :n_draw]
            sampled_values = a_arr[sampled_indices]  # shape: (n_permutations, n_draw)

            # vectorized binning: map each sampled value to a bin index, then
            # tally counts per row without a Python-level per-row histogram call
            bin_idx = numpy.searchsorted(bin_edges, sampled_values, side='right') - 1
            valid = (bin_idx >= 0) & (bin_idx < n_bins)

            # flatten (row, bin) pairs and use bincount on a combined index to
            # get per-row, per-bin counts in one vectorized pass
            row_idx = numpy.repeat(numpy.arange(n_permutations), n_draw).reshape(n_permutations, n_draw)
            flat_bin = numpy.where(valid, row_idx * n_bins + bin_idx, -1).ravel()
            flat_bin = flat_bin[flat_bin >= 0]

            counts_flat = numpy.bincount(flat_bin, minlength=n_permutations * n_bins)
            null_hists += counts_flat.reshape(n_permutations, n_bins)

        # ---- p-values, expected counts, fold enrichment ----
        p_values = (numpy.sum(null_hists >= observed_hist, axis=0) + 1) / (n_permutations + 1)
        expected_counts = null_hists.mean(axis=0)
        fold_enrichment = numpy.divide(
            observed_hist, expected_counts,
            out=numpy.full(n_bins, numpy.nan), where=expected_counts > 0
        )

        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        results = pandas.DataFrame({
            "bin_center": bin_centers,
            "bin_start": bin_edges[:-1],
            "bin_end": bin_edges[1:],
            "observed_count": observed_hist,
            "expected_count": expected_counts,
            "fold_enrichment": fold_enrichment,
            "p_value": p_values,
        })

        _, q_values, _, _ = multipletests(results["p_value"], method="fdr_bh")
        results["q_value"] = q_values

        return results, null_hists

    results, null_hists = gene_level_permutation_test_fast(
        df, foreground_key="m6a", background_key="background_adenosines", distance=DISTANCE, n_permutations=1000 # ran 100,000 permutations for the final analysis, but this is slow for testing
    )
    _, q_values, _, _ = multipletests(results["p_value"], method="fdr_bh")
    results["q_value"] = q_values

    significant = results[results["q_value"] < 0.05].sort_values("fold_enrichment", ascending=False)
    print(significant.to_string(index=False))
    print("len non-significant: {}".format(len(results[results["q_value"] >= 0.05])))
    print("len significant: {}".format(len(results[results["q_value"] < 0.05])))

    d_num_pam_sites_hist = {}
    d_num_pam_sites = {}
    pam_count_keys = ["{}_count".format(x) for x in keys]
    for k in pam_count_keys:
        v = df[k].to_list()
        d_num_pam_sites_hist[k] = [v.count(i) for i in range(DISTANCE * 2)]
        d_num_pam_sites[k] = sorted(v)

    # print(d_coverages)
    # print(d_total_offsets)

    d_colors = {
        'NGG': 'green',
        'NGG_count': 'green',
        'TTTN': 'blue',
        'TTTN_count': 'blue',
        'm6a': 'purple',
        'm6a_count': 'purple',
        'annotation_three_p': 'green',
        'annotation_three_p_count': 'green'
    }

    import hashlib

    # Predefined palette of colors
    palette = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
            "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe"]

    def get_color(label):
        if label not in d_colors:
            # Use a hash of the label to get a deterministic index in the palette
            idx = int(hashlib.md5(label.encode()).hexdigest(), 16) % len(palette)
            d_colors[label] = palette[idx]
        return d_colors[label]

    fig, axes = plt.subplots()

    d_labels = {
        'annotation_three_p': 'annotated 3\' end',
        'm6a': 'canonical m6A',
        "NGG": 'NGG',
        "TTTN": 'TTTN'
    }

    def get_label(label):
        if label not in d_labels:
            # Use a hash of the label to get a deterministic index in the palette
            return label
        return d_labels[label]

    for k, v in d_offset_hists.items():
        if PLOT_COUNTS:
            normalised_v_by_reference_count = v
        else:
            normalised_v_by_reference_count = [x / len(df) * 100 for x in v]
        # axes.plot(x_ticks, v, label=k, color=d_colors[k])
        axes.plot(x_ticks, normalised_v_by_reference_count, label=k, color=get_color(k))
        axes.fill_between(x_ticks, normalised_v_by_reference_count, alpha=0.2, color=get_color(k))

    # TODO HACK
    CUSTOM_Y_MAX = None

    axes.axvline(x=0, color='grey', ls="--", linewidth=1.0)
    # axes.set_ylabel('count')
    axes.set_ylabel('site frequency (%)')

    axes.set_xlabel('offset from {} (nt)'.format(REFERENCE_LABEL))
    axes.set_xlim(xmin=-DISTANCE, xmax=DISTANCE)
    axes.set_ylim(ymin=0)

    if CUSTOM_Y_MAX:
        axes.set_ylim(ymin=0, ymax=CUSTOM_Y_MAX)
    axes.legend()
    plt.legend(loc="upper right")

    if OUTPUT:
        RELATIVE_OUTFILE = OUTPUT.replace(".{}".format(OUTPUT_FORMAT), "_offset_hist.{}".format(OUTPUT_FORMAT))
        plt.savefig(RELATIVE_OUTFILE, dpi=300, transparent=True, format=OUTPUT_FORMAT)

    # KDE plot
    fig, axes = plt.subplots()

    for k, v in d_offset_kdes.items():
        axes.plot(x_ticks, v, label=get_label(k), color=get_color(k))
        axes.fill_between(x_ticks, v, alpha=0.2, color=get_color(k))

    axes.axvline(x=0, color='grey', ls="--", linewidth=1.0)
    axes.set_ylabel('count')
    axes.set_xlabel('offset from {} (nt)'.format(REFERENCE_LABEL))
    axes.set_xlim(xmin=-DISTANCE, xmax=DISTANCE)
    axes.set_ylim(ymin=0)

    # if CUSTOM_Y_MAX:
        # axes.set_ylim(ymin=0, ymax=CUSTOM_Y_MAX)
    axes.legend()
    plt.legend(loc="upper right")

    print("reference count: {}".format(len(df)))
    for k, v in d_total_offsets.items():
        print("{} within range ({}) count: {}".format(k, DISTANCE, len(v)))

    # plot distribution of pam site count around each reference point
    # HIST OF PAM COUNT PER GENE
    if OUTPUT:
        DENSITY_OUTFILE = OUTPUT.replace(".{}".format(OUTPUT_FORMAT), "_offset_density.{}".format(OUTPUT_FORMAT))
        plt.savefig(DENSITY_OUTFILE, dpi=300, transparent=True, format=OUTPUT_FORMAT)

    fig, axes = plt.subplots()

    max_pam_count = 0
    for k, v in d_num_pam_sites.items():
        if max(v) > max_pam_count:
            max_pam_count = max(v)

    for k, v in d_num_pam_sites_hist.items():
        print(max_pam_count)
        max_pam_sites_range = int(max_pam_count * 1.1)
        max_pam_sites_range = 100
        num_gene_x_ticks = list(range(max_pam_sites_range))
        axes.plot(num_gene_x_ticks, v[:max_pam_sites_range], label=k, color=get_color(k))
        print(v[:max_pam_sites_range])
        axes.fill_between(num_gene_x_ticks, v[:max_pam_sites_range], alpha=0.2, color=get_color(k))

    axes.set_ylabel('count (genes)')
    axes.set_xlabel('number of PAM sites')
    # axes.set_xscale('log')
    axes.set_xlim(xmin=0, xmax=max_pam_sites_range)
    axes.set_ylim(ymin=0)


    # TODO HACK
    if CUSTOM_Y_MAX:
        axes.set_ylim(ymin=0, ymax=1200)
    axes.legend()
    plt.legend(loc="upper right")

    if OUTPUT:
        PLOT_OUTFILE = OUTPUT.replace(".{}".format(OUTPUT_FORMAT), "_count_per_reference.{}".format(OUTPUT_FORMAT))
        plt.savefig(PLOT_OUTFILE, dpi=300, transparent=True, format=OUTPUT_FORMAT)


    plt.show()


