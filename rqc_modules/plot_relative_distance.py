import pandas
import ast
import scipy.stats
import matplotlib.pyplot as plt


def plot_relative_distance(args):
    DISTANCE = args.distance
    REFERENCE_LABEL = args.label
    PLOT_COUNTS = False  # TODO: make this an argument
    INPUT = args.input
    OUTPUT = args.output

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
        df[k] = df[k].apply(lambda s: list(ast.literal_eval(s)))

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
                d_total_offsets[k] = coverages[k]
            else:
                d_total_offsets[k] += coverages[k]

    for k, v in d_total_offsets.items():
        d_offset_hists[k] = [v.count(i) for i in range(-DISTANCE, DISTANCE)]

        kernel = scipy.stats.gaussian_kde(v)
        kde = kernel(x_ticks)

        d_offset_kdes[k] = kde

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
        'TTTN_count': 'blue'
    }

    fig, axes = plt.subplots()

    for k, v in d_offset_hists.items():
        if PLOT_COUNTS:
            normalised_v_by_reference_count = v
        else:
            normalised_v_by_reference_count = [x / len(df) * 100 for x in v]
        # axes.plot(x_ticks, v, label=k, color=d_colors[k])
        axes.plot(x_ticks, normalised_v_by_reference_count, label=k, color=d_colors[k])
        axes.fill_between(x_ticks, normalised_v_by_reference_count, alpha=0.2, color=d_colors[k])

    # TODO HACK
    CUSTOM_Y_MAX = 60

    axes.axvline(x=0, color='grey', label=REFERENCE_LABEL, ls="--", linewidth=1.0)
    # axes.set_ylabel('count')
    axes.set_ylabel('site frequency (%)')

    axes.set_xlabel('offset from {} (nt)'.format(REFERENCE_LABEL))
    axes.set_xlim(xmin=-DISTANCE, xmax=DISTANCE)
    axes.set_ylim(ymin=0)

    if CUSTOM_Y_MAX:
        axes.set_ylim(ymin=0, ymax=CUSTOM_Y_MAX)
    axes.legend()
    plt.legend(loc="upper right")

    # KDE plot
    # fig, axes = plt.subplots()

    # for k, v in d_offset_kdes.items():
    #     axes.plot(x_ticks, v, label=k, color=d_colors[k])
    #     axes.fill_between(x_ticks, v, alpha=0.2, color=d_colors[k])

    # axes.axvline(x=0, color='grey', label=REFERENCE_LABEL, ls="--", linewidth=1.0)
    # axes.set_ylabel('count')
    # axes.set_xlabel('offset from {} (nt)'.format(REFERENCE_LABEL))
    # axes.set_xlim(xmin=-DISTANCE, xmax=DISTANCE)
    # axes.set_ylim(ymin=0)

    # # if CUSTOM_Y_MAX:
    #     # axes.set_ylim(ymin=0, ymax=CUSTOM_Y_MAX)
    # axes.legend()
    # plt.legend(loc="upper right")

    print("reference count: {}".format(len(df)))
    for k, v in d_total_offsets.items():
        print("{} within range ({}) count: {}".format(k, DISTANCE, len(v)))

    # plot distribution of pam site count around each reference point
    # HIST OF PAM COUNT PER GENE
    if OUTPUT:
        plt.savefig("plot_relative_distance_offsets_{}".format(OUTPUT), format=OUTPUT_FORMAT)

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
        axes.plot(num_gene_x_ticks, v[:max_pam_sites_range], label=k, color=d_colors[k])
        print(v[:max_pam_sites_range])
        axes.fill_between(num_gene_x_ticks, v[:max_pam_sites_range], alpha=0.2, color=d_colors[k])

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
        plt.savefig("plot_relative_distance_by_gene_{}".format(OUTPUT), format=OUTPUT_FORMAT)


    plt.show()


