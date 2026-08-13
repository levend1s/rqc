import re
import sys

import pandas
import pysam
import numpy

from sklearn.preprocessing import MultiLabelBinarizer

from rqc_modules.constants import PYSAM_MOD_TUPLES
from rqc_modules.utils import process_input_files, process_annotation_file

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform

from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt
from collections import Counter

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import matplotlib

import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

RANDOM_SEED = 42
numpy.random.seed(RANDOM_SEED)

def timeit(func):
    @wraps(func)
    def _wrapped(*args, **kwargs):
        t0 = time.perf_counter()
        res = func(*args, **kwargs)
        t1 = time.perf_counter()
        print(f"{func.__name__} took {t1-t0:.3f}s")
        return res
    return _wrapped

def merge_small_clusters_tree_aware(labels, Z, min_size, weights=None):
    """
    Tree-aware version: a small cluster is only ever merged into its
    nearest neighbor AS DEFINED BY THE DENDROGRAM (its sibling subtree
    at the point it first joins something else) — never into a cluster
    on a completely different branch, even if that cluster happens to
    be closer by raw distance.
    """
    labels = labels.copy()
    n = len(labels)  # number of leaves == number of unique patterns

    if weights is None:
        weights = {i: 1 for i in range(n)}

    # --- Precompute tree structure from Z (iterative, no recursion) ---
    children = {}
    for i in range(len(Z)):
        l, r = int(Z[i, 0]), int(Z[i, 1])
        children[n + i] = (l, r)

    leaves_under = [None] * (2 * n - 1)
    for i in range(n):
        leaves_under[i] = frozenset([i])
    for i in range(len(Z)):
        l, r = int(Z[i, 0]), int(Z[i, 1])
        leaves_under[n + i] = leaves_under[l] | leaves_under[r]

    parent = [None] * (2 * n - 1)
    for i in range(len(Z)):
        l, r = int(Z[i, 0]), int(Z[i, 1])
        node = n + i
        parent[l] = node
        parent[r] = node

    # Map each exact leaf-set to the node representing it (only valid for
    # nodes that ARE pure subtrees — every node in a tree qualifies).
    leafset_to_node = {leaves_under[node]: node for node in range(2 * n - 1)}

    def cluster_weight_map():
        w = {}
        for cl in numpy.unique(labels):
            idx = numpy.where(labels == cl)[0]
            w[cl] = sum(weights.get(i, 1) for i in idx)
        return w

    while True:
        cluster_weights = cluster_weight_map()
        small_clusters = [cl for cl, w in cluster_weights.items() if w < min_size]

        if not small_clusters:
            break
        if len(cluster_weights) <= 1:
            break

        small_cl = min(small_clusters, key=lambda cl: cluster_weights[cl])
        idx_small = numpy.where(labels == small_cl)[0]

        node = leafset_to_node.get(frozenset(idx_small.tolist()))

        if node is None:
            # Small cluster's members no longer form a pure subtree
            # (can happen after a few merges). Fall back to nearest
            # currently-existing ancestor subtree that fully contains it.
            candidate_nodes = [
                nd for nd, leafset in enumerate(leaves_under)
                if set(idx_small.tolist()) <= leafset
            ]
            node = min(candidate_nodes, key=lambda nd: len(leaves_under[nd])) if candidate_nodes else None

        target_label = None
        if node is not None and parent[node] is not None:
            parent_node = parent[node]
            l, r = children[parent_node]
            sibling = r if l == node else l
            sibling_leaves = leaves_under[sibling]

            label_weights = {}
            for i in sibling_leaves:
                lbl = labels[i]
                if lbl == small_cl:
                    continue  # ignore any of the small cluster's own leaves if present
                label_weights[lbl] = label_weights.get(lbl, 0) + weights.get(i, 1)

            if label_weights:
                target_label = max(label_weights, key=label_weights.get)

        if target_label is None:
            # No usable sibling (e.g. hit the root, or sibling is entirely
            # the same small cluster) — nothing left to merge into.
            break

        labels[idx_small] = target_label

    return labels


from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import itertools
import seaborn as sns


def find_correlated_events_multivariate(rows, target_prefixes=None, predictor_prefixes=None,
                                          min_count=5, fdr_threshold=0.05, show_heatmap=False):
    """
    rows : list[list[str]]
        Per-read token lists that have ALREADY been filtered upstream
        (e.g. by min_feature_freq in run_pairwise_clustering). This
        function does not re-filter by frequency — it just builds the
        binary matrix directly from what's passed in.

    min_count : int
        Minimum reads a token must appear in to be usable as a target
        or predictor in the regression (guards against unstable fits on
        near-constant columns even after upstream filtering). Set to 0
        or 1 to skip this safety check entirely.
    """
    mlb_full = MultiLabelBinarizer()
    X_full = mlb_full.fit_transform(rows)
    tokens = numpy.array(mlb_full.classes_)

    n_reads = len(rows)
    token_counts = X_full.sum(axis=0)

    # --- print every token and its count, sorted descending ---
    print(f"\nToken counts (n_reads={n_reads}):")
    order = numpy.argsort(-token_counts)
    for idx in order:
        tok = tokens[idx]
        cnt = token_counts[idx]
        print(f"  {tok:40s} {cnt:6d}  ({cnt/n_reads:.2%})")
    print()

    target_idx = [i for i, t in enumerate(tokens) if target_prefixes is None or t.startswith(target_prefixes)]

    results = []
    for ti in target_idx:
        target_tok = tokens[ti]
        y = X_full[:, ti]

        predictor_idx = [
            i for i in range(len(tokens))
            if i != ti and (predictor_prefixes is None or tokens[i].startswith(predictor_prefixes))
        ]
        if len(predictor_idx) == 0:
            continue

        Xp = X_full[:, predictor_idx]
        pred_toks = tokens[predictor_idx]

        # drop predictors with ~zero variance among these reads (all 0s or all 1s)
        var_mask = (Xp.sum(axis=0) >= min_count) & (Xp.sum(axis=0) <= len(y) - min_count)
        if var_mask.sum() == 0:
            continue
        Xp = Xp[:, var_mask]
        pred_toks = pred_toks[var_mask]

        Xp_const = sm.add_constant(Xp)

        try:
            model = sm.Logit(y, Xp_const)
            fit = model.fit(disp=0, method="newton", maxiter=100)
        except Exception:
            try:
                fit = model.fit_regularized(disp=0, alpha=0.1, maxiter=200)
            except Exception:
                continue

        params = fit.params[1:]
        pvals = getattr(fit, "pvalues", None)
        pvals = pvals[1:] if pvals is not None else numpy.full(len(params), numpy.nan)

        for pt, coef, p in zip(pred_toks, params, pvals):
            results.append({
                "target": target_tok,
                "predictor": pt,
                "coef": coef,
                "odds_ratio": numpy.exp(coef),
                "p_value": p,
            })

    if not results:
        print("No valid target/predictor combinations — nothing to test.")
        return pandas.DataFrame(), pandas.DataFrame()

    res_df = pandas.DataFrame(results)
    valid_p = res_df["p_value"].notna()
    res_df.loc[valid_p, "p_adj"] = multipletests(res_df.loc[valid_p, "p_value"], method="fdr_bh")[1]
    res_df = res_df.sort_values("p_adj")

    significant = res_df[res_df["p_adj"] < fdr_threshold].sort_values("coef")

    if show_heatmap and len(res_df) > 0:
        pivot = res_df.pivot_table(index="target", columns="predictor", values="coef")
        plt.figure(figsize=(max(6, pivot.shape[1] * 0.4), max(5, pivot.shape[0] * 0.4)))
        sns.heatmap(pivot, cmap="RdBu_r", center=0, xticklabels=True, yticklabels=True)
        plt.title("Adjusted association (logistic regression coefficient)\ncontrolling for all other features")
        plt.tight_layout()
        plt.show()

    return res_df, significant

# https://uc-r.github.io/hc_clustering
def run_pairwise_clustering(df, feature_cols, min_support, distance_threshold=0.1, min_feature_freq=0.01, show_dendrogram=False):
    rows = []

    for _, row in df.iterrows():
        tokens = []
        for col in feature_cols:
            if col == "introns":
                for value in row.get(col, []) or []:
                    start, end = value
                    tokens.append(f"{col}_{start}-{end}")
            else:
                for value in row.get(col, []) or []:
                    tokens.append(f"{col}_{value}")
        rows.append(tokens)

    # --- Drop rare features: keep only tokens present in >= min_feature_freq of reads ---
    min_token_rows = int(numpy.ceil(min_feature_freq * len(df)))
    token_row_counts = Counter()
    for tokens in rows:
        for t in set(tokens):  # count once per row, even if repeated within a row
            token_row_counts[t] += 1

    keep_tokens = {t for t, c in token_row_counts.items() if c >= min_token_rows}
    n_dropped = len(token_row_counts) - len(keep_tokens)
    # print(f"Dropping {n_dropped} of {len(token_row_counts)} tokens below {min_feature_freq:.1%} frequency "
    #       f"(< {min_token_rows} of {len(df)} reads)")

    rows = [[t for t in tokens if t in keep_tokens] for tokens in rows]

    print("Top tokens after filtering:")
    for token, count in token_row_counts.most_common(10):
        if token in keep_tokens:
            print(f"  {token}: {count} ({count/len(df):.2%})")

    # --- optional: test for correlated / mutually exclusive feature pairs ---
    check_correlations = True
    # NOTE only works for pairs of features
    if check_correlations:
        corr_all, corr_sig = find_correlated_events_multivariate(
            rows,
            target_prefixes="introns",
            predictor_prefixes="m6A",
            min_count=min_token_rows,   # reuse the SAME threshold as the rest of the function
            fdr_threshold=0.05,
            show_heatmap=True,
        )

    # --- Deduplicate: cluster unique patterns, weighted by how many rows share them ---
    pattern_key = [tuple(sorted(set(tokens))) for tokens in rows]
    unique_patterns = sorted(set(pattern_key))
    pattern_to_uid = {p: i for i, p in enumerate(unique_patterns)}
    row_to_uid = numpy.array([pattern_to_uid[p] for p in pattern_key])

    unique_rows = [list(p) for p in unique_patterns]

    mlb = MultiLabelBinarizer()
    X = mlb.fit_transform(unique_rows)
    dists = pdist(X, metric="jaccard")
    Z = linkage(dists, method="average")

    unique_labels_arr = fcluster(Z, t=distance_threshold, criterion="distance")

    # --- Enforce minimum cluster size (in terms of ORIGINAL row counts, not unique patterns) ---
    pattern_counts = Counter(row_to_uid)  # how many original rows each unique pattern represents

    unique_labels_arr = merge_small_clusters_tree_aware(
        unique_labels_arr, Z, min_support, weights=pattern_counts
    )

    # Map unique-pattern cluster labels back to every original row
    labels = unique_labels_arr[row_to_uid]

    unique_labels = sorted(set(labels))
    relabel_map = {old: new for new, old in enumerate(unique_labels)}
    labels = numpy.array([relabel_map[lbl] for lbl in labels])
    unique_labels_arr = numpy.array([relabel_map[lbl] for lbl in unique_labels_arr])   # add this

    # --- build descriptive cluster names from top features ---
    top_n = 10

    def sanitize(tok):
        # keep names filesystem/column-friendly: no spaces, slashes, etc.
        return re.sub(r"[^0-9a-zA-Z_]+", "-", tok)

    cluster_names = {}
    cluster_summary = []
    for cluster_id in sorted(set(labels)):
        idx = [i for i, lbl in enumerate(labels) if lbl == cluster_id]
        counts = Counter()
        for i in idx:
            counts.update(set(rows[i]))
        top = counts.most_common(top_n)
        top_str = ", ".join(f"{tok} ({c}/{len(idx)})" for tok, c in top)

        top_tokens_clean = [sanitize(tok) for tok, _ in top]
        cluster_name = "_".join([f"cluster{cluster_id}"] + top_tokens_clean)
        cluster_names[cluster_id] = cluster_name

        cluster_summary.append({
            "cluster": cluster_name,
            "n_samples": len(idx),
            "top_positions": top_str,
        })

    df_out = df.copy()
    df_out["combined_primary_itemset"] = [cluster_names[int(lbl)] for lbl in labels]

    # ------------------ plot dendrogram ---------------------
    if show_dendrogram:
        for i, count in pattern_counts.most_common(10):
            label = ", ".join(unique_rows[i])
            print(f"{count:4d}  {label}")
        n_leaves = X.shape[0]

        if n_leaves < 2:
            print(f"Skipping dendrogram: only {n_leaves} unique pattern(s)")
        else:
            sys.setrecursionlimit(max(1000, n_leaves + 10))

            # unique_labels_arr[i] = cluster id for leaf i (i.e. for unique_rows[i])
            # Build a color for each cluster id
            cluster_ids_sorted = sorted(set(unique_labels_arr))
            cmap = matplotlib.colormaps["tab20"].resampled(len(cluster_ids_sorted))            
            cluster_color_map = {
                cid: cm.colors.to_hex(cmap(i)) for i, cid in enumerate(cluster_ids_sorted)
            }
            default_color = "#808080"  # gray for links that span multiple clusters

            # Map each node in Z (leaves 0..n-1, then merged nodes n..2n-2) to the
            # set of original leaf cluster ids under it, so we know whether a link
            # is "pure" (all leaves belong to one cluster) or mixed.
            node_cluster_ids = {}
            for i in range(n_leaves):
                node_cluster_ids[i] = {unique_labels_arr[i]}

            for i, (left, right, dist, count) in enumerate(Z):
                node_id = n_leaves + i
                left_ids = node_cluster_ids[int(left)]
                right_ids = node_cluster_ids[int(right)]
                node_cluster_ids[node_id] = left_ids | right_ids

            def link_color_func(node_id):
                ids = node_cluster_ids[node_id]
                if len(ids) == 1:
                    return cluster_color_map[next(iter(ids))]
                return default_color

            leaf_labels = [", ".join(unique_rows[i]) for i in range(n_leaves)]

            ddata = dendrogram(
                Z,
                labels=leaf_labels,
                # no_labels=True,
                link_color_func=link_color_func,
            )
            plt.ylabel("Jaccard distance")
            plt.axhline(y=distance_threshold, color="red", linestyle="--", label="distance threshold")

            ax = plt.gca()
            leaf_order = ddata["leaves"]
            xticklabels = ax.get_xticklabels()
            for tick_label, orig_idx in zip(xticklabels, leaf_order):
                tick_label.set_color(cluster_color_map[unique_labels_arr[orig_idx]])

            for pos, orig_idx in enumerate(leaf_order):
                x = 5 + 10 * pos
                ax.plot(x, 0, marker="s", markersize=4,
                        color=cluster_color_map[unique_labels_arr[orig_idx]],
                        clip_on=False)

            # --- legend ---
            cluster_row_counts = Counter()
            for i, cid in enumerate(unique_labels_arr):
                cluster_row_counts[cid] += pattern_counts.get(i, 0)

            legend_handles = [
                mpatches.Patch(
                    color=cluster_color_map[cid],
                    label=f"(n={cluster_row_counts[cid]}) {cluster_names.get(cid, f'cluster_{cid}')}"
                )
                for cid in cluster_ids_sorted
            ]
            ax.legend(
                handles=legend_handles,
                loc="lower left",
                bbox_to_anchor=(0, 1.02),
                fontsize=8,
                title="Cluster",
            )

            # blue_leaf_positions = [pos for pos, orig_idx in enumerate(ddata["leaves"]) if unique_labels_arr[orig_idx] == 3]
            # print(blue_leaf_positions)

            plt.xlabel("")
            plt.tight_layout()  # reserve right 15% of figure for the legend

            plt.show()

    return df_out, pandas.DataFrame(cluster_summary)

def gather_read_entries_for_region(
    row,
    bam_labels,
    input_files,
    bam_handles,
    COVERAGE_PADDING,
    MIN_DELETION_LENGTH,
    MODS,
    PYSAM_MOD_THRESHOLD,
):
    """
    Extract read entries from all BAM labels for a given annotation `row`.

    Returns a list of read_entry lists suitable for constructing the
    `read_table` DataFrame (same format as the inline loop previously
    in `cluster_transcripts`).
    """
    read_entries = []

    for label in bam_labels:
        samfile_path = input_files[label]["path"]
        samfile = bam_handles[label]

        READS_IN_REGION = list(
            samfile.fetch(
                contig=row["seq_id"],
                start=row["start"] - COVERAGE_PADDING,
                stop=row["end"] + COVERAGE_PADDING,
            )
        )

        # filter reads
        for r in READS_IN_REGION:
            read_is_forward = r.is_forward

            if r.is_secondary or r.is_supplementary:
                continue
            if (row["strand"] == "+" and not read_is_forward) or (
                row["strand"] == "-" and read_is_forward
            ):
                continue

            mod_positions = []

            ref_pos = r.get_reference_positions(full_length=True)
            mb = r.modified_bases

            read_start = r.reference_start
            read_end = r.reference_end
            read_qualities = r.query_qualities
            read_cigartuples = r.cigartuples
            read_name = r.query_name
            read_query_length = r.query_length

            for mod in MODS:
                if read_is_forward:
                    pysam_mod_tuple_code = "{}_for".format(mod)
                else:
                    pysam_mod_tuple_code = "{}_rev".format(mod)

                genomic_mod_positions = []
                if pysam_mod_tuple_code in PYSAM_MOD_TUPLES:
                    mods_probs = mb.get(PYSAM_MOD_TUPLES[pysam_mod_tuple_code]) if mb else None

                    if mods_probs:
                        # keep only mod positions which are above mod prob threshold
                        read_mod_positions = [x[0] for x in mods_probs if x[1] >= PYSAM_MOD_THRESHOLD]
                        genomic_mod_positions = [ref_pos[mod] for mod in read_mod_positions if ref_pos[mod] is not None]

                mod_positions.append(genomic_mod_positions)

            # introns
            introns = []
            ref_pos_iter = read_start

            for op, length in read_cigartuples:
                if op == 2:  # D =
                    if length >= MIN_DELETION_LENGTH:
                        introns.append((ref_pos_iter, ref_pos_iter + length))

                if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
                    ref_pos_iter += length

            # phred quality
            avg_quality = (sum(read_qualities) / len(read_qualities)) if read_qualities else 0

            poly_a_length = 0
            if r.has_tag("pt:i"):
                poly_a_length = r.get_tag("pt:i")

            read_entry = [
                read_name,
                label,
                # store the threshold used for this run (kept for provenance)
                PYSAM_MOD_THRESHOLD / 256.0,
                samfile_path,
                row["seq_id"],
                row["ID"],
                read_start,
                read_end,
                "+" if read_is_forward else "-",
                read_query_length,
                poly_a_length,
                avg_quality,
                introns,
            ]

            read_entry += mod_positions

            read_entries.append(read_entry)

    return read_entries

    # continuous variable encoding for feature clusters may not work, but for single categories (just m6A) could order by the midpoint/max of modification for position information along the transcript

def cluster_transcripts(args):
    INPUT = args.input
    ANNOTATION_FILE = args.annotation
    IDS = args.ids
    FEATURE_TYPE = args.type
    COVERAGE_PADDING = args.padding
    MOD_PROB_THRESHOLD = args.mod_prob_threshold
    OUTPUT_FILE = args.outfile
    MIN_DELETION_LENGTH = args.min_deletion_length
    CLUSTER_COLS = args.cluster_cols.split(',') if args.cluster_cols is not None else None
    MIN_CLUSTER_PERC = args.min_cluster_percent
    MIN_CLUSTER_SIZE_BULK = args.min_cluster_size_bulk
    DISTANCE_THRESHOLD = args.distance_threshold
    MIN_FEATURE_FREQ = args.min_feature_freq
    SHOW_DENDROGRAM = args.show_dendrogram
    MINIMUM_READS_TO_PROCESS = MIN_CLUSTER_SIZE_BULK
    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)
    MODS = [m for m in CLUSTER_COLS if m != "introns"]

    # -------------------- begin processing -------------------- # 

    # load annotation file and find indexes for all parent children
    gff_df = process_annotation_file(ANNOTATION_FILE)
    if COVERAGE_PADDING:
        gff_df["type"] = gff_df["type"].cat.add_categories(["{}bp".format(COVERAGE_PADDING)])

    if FEATURE_TYPE:
        matches = gff_df[gff_df['type'] == FEATURE_TYPE]
    else:
        matches = gff_df[gff_df['ID'].isin(IDS)]
    if matches.empty:
        print("ERROR: no matches found for type {} and ids {}".format(FEATURE_TYPE, IDS))

    print(matches)

    read_table_header = [
        "read_id",
        "label",
        "mod_prob_threshold",
        "bamfile_path",
        "contig",
        "ID",
        "read_start",
        "read_end",
        "read_strand",
        "read_length",
        "poly_a_length",
        "average_quality",
        "introns"
    ]

    for mod in MODS:
        read_table_header.append("{}".format(mod))

    input_files = process_input_files(INPUT)

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']

    # Open once: label -> handle
    bam_handles = {
        label: pysam.AlignmentFile(input_files[label]['path'], "rb")
        for label in bam_labels
    }

    all_count_tables = []

    try:
        # ------------------- MAIN READ AND FEATURE EXTRACTION LOOP FROM BAMFILES ------------------- #
        for _, row in matches.iterrows():
            print("processing {}...".format(row["ID"]))
            read_entries = []

            read_entries = gather_read_entries_for_region(
                row,
                bam_labels,
                input_files,
                bam_handles,
                COVERAGE_PADDING,
                MIN_DELETION_LENGTH,
                MODS,
                PYSAM_MOD_THRESHOLD,
            )

            read_table = pandas.DataFrame(read_entries, columns=read_table_header)

            # ------------------- CLUSTER PROCESSING (OPTIMIZED) ------------------- #
            df = read_table.reset_index(drop=True).copy()
            num_total_reads = len(df)

            try:
                if num_total_reads < MINIMUM_READS_TO_PROCESS:
                    raise ValueError("not enough reads to process!")

                if MIN_CLUSTER_PERC is not None:
                    min_cluster_size = max(MIN_CLUSTER_SIZE_BULK, int(numpy.ceil(MIN_CLUSTER_PERC * num_total_reads / len(bam_labels))))
                    print("USING MIN_CLUSTER_PERC {}: min_cluster_size = {}".format(MIN_CLUSTER_PERC, min_cluster_size))
                else:
                    min_cluster_size = MIN_CLUSTER_SIZE_BULK
                    print("USING MIN_CLUSTER_SIZE: min_cluster_size = {}".format(min_cluster_size))


                # So for a gene you can cluster by continuous variables (euclidean distance) or by categorical features (Jaccard distance) or by a combination of both (e.g. weighted sum of distances). The latter is experimental and may not work well, but it is possible to implement.
                df_clustered, closed_sets = run_pairwise_clustering(
                    df, CLUSTER_COLS, min_cluster_size, DISTANCE_THRESHOLD, MIN_FEATURE_FREQ, SHOW_DENDROGRAM
                )

                df_clustered["cluster"] = df_clustered["combined_primary_itemset"]
                
        
                # ---------- optional UMAP only for visualization ----------
                if len(matches) == 1:
                    print("generating UMAP embedding (viz only)...")

                    # emb = umap.UMAP(
                    #     n_neighbors=20,
                    #     min_dist=0.3,
                    #     metric="euclidean"
                    # ).fit_transform(df_clustered)
            
                    # df_clustered["umap_x"] = emb[:, 0]
                    # df_clustered["umap_y"] = emb[:, 1]
                    # print("DONE")
        
                    # ---------- write ----------
                    UMAP_OUTPUT_FILE = "{}.umap".format(OUTPUT_FILE)
                    df_clustered.to_csv(UMAP_OUTPUT_FILE, sep="\t", index=False)
                    print("UMAP_OUTPUT_FILE")
                    print(f"Done. Wrote: {UMAP_OUTPUT_FILE}")

                ct = (
                    df_clustered
                    .assign(row_key=lambda d: d["ID"].astype(str) + "_" + d["cluster"].astype(str))
                    .groupby(["row_key", "label"])
                    .size()
                    .unstack(fill_value=0)
                )
            except ValueError as e:
                print("ERROR:", e)
                # fallback: one row ID_clusterNA with per-label counts
                if df is None or df.empty:
                    ct = pandas.DataFrame(index=[f"{row['ID']}_clusterNA"])
                else:
                    ct = (
                        df.groupby("label")
                        .size()
                        .to_frame(name="n")
                        .T
                    )
                    ct.index = [f"{row['ID']}_clusterNA"]
            all_count_tables.append(ct)

            # if row["seq_id"] == "Pf3D7_02_v3":
            #     break
            
    finally:
        # Always close all handles
        for fh in bam_handles.values():
            fh.close()

    # after loop: main combined table
    main_count_table = (
        pandas.concat(all_count_tables, axis=0)
        .fillna(0)
        .astype(int)
        .sort_index()
    )

    # ---------- write ----------
    out_df = main_count_table.rename_axis("ID_cluster").reset_index()
    out_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"Done. Wrote: {OUTPUT_FILE}")
