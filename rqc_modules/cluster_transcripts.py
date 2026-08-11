import re

import pandas
import pysam
import numpy

from sklearn.feature_extraction import DictVectorizer
from sklearn.cluster import HDBSCAN
from sklearn.preprocessing import MultiLabelBinarizer

from rqc_modules.constants import PYSAM_MOD_TUPLES
from rqc_modules.utils import process_input_files, process_annotation_file

# import umap
import ast

from mlxtend.frequent_patterns import fpgrowth
from mlxtend.preprocessing import TransactionEncoder

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform

from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt
from collections import Counter

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


def merge_small_clusters(labels, dist_matrix, min_size):
    labels = labels.copy()
    counts = pandas.Series(labels).value_counts()
    small_clusters = counts[counts < min_size].index.tolist()

    for small_cl in small_clusters:
        idx_small = numpy.where(labels == small_cl)[0]
        if len(idx_small) == 0:
            continue  # already absorbed by a previous merge

        # Find the nearest *other* cluster (avg distance to each other cluster's points)
        other_labels = numpy.unique(labels[labels != small_cl])
        if len(other_labels) == 0:
            break  # only one cluster left, nothing to merge into

        best_cl, best_dist = None, numpy.inf
        for cl in other_labels:
            idx_other = numpy.where(labels == cl)[0]
            avg_dist = dist_matrix[numpy.ix_(idx_small, idx_other)].mean()
            if avg_dist < best_dist:
                best_dist, best_cl = avg_dist, cl

        labels[idx_small] = best_cl

    return labels

def run_pairwise_clustering(df, feature_cols, min_support):
    min_rows = int(numpy.ceil(min_support * len(df)))
    rows = []
    for _, row in df.iterrows():
        tokens = []

        for col in feature_cols:
            # handle introns column: list of (start, end) tuples
            if col == "introns":
                for value in row.get(col, []) or []:
                    start, end = value
                    tokens.append(f"{col}_{start}-{end}")
            else:
                for value in row.get(col, []) or []:
                    tokens.append(f"{col}_{value}")

        rows.append(tokens)

    # count how many rows each token appears in
    token_row_counts = Counter()
    for tokens in rows:
        for t in set(tokens):  # set() so a token counts once per row, even if repeated
            token_row_counts[t] += 1

    # keep only tokens meeting the row-count threshold
    keep_tokens = {t for t, c in token_row_counts.items() if c >= min_rows}
    rows = [[t for t in tokens if t in keep_tokens] for tokens in rows]

    mlb = MultiLabelBinarizer()
    X = mlb.fit_transform(rows)
    dists = pdist(X, metric="jaccard")
    Z = linkage(dists, method="average")

    # or look at the biggest jumps directly
    merge_dists = Z[:, 2]
    diffs = numpy.diff(merge_dists)
    knee_idx = numpy.argmax(diffs)  # index of biggest jump
    suggested_threshold = merge_dists[knee_idx]
    print(f"Suggested threshold from knee method: {suggested_threshold:.3f}")
    # suggested_threshold = 0.5  # HACK override for now, since the knee method is not working well

    labels = fcluster(Z, t=suggested_threshold, criterion="distance")

    # --- Enforce minimum cluster size by merging small clusters into nearest cluster ---
    dist_matrix = squareform(dists)  # full pairwise distance matrix
    labels = labels.copy()

    # labels = merge_small_clusters(labels, dist_matrix, min_rows)

    unique_labels = sorted(set(labels))
    relabel_map = {old: new for new, old in enumerate(unique_labels)}
    labels = numpy.array([relabel_map[lbl] for lbl in labels])

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
    print(pandas.DataFrame(cluster_summary).to_string(index=False))
    n_samples = X.shape[0]
    n_unique = len(numpy.unique(X, axis=0))
    import sys
    if n_samples < 2 or n_unique < 2:
        print(f"Skipping dendrogram: only {n_unique} unique pattern(s) across {n_samples} samples")
        ddata = None
    else:
        sys.setrecursionlimit(max(1000, n_samples + 10))  # avoid recursion limit error for large dendrograms
        ddata = dendrogram(Z, no_labels=True)
        plt.ylabel("Jaccard distance")
        plt.axhline(y=suggested_threshold, color="red", linestyle="--", label="suggested threshold")

        # scipy's dendrogram places leaves at x = 5, 15, 25, ... in the ORIGINAL
        # left-to-right order given by ddata["leaves"] (a list of original indices)
        leaf_order = ddata["leaves"]
        leaf_x = {orig_idx: 5 + 10 * pos for pos, orig_idx in enumerate(leaf_order)}

        # find the x-range (min/max) spanned by each cluster's leaves
        cluster_span = {}
        for orig_idx, x in leaf_x.items():
            cid = labels[orig_idx]
            cluster_span.setdefault(cid, []).append(x)

        ax = plt.gca()
        for cid, xs in cluster_span.items():
            mid_x = (min(xs) + max(xs)) / 2
            ax.text(mid_x, -0.02 * max(Z[:, 2]), f"cluster_{cid}",
                    ha="center", va="top", fontsize=9, rotation=90)

        plt.xlabel("")
        plt.tight_layout()
        plt.show()

    return df_out, None

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
    MIN_CLUSTER_PERC = args.min_cluster_perc

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)
    MODS = [m for m in CLUSTER_COLS if m != "introns"]

    MINIMUM_READS_TO_PROCESS = 40
    MIN_CLUSTER_SIZE_IN_SAMPLE = 10

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

                min_cluster_size = max(MIN_CLUSTER_SIZE_IN_SAMPLE * len(bam_labels), int(numpy.ceil(MIN_CLUSTER_PERC * num_total_reads / len(bam_labels))))
                min_support = min_cluster_size / num_total_reads

                # So for a gene you can cluster by continuous variables (euclidean distance) or by categorical features (Jaccard distance) or by a combination of both (e.g. weighted sum of distances). The latter is experimental and may not work well, but it is possible to implement.
                df_clustered, closed_sets = run_pairwise_clustering(
                    df, CLUSTER_COLS, min_support
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

            if row["seq_id"] == "Pf3D7_02_v3":
                break
            
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
