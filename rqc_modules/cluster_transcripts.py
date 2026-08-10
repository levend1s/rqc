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

def _parse_positions(value):
    """
    Row entries arrive as either an actual list/set, or a string
    representation like "[1000, 1010, 1020]" (e.g. read from a CSV).
    Normalize everything to a python set of ints so downstream set
    operations (intersection, subset checks) are cheap and consistent.
    """
    if isinstance(value, (list, set, tuple)):
        return set(value)
    if isinstance(value, str):
        if not value.strip():
            return set()
        return set(ast.literal_eval(value))
    if value is None or (isinstance(value, float) and pandas.isna(value)):
        return set()
    raise TypeError(f"Unrecognized position value type: {type(value)}")

def _intron_to_tokens(value):
    """Convert an introns value (list/tuple of (start,end) pairs or string) into
    a list of unique string tokens "start-end" suitable for mining/matching.
    """
    if isinstance(value, str):
        if not value.strip():
            return []
        parsed = ast.literal_eval(value)
    elif value is None or (isinstance(value, float) and pandas.isna(value)):
        return []
    else:
        parsed = value

    tokens = []
    for t in parsed:
        if isinstance(t, (list, tuple)) and len(t) >= 2:
            try:
                s = int(t[0])
                e = int(t[1])
                tokens.append(f"{s}-{e}")
            except Exception:
                tokens.append(str(t))
        else:
            tokens.append(str(t))
    # deduplicate while preserving order
    return list(dict.fromkeys(tokens))

def _mine_closed_itemsets(transactions, min_support):
    """
    Given a list of sets (one set of positions per row, for a single
    category), mine frequent itemsets with FP-Growth and reduce them to
    *closed* itemsets.

    An itemset X is closed if no proper superset of X has the same
    support. This is a middle ground between "all frequent itemsets"
    (redundant - lots of subsets carry identical information) and
    "maximal itemsets" (lossy - can discard a dominant, more general
    pattern just because a rarer superset happened to still clear the
    support threshold). Closed itemsets preserve every distinct support
    "step" in the lattice, so no genuinely different pattern gets
    silently merged away.

    Returns a DataFrame with columns: itemset (frozenset of positions),
    support (float fraction of rows), count (absolute row count).
    """
    if not any(transactions):
        return pandas.DataFrame(columns=["itemset", "support", "count"])

    encoder = TransactionEncoder()
    one_hot = encoder.fit(transactions).transform(transactions)
    one_hot_df = pandas.DataFrame(one_hot, columns=encoder.columns_)

    frequent = fpgrowth(one_hot_df, min_support=min_support, use_colnames=True)
    if frequent.empty:
        return pandas.DataFrame(columns=["itemset", "support", "count"])

    # Sort by itemset size descending so we can check "does any larger,
    # already-confirmed itemset contain this one with equal support"
    # without recomputing supersets from scratch each time.
    frequent = frequent.sort_values(
        by="itemsets", key=lambda s: s.map(len), ascending=False
    ).reset_index(drop=True)

    itemsets = frequent["itemsets"].tolist()
    supports = frequent["support"].tolist()

    is_closed = [True] * len(itemsets)
    for i, (small_set, small_support) in enumerate(zip(itemsets, supports)):
        for j, (big_set, big_support) in enumerate(zip(itemsets, supports)):
            if i == j:
                continue
            if len(big_set) > len(small_set) and small_set < big_set:
                if big_support == small_support:
                    is_closed[i] = False
                    break

    closed = frequent[is_closed].copy()
    closed["itemset"] = closed["itemsets"].apply(frozenset)
    closed["count"] = (closed["support"] * len(transactions)).round().astype(int)
    return closed[["itemset", "support", "count"]].sort_values(
        "support", ascending=False
    ).reset_index(drop=True)


    """
    Run closed-itemset mining independently for each category in
    mod_cols (kept separate on purpose - see prior discussion: combining
    categories into one item universe blows up the search space and
    drowns cross-category patterns in intra-category noise).

    Returns: dict {col_name: closed_itemsets_df}
    """
    results = {}
    for col in mod_cols:
        if col == "introns":
            transactions = df[col].apply(_intron_to_tokens).tolist()
        else:
            transactions = df[col].apply(_parse_positions).tolist()

        results[col] = _mine_closed_itemsets(transactions, min_support)
    return results

def _build_combined_transactions(df, mod_cols):
    """Build combined transactions across multiple columns.

    Each item is prefixed with the column name to keep identifiers
    unique (e.g. 'm6A_1000', 'introns_725853-726035'). Returns a list
    of lists (one list of tokens per row) suitable for the
    TransactionEncoder / _mine_closed_itemsets input.
    """
    transactions = []
    for _, row in df.iterrows():
        tokens = []
        for col in mod_cols:
            if col == "introns":
                items = _intron_to_tokens(row.get(col, []))
            else:
                # _parse_positions returns a set of ints/positions
                items = _parse_positions(row.get(col, []))

            for it in items:
                tokens.append(f"{col}_{it}")

        # deduplicate while preserving order
        transactions.append(list(dict.fromkeys(tokens)))
    return transactions

def mine_closed_itemsets_combined(df, mod_cols, min_support):
    """Mine closed itemsets on the combined token universe from all
    `mod_cols` together (not per-category). Items are column-prefixed so
    cross-category itemsets are possible.

    Returns a DataFrame with columns: itemset (frozenset of tokens),
    support (fraction), count (absolute rows).
    """
    transactions = _build_combined_transactions(df, mod_cols)
    # print("Mining combined closed itemsets over cols {} (rows={}, min_support={})".format(
    #     mod_cols, len(transactions), min_support
    # ))
    return _mine_closed_itemsets(transactions, min_support)

def compute_itemset_rows(df, mod_cols, closed, combined=True):
    """Compute mapping label -> set(row_indices) for itemsets.

    - If `combined` is True, `closed` is expected to be a single
      DataFrame (as returned by `mine_closed_itemsets_combined`) where
      each `itemset` is a frozenset of column-prefixed tokens.
    - If `combined` is False, `closed` is expected to be a dict
      {col: closed_df} as returned by `mine_closed_itemsets_per_category`.

    Returns a dict mapping textual label -> set(row_indices).
    """
    itemset_rows = {}
    if closed is None:
        return itemset_rows

    if combined:
        if closed.empty:
            return itemset_rows

        labels = ["-".join(sorted(x)) for x in closed["itemset"].tolist()]
        sets = [set(x) for x in closed["itemset"].tolist()]

        # initialize
        for lbl in labels:
            itemset_rows[lbl] = set()

        for idx, row in df.iterrows():
            tokens = []
            for col in mod_cols:
                if col == "introns":
                    items = _intron_to_tokens(row.get(col, []))
                else:
                    items = _parse_positions(row.get(col, []))

                for it in items:
                    tokens.append(f"{col}_{it}")

            pos_set = set(tokens)
            for lbl, s in zip(labels, sets):
                if s <= pos_set:
                    itemset_rows[lbl].add(idx)

        return itemset_rows

    # per-category dict input
    if isinstance(closed, dict):
        for col, closed_df in closed.items():
            if closed_df is None or closed_df.empty:
                continue
            for _, r in closed_df.iterrows():
                label = "-".join(f"{col}_{p}" for p in sorted(r["itemset"]))
                itemset_rows[label] = set()

        for idx, row in df.iterrows():
            for col in mod_cols:
                if col == "introns":
                    items = _intron_to_tokens(row.get(col, []))
                else:
                    items = _parse_positions(row.get(col, []))

                pos_set = set(items)
                closed_df = closed.get(col, pandas.DataFrame())
                if closed_df is None or closed_df.empty:
                    continue
                for _, r in closed_df.iterrows():
                    itemset = set(r["itemset"])
                    if itemset and itemset <= pos_set:
                        label = "-".join(f"{col}_{p}" for p in sorted(r["itemset"]))
                        itemset_rows[label].add(idx)

        return itemset_rows

    return itemset_rows

def build_itemset_membership_combined(df, mod_cols, closed_df):
    """Build combined-itemset membership for every row using a single
    closed-itemset dataframe mined over the column-prefixed token
    universe (see `mine_closed_itemsets_combined`). Adds two columns
    to the returned DataFrame:

      - `combined_itemsets`: list of matched itemset labels (joined tokens)
      - `combined_primary_itemset`: the single chosen primary label

    The primary selection prefers larger itemsets and breaks ties by
    higher support (same rule as the per-category version).
    """
    df = df.copy()

    closed = closed_df if closed_df is not None else pandas.DataFrame()
    itemsets = (
        sorted(closed["itemset"].tolist(), key=len, reverse=True)
        if not closed.empty
        else []
    )

    support_map = {}
    size_map = {}
    if not closed.empty:
        for _, r in closed.iterrows():
            # label is the dash-joined sorted tokens (tokens already include col_ prefix)
            label = "-".join(sorted(r["itemset"]))
            support_map[label] = float(r.get("support", 0.0))
            size_map[label] = len(r["itemset"]) if "itemset" in r else 0

    row_labels = []
    row_primary = []

    for _, row in df.iterrows():
        tokens = []
        for col in mod_cols:
            if col == "introns":
                items = _intron_to_tokens(row.get(col, []))
            else:
                items = _parse_positions(row.get(col, []))

            for it in items:
                tokens.append(f"{col}_{it}")

        # preserve order, deduplicate
        tokens = list(dict.fromkeys(tokens))

        if not tokens:
            row_labels.append(["__empty__"])
            row_primary.append("__empty__")
            continue

        pos_set = set(tokens)
        matches = [
            "-".join(sorted(itemset))
            for itemset in itemsets
            if itemset <= pos_set
        ]

        if matches:
            primary = max(
                matches, key=lambda lbl: (size_map.get(lbl, 0), support_map.get(lbl, 0.0))
            )
            row_labels.append(matches)
            row_primary.append(primary)
        else:
            row_labels.append(["__noise__"])
            row_primary.append("__noise__")

    df["combined_itemsets"] = row_labels
    df["combined_primary_itemset"] = row_primary
    return df

# @timeit
def run_frequent_pattern_analysis_backup(df, mod_cols, min_support):
    """
    Full pipeline:
      1. Mine closed frequent itemsets independently per category.
      2. Tag each row with every closed itemset it satisfies, per category
         (multi-hot, since a row can match more than one incomparable
         itemset).
      3. Compute cross-category lift between itemsets from different
         categories, so "does mod1 pattern A relate to mod2 pattern B"
         is answered with a proper independence-adjusted measure instead
         of a raw joint count.

    min_percent: minimum % of rows an itemset must appear in to be
    considered frequent (same semantics as the original min_reads_cluster
    threshold).

    Returns: (df_with_itemset_columns, closed_itemsets_by_col, lift_df)
    """
    closed_itemsets_by_col = mine_closed_itemsets_combined(df, mod_cols, min_support)  # for debugging / inspection only
    print(closed_itemsets_by_col)

    df_out = build_itemset_membership_combined(df, mod_cols, closed_itemsets_by_col)
    # an issue with rRNA 18S is that there would be many itemsets, but after finding membership many of those get depleted in

    # TODO: we need to figure out how to merge very low count itemsets into the most similar higher count itemset (e.g., by Jaccard overlap) so that we don't have a combinatorial explosion of itemsets with very low counts. This is a known issue with frequent pattern mining, and we need to implement a merging strategy to handle this.

    cluster_counts = df_out["combined_primary_itemset"].value_counts()
    print(cluster_counts)

    # For large datasets, run UMAP (jaccard) -> HDBSCAN to discover clusters
    try:
        nreads = len(df)
    except Exception:
        nreads = df.shape[0] if hasattr(df, "shape") else 0

    # TODO: this logic should be based on length of closed_itemsets
    if nreads > 5000:
        try:
            print(f"Running UMAP+HDBSCAN on {nreads} reads (min_support={min_support})")

            # build binary vectors from combined_itemsets
            vec = DictVectorizer(sparse=True)
            rows = []
            for lst in df_out.get("combined_itemsets", []):
                if not lst or (len(lst) == 1 and lst[0] in ("__empty__", "__noise__")):
                    rows.append({})
                else:
                    rows.append({lbl: 1 for lbl in lst})

            X = vec.fit_transform(rows)

            if X.shape[1] == 0:
                print("No itemset features found; skipping UMAP+HDBSCAN")
            else:
                min_cluster_size = int(numpy.ceil(min_support * nreads))
                min_samples = min(10, max(1, int(numpy.ceil(min_cluster_size / 10.0))))

                # UMAP embedding with Jaccard for binary set similarity
                reducer = umap.UMAP(n_components=10, metric="jaccard", random_state=42)
                emb = reducer.fit_transform(X)

                # HDBSCAN on embedding (euclidean)
                clusterer = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, metric="euclidean")
                labels = clusterer.fit_predict(emb)

                df_out["hdbscan_cluster"] = labels.astype(str)
                print("HDBSCAN clusters:", numpy.unique(labels, return_counts=True))

        except Exception as e:
            print("UMAP+HDBSCAN unavailable or failed:", str(e))
            print("Skipping HDBSCAN clustering. Install umap-learn and hdbscan to enable this.")

    # Merge small primary itemsets into larger popular itemsets.
    # Small = count < min_count (min_support * nrows).
    nrows = len(df)
    min_count = int(numpy.ceil(min_support * nrows))
    small_labels = [lab for lab, c in cluster_counts.items() if c < min_count]
    big_labels = [lab for lab, c in cluster_counts.items() if c >= min_count]
    print(f"Small labels (count < {min_count}): {small_labels}")
    print(f"Big labels (count >= {min_count}): {big_labels}")

    merged_mapping = {}
    removed_labels = []
    if small_labels:
        # prepare closed DF with textual labels
        if closed_itemsets_by_col is None or closed_itemsets_by_col.empty:
            print("No closed itemsets available to merge into; skipping merge for small labels.")
        else:
            closed_tmp = closed_itemsets_by_col.copy()
            closed_tmp["label"] = closed_tmp["itemset"].apply(lambda s: "-".join(sorted(s)))

            # popular candidates are those already above threshold
            popular = closed_tmp[closed_tmp["count"] >= min_count].copy()
            if popular.empty:
                print("No popular itemsets (count >= min_count); mapping small labels to __noise__")
                for lab in small_labels:
                    merged_mapping[lab] = "__noise__"
                    removed_labels.append(lab)
            else:
                # choose the global target as the largest itemset (tie-break by support)
                popular["size"] = popular["itemset"].apply(len)
                popular = popular.sort_values(["size", "support"], ascending=[False, False]).reset_index(drop=True)
                target_label = "-".join(sorted(popular.loc[0, "itemset"]))

                # accumulate increments: move counts of each small label into target
                total_inc = 0
                for lab in small_labels:
                    c = int(cluster_counts.get(lab, 0))
                    if c <= 0:
                        # still record removal but no increment
                        merged_mapping[lab] = target_label
                        removed_labels.append(lab)
                        continue
                    merged_mapping[lab] = target_label
                    total_inc += c
                    removed_labels.append(lab)

                # apply increment to target count and recompute support
                if total_inc > 0:
                    closed_tmp.loc[closed_tmp["label"] == target_label, "count"] += int(total_inc)
                    closed_tmp["support"] = closed_tmp["count"] / float(nrows)

                # drop removed small itemsets from closed_tmp
                closed_tmp = closed_tmp[~closed_tmp["label"].isin(removed_labels)].reset_index(drop=True)
                # keep only the canonical columns
                closed_itemsets_by_col = closed_tmp[["itemset", "support", "count"]].copy()

    # Add merged primary column to output (maps small labels -> chosen popular label)
    def _merge_map(lbl):
        return merged_mapping.get(lbl, lbl)

    df_out["combined_primary_itemset"] = df_out["combined_primary_itemset"].apply(_merge_map)

    return df_out, closed_itemsets_by_col

from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt
from collections import Counter

def run_frequent_pattern_analysis(df, mod_cols, min_support):
    min_rows = int(numpy.ceil(min_support * len(df)))
    print(min_rows)
    rows = []
    for _, row in df.iterrows():
        tokens = []
        for col in mod_cols:
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
    labels = fcluster(Z, t=max(1, min(10, len(df))), criterion="maxclust")

    df_out = df.copy()
    df_out["combined_primary_itemset"] = [f"cluster_{int(lbl)}" for lbl in labels]

    # --- summarize which positions define each cluster ---
    top_n = 3
    cluster_summary = []
    for cluster_id in sorted(set(labels)):
        idx = [i for i, lbl in enumerate(labels) if lbl == cluster_id]
        counts = Counter()
        for i in idx:
            counts.update(set(rows[i]))
        top = counts.most_common(top_n)
        top_str = ", ".join(f"{tok} ({c}/{len(idx)})" for tok, c in top)
        cluster_summary.append({
            "cluster": f"cluster_{cluster_id}",
            "n_samples": len(idx),
            "top_positions": top_str,
        })

    print(pandas.DataFrame(cluster_summary).to_string(index=False))

    # --- plot dendrogram ---
    # ddata = dendrogram(Z, no_labels=True)
    # plt.ylabel("Jaccard distance")

    # # scipy's dendrogram places leaves at x = 5, 15, 25, ... in the ORIGINAL
    # # left-to-right order given by ddata["leaves"] (a list of original indices)
    # leaf_order = ddata["leaves"]
    # leaf_x = {orig_idx: 5 + 10 * pos for pos, orig_idx in enumerate(leaf_order)}

    # # find the x-range (min/max) spanned by each cluster's leaves
    # cluster_span = {}
    # for orig_idx, x in leaf_x.items():
    #     cid = labels[orig_idx]
    #     cluster_span.setdefault(cid, []).append(x)

    # ax = plt.gca()
    # for cid, xs in cluster_span.items():
    #     mid_x = (min(xs) + max(xs)) / 2
    #     ax.text(mid_x, -0.02 * max(Z[:, 2]), f"cluster_{cid}",
    #              ha="center", va="top", fontsize=9, rotation=90)

    # plt.xlabel("")
    # plt.tight_layout()
    # plt.show()

    

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
                df_clustered, closed_sets = run_frequent_pattern_analysis(
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
                    .assign(row_key=lambda d: d["ID"].astype(str) + "_cluster" + d["cluster"].astype(str))
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
