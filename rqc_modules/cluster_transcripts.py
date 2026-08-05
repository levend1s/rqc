import pandas
import pysam
import numpy

from sklearn.feature_extraction import DictVectorizer
from sklearn.preprocessing import normalize
from sklearn.cluster import HDBSCAN

import scipy.sparse as sp

from rqc_modules.constants import PYSAM_MOD_TUPLES
from rqc_modules.utils import process_input_files, process_annotation_file

import ast
from itertools import combinations

import pandas as pd
from mlxtend.frequent_patterns import fpgrowth
from mlxtend.preprocessing import TransactionEncoder

def append_features_present_in_gt_percent(
    df: pandas.DataFrame,
    feature_columns: list[str],
    threshold: 0,
    new_suffix: str = "_gt1pct",
    verbose: bool = True,
):

    out = df.copy()
    kept_features_by_col = {}
    kept_counts_by_col = {}

    for col in feature_columns:
        exploded = (
            out[col]
            .apply(lambda x: x if isinstance(x, (list, tuple, set)) else [])
            .apply(lambda x: list(set(x)))  # de-duplicate within read
            .explode()
            .dropna()
        )

        read_counts = exploded.value_counts()  # feature -> #reads containing feature
        kept_counts = read_counts[read_counts > threshold].sort_values(ascending=False)

        keep_set = set(kept_counts.index)
        kept_features_by_col[col] = keep_set
        kept_counts_by_col[col] = kept_counts.to_dict()

        new_col = f"{col}{new_suffix}"
        out[new_col] = out[col].apply(
            lambda x: [f for f in (x if isinstance(x, (list, tuple, set)) else []) if f in keep_set]
        )

        # print(f"\nColumn: {col}")
        # print(f"Threshold: > {threshold:.2f} reads ({min_percent}%)")
        # if len(kept_counts) == 0:
        #     print("No features passed threshold.")
        # else:
        #     print("Kept features and counts:")
        #     for feat, cnt in kept_counts.items():
        #         print(f"  {feat}: {cnt}")

    return out, kept_features_by_col, kept_counts_by_col

"""
Frequent-pattern analysis for categorical modification positions.

Replaces the earlier "clustering" approach with proper frequent itemset
mining (FP-Growth), using *closed* itemsets rather than maximal itemsets
to avoid burying dominant sub-patterns, and comparing patterns across
categories (mod_cols) using lift rather than raw co-occurrence counts.
"""


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

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
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return set()
    raise TypeError(f"Unrecognized position value type: {type(value)}")


# ---------------------------------------------------------------------------
# Step 1: per-category closed itemset mining
# ---------------------------------------------------------------------------

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
        return pd.DataFrame(columns=["itemset", "support", "count"])

    encoder = TransactionEncoder()
    one_hot = encoder.fit(transactions).transform(transactions)
    one_hot_df = pd.DataFrame(one_hot, columns=encoder.columns_)

    frequent = fpgrowth(one_hot_df, min_support=min_support, use_colnames=True)
    if frequent.empty:
        return pd.DataFrame(columns=["itemset", "support", "count"])

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


def mine_closed_itemsets_per_category(df, mod_cols, min_support):
    """
    Run closed-itemset mining independently for each category in
    mod_cols (kept separate on purpose - see prior discussion: combining
    categories into one item universe blows up the search space and
    drowns cross-category patterns in intra-category noise).

    Returns: dict {col_name: closed_itemsets_df}
    """
    results = {}
    for col in mod_cols:
        transactions = df[col].apply(_parse_positions).tolist()
        results[col] = _mine_closed_itemsets(transactions, min_support)
    return results


# ---------------------------------------------------------------------------
# Step 2: multi-hot row membership (a row can satisfy several closed
# itemsets in the same category at once - there is no single "the"
# maximal itemset per row, so we don't force a single-label assignment)
# ---------------------------------------------------------------------------

def build_itemset_membership(df, mod_cols, closed_itemsets_by_col):
    """
    For every row and every category, record which closed itemsets that
    row satisfies (its position set is a superset of the itemset).

    Adds one column per category to df:
        "{col}_itemsets" -> list of matched itemset labels, e.g.
            ["mod1_1000-mod1_1010", "mod1_2000"]
        or ["__none__"] if the row has features in that category but
        none of them form a frequent closed pattern (this is the
        "noise" case), or ["__empty__"] if the row had no features at
        all in that category. Keeping these as explicit labels (rather
        than dropping the row) keeps counts consistent downstream.
    """
    df = df.copy()

    for col in mod_cols:
        closed_df = closed_itemsets_by_col.get(col, pd.DataFrame())
        # Largest itemsets first purely so labels read as "most specific
        # pattern first" if you inspect them; matching itself doesn't
        # require any particular order since we keep *all* matches.
        itemsets = sorted(
            closed_df["itemset"].tolist(), key=len, reverse=True
        ) if not closed_df.empty else []

        row_labels = []
        for positions in df[col].apply(_parse_positions):
            if not positions:
                row_labels.append(["__empty__"])
                continue

            matches = [
                "-".join(f"{col}_{p}" for p in sorted(itemset))
                for itemset in itemsets
                if itemset <= positions
            ]
            row_labels.append(matches if matches else ["__noise__"])

        df[f"{col}_itemsets"] = row_labels

    return df


# ---------------------------------------------------------------------------
# Step 3: cross-category comparison via lift, not raw joint counts
# ---------------------------------------------------------------------------

def cross_category_lift(df, mod_cols):
    """
    Compare itemsets *across* categories using lift:

        lift(A, B) = P(A and B) / (P(A) * P(B))

    Raw joint counts are misleading here because categories can have very
    different base rates (an itemset with 90% support vs one with 8%
    support) - a high joint count might just reflect the common one being
    common, not a real relationship. Lift > 1 means the two patterns
    co-occur more than you'd expect if they were independent; lift < 1
    means they co-occur less than expected (they tend to exclude each
    other); lift ~= 1 means no real relationship.

    Only compares itemsets between *different* categories (within-category
    relationships are already captured by the closed itemsets themselves).
    Rows are exploded per category so a row contributes to every itemset
    it matches (consistent with the multi-hot membership from step 2).

    Returns a DataFrame: col_a, itemset_a, col_b, itemset_b, support_a,
    support_b, joint_support, lift, joint_count.
    """
    n_rows = len(df)
    if n_rows == 0:
        return pd.DataFrame()

    # membership[col] -> {itemset_label: boolean pandas Series over rows}
    membership = {}
    for col in mod_cols:
        label_lists = df[f"{col}_itemsets"]
        labels_in_col = sorted(set(l for labels in label_lists for l in labels))
        membership[col] = {
            label: label_lists.apply(lambda labels, lbl=label: lbl in labels)
            for label in labels_in_col
            if label not in ("__empty__", "__noise__")
        }

    rows = []
    for col_a, col_b in combinations(mod_cols, 2):
        for label_a, mask_a in membership[col_a].items():
            support_a = mask_a.mean()
            for label_b, mask_b in membership[col_b].items():
                support_b = mask_b.mean()
                joint_mask = mask_a & mask_b
                joint_support = joint_mask.mean()
                joint_count = int(joint_mask.sum())
 
                denom = support_a * support_b
                lift = (joint_support / denom) if denom > 0 else float("nan")
 
                rows.append({
                    "col_a": col_a,
                    "itemset_a": label_a,
                    "col_b": col_b,
                    "itemset_b": label_b,
                    "support_a": support_a,
                    "support_b": support_b,
                    "joint_support": joint_support,
                    "joint_count": joint_count,
                    "lift": lift,
                })
 
    # Guard against the empty case (only one category passed in, or no
    # category had any itemset clear the support threshold) so we return
    # a properly-shaped empty frame instead of crashing on sort_values.
    columns = [
        "col_a", "itemset_a", "col_b", "itemset_b",
        "support_a", "support_b", "joint_support", "joint_count", "lift",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)
 
    return pd.DataFrame(rows)[columns].sort_values(
        "lift", ascending=False
    ).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def run_frequent_pattern_analysis(df, mod_cols, min_percent=1.0):
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
    min_support = min_percent / 100.0

    closed_itemsets_by_col = mine_closed_itemsets_per_category(df, mod_cols, min_support)
    print("hey!")
    print(closed_itemsets_by_col)
    df_out = build_itemset_membership(df, mod_cols, closed_itemsets_by_col)
    lift_df = cross_category_lift(df_out, mod_cols)

    return df_out, closed_itemsets_by_col, lift_df

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
    CLUSTER_INTRONS = args.cluster_introns
    CLUSTER_MODS = args.cluster_mods.split(',') if args.cluster_mods is not None else None

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)
    MODS = ['m6A', 'm5C', 'pseU', 'm6A_inosine']

    NUM_SAMPLES = 4
    MIN_CLUSTER_SIZE = 10 * NUM_SAMPLES
    MAX_CLUSTER_SIZE = 200 * NUM_SAMPLES
    MIN_CLUSTER_PERC = 0.01

    MINIMUM_READS_TO_PROCESS = 40

    RANDOM_SEED = 42
    numpy.random.seed(RANDOM_SEED)


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
        read_table_header.append("{}_positions".format(mod))

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

            for label in bam_labels:
                samfile_path = input_files[label]['path']
                samfile = bam_handles[label]

                READS_IN_REGION = list(samfile.fetch(
                                    contig=row['seq_id'], 
                                    start=row['start']-COVERAGE_PADDING, 
                                    stop=row['end']+COVERAGE_PADDING
                                ))
                # filter reads
                for i in range(len(READS_IN_REGION)):
                    r = READS_IN_REGION[i]

                    read_is_forward = r.is_forward

                    if r.is_secondary or r.is_supplementary:
                        continue
                    if (row["strand"] == "+" and not read_is_forward) or (row["strand"] == "-" and read_is_forward):
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
                            pysam_mod_tuple_code = '{}_for'.format(mod)
                        else:
                            pysam_mod_tuple_code = '{}_rev'.format(mod)

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
                    ref_pos = read_start

                    for op, length in read_cigartuples:
                        if op == 2:  # D = 
                            if length >= MIN_DELETION_LENGTH:
                                introns.append((ref_pos, ref_pos + length))

                        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
                            ref_pos += length

                    # phred quality
                    avg_quality = (
                        sum(read_qualities) / len(read_qualities)
                        if read_qualities
                        else 0
                    )                

                    poly_a_length = 0

                    if r.has_tag('pt:i'):
                        poly_a_length = r.get_tag('pt:i')

                    read_entry = [
                        read_name,
                        label,
                        MOD_PROB_THRESHOLD,
                        samfile_path,
                        row['seq_id'],
                        row['ID'],
                        read_start,
                        read_end,
                        '+' if read_is_forward else '-',
                        read_query_length,
                        poly_a_length,
                        avg_quality,
                        introns,
                    ]

                    read_entry += mod_positions

                    read_entries.append(read_entry)
                    
            read_table = pandas.DataFrame(read_entries, columns=read_table_header)

            # ------------------- CLUSTER PROCESSING (OPTIMIZED) ------------------- #
            df = read_table.reset_index(drop=True).copy()
            num_total_reads = len(df)

            try:
                if num_total_reads < MINIMUM_READS_TO_PROCESS:
                    raise ValueError("not enough reads to process!")
    
                min_cluster_size = min(
                    MAX_CLUSTER_SIZE,
                    max(MIN_CLUSTER_SIZE, int(numpy.ceil(MIN_CLUSTER_PERC * num_total_reads / NUM_SAMPLES)))
                )

                mod_col_names = []
                if CLUSTER_MODS is not None:
                    mod_col_names = ["{}_positions".format(m) for m in CLUSTER_MODS]


                df_clustered, closed_sets, lift = run_frequent_pattern_analysis(
                    df, mod_col_names, min_percent=1.0
                )

                # print(df_clustered)
                # print(closed_sets)
                # print(lift)
        
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
