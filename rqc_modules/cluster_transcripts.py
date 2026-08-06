import pandas
import pysam
import numpy

from sklearn.feature_extraction import DictVectorizer
from sklearn.preprocessing import normalize
from sklearn.cluster import HDBSCAN, DBSCAN

import scipy.sparse as sp

from rqc_modules.constants import PYSAM_MOD_TUPLES
from rqc_modules.utils import process_input_files, process_annotation_file

import ast
from itertools import combinations

import pandas as pd
from mlxtend.frequent_patterns import fpgrowth
from mlxtend.preprocessing import TransactionEncoder

import time
from functools import wraps

def timeit(func):
    @wraps(func)
    def _wrapped(*args, **kwargs):
        t0 = time.perf_counter()
        res = func(*args, **kwargs)
        t1 = time.perf_counter()
        print(f"{func.__name__} took {t1-t0:.3f}s")
        return res
    return _wrapped

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


def _intron_to_tokens(value):
    """Convert an introns value (list/tuple of (start,end) pairs or string) into
    a list of unique string tokens "start-end" suitable for mining/matching.
    """
    if isinstance(value, str):
        if not value.strip():
            return []
        parsed = ast.literal_eval(value)
    elif value is None or (isinstance(value, float) and pd.isna(value)):
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


def merge_rare_itemsets(closed_df, total_rows, min_count=None, min_percent=None, jaccard_threshold=0.5):
    """Merge low-count itemsets into similar popular itemsets using Jaccard overlap.

    This is a single-pass merge: rare itemsets (count < min_count) are
    assigned to the best popular itemset if their Jaccard >= threshold.
    Returns (merged_df, mapping) where mapping maps rare_label -> target_label or '__noise__'.
    """
    if closed_df is None or closed_df.empty:
        return closed_df, {}

    if min_count is None:
        if min_percent is None:
            raise ValueError("Either min_count or min_percent must be provided")
        min_count = int(numpy.ceil(min_percent * total_rows))

    def label_of(itemset):
        return "-".join(sorted(itemset))

    closed = closed_df.copy().reset_index(drop=True)
    closed["label"] = closed["itemset"].map(label_of)

    popular = closed[closed["count"] >= min_count].copy()
    rare = closed[closed["count"] < min_count].copy()

    mapping = {}
    if popular.empty:
        for _, r in rare.iterrows():
            mapping[label_of(r["itemset"])] = "__noise__"
        return closed[closed["count"] >= min_count], mapping

    pop_sets = popular["itemset"].tolist()
    pop_labels = popular["label"].tolist()

    for _, r in rare.iterrows():
        r_set = set(r["itemset"])
        best_score = 0.0
        best_idx = None
        for i, p in enumerate(pop_sets):
            inter = len(r_set & set(p))
            union = len(r_set | set(p))
            score = inter / union if union > 0 else 0.0
            if score > best_score:
                best_score = score
                best_idx = i

        r_label = label_of(r["itemset"])
        if best_score >= jaccard_threshold and best_idx is not None:
            target_label = pop_labels[best_idx]
            mapping[r_label] = target_label
            popular.loc[popular["label"] == target_label, "count"] += int(r["count"])
        else:
            mapping[r_label] = "__noise__"

    popular["support"] = popular["count"] / float(total_rows)
    merged = popular[["itemset", "support", "count"]].copy()
    merged = merged.sort_values("support", ascending=False).reset_index(drop=True)
    return merged, mapping


def merge_until_min_support(closed_df, total_rows, min_support, jaccard_threshold=0.5, max_iter=10):
    """Iteratively merge rare itemsets into popular ones until every
    remaining itemset has count >= min_count or no merges occur.

    Returns (merged_closed_df, mapping_total).
    """
    if closed_df is None or closed_df.empty:
        return closed_df, {}

    min_count = int(numpy.ceil(min_support * total_rows))
    current = closed_df.copy().reset_index(drop=True)
    mapping_total = {}

    for _ in range(max_iter):
        popular = current[current["count"] >= min_count].copy()
        rare = current[current["count"] < min_count].copy()
        if rare.empty:
            break
        if popular.empty:
            for _, r in rare.iterrows():
                mapping_total["-".join(sorted(r["itemset"]))] = "__noise__"
            break

        pop_sets = popular["itemset"].tolist()
        pop_labels = ["-".join(sorted(p)) for p in pop_sets]

        merged_any = False
        for _, r in rare.iterrows():
            r_set = set(r["itemset"])
            best_score = 0.0
            best_idx = None
            for i, p in enumerate(pop_sets):
                inter = len(r_set & set(p))
                union = len(r_set | set(p))
                score = inter / union if union > 0 else 0.0
                if score > best_score:
                    best_score = score
                    best_idx = i

            r_label = "-".join(sorted(r["itemset"]))
            if best_score >= jaccard_threshold and best_idx is not None:
                target_label = pop_labels[best_idx]
                mapping_total[r_label] = target_label
                popular.loc[popular["label"] == target_label, "count"] += int(r["count"])
                merged_any = True
            else:
                mapping_total[r_label] = "__noise__"

        if not merged_any:
            break

        popular["support"] = popular["count"] / float(total_rows)
        # keep unmerged rares
        remaining = []
        for _, r in rare.iterrows():
            r_label = "-".join(sorted(r["itemset"]))
            if mapping_total.get(r_label) == "__noise__":
                remaining.append(r)

        remaining_df = pd.DataFrame(remaining)
        if remaining_df.empty:
            current = popular[["itemset", "support", "count"]].copy().reset_index(drop=True)
        else:
            current = pd.concat([popular[["itemset", "support", "count"]], remaining_df[["itemset", "support", "count"]]], ignore_index=True)

    if not (current is None) and not current.empty:
        current = current.sort_values("support", ascending=False).reset_index(drop=True)

    return current, mapping_total


# ---------------------------------------------------------------------------
# Itemset -> rows mapping and greedy max-coverage pruning
# ---------------------------------------------------------------------------
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
                closed_df = closed.get(col, pd.DataFrame())
                if closed_df is None or closed_df.empty:
                    continue
                for _, r in closed_df.iterrows():
                    itemset = set(r["itemset"])
                    if itemset and itemset <= pos_set:
                        label = "-".join(f"{col}_{p}" for p in sorted(r["itemset"]))
                        itemset_rows[label].add(idx)

        return itemset_rows

    return itemset_rows



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

        # Build maps from textual label -> support and -> itemset size so
        # we can deterministically pick a primary itemset: prefer larger
        # itemsets, tie-break by higher support.
        support_map = {}
        size_map = {}
        if not closed_df.empty:
            for _, r in closed_df.iterrows():
                label = "-".join(f"{col}_{p}" for p in sorted(r["itemset"]))
                support_map[label] = float(r["support"]) if "support" in r else 0.0
                size_map[label] = len(r["itemset"]) if "itemset" in r else 0

        row_labels = []
        row_primary = []
        # For introns, normalize positions to the same string tokens used
        # during mining so set containment checks succeed.
        if col == "introns":
            positions_series = df[col].apply(_intron_to_tokens)
        else:
            positions_series = df[col].apply(_parse_positions)

        for positions in positions_series:
            if not positions:
                row_labels.append(["__empty__"])
                row_primary.append("__empty__")
                continue

            pos_set = set(positions)
            matches = [
                "-".join(f"{col}_{p}" for p in sorted(itemset))
                for itemset in itemsets
                if itemset <= pos_set
            ]

            if matches:
                # Choose the primary itemset by largest size (most specific)
                # and break ties by higher support.
                primary = max(
                    matches,
                    key=lambda lbl: (size_map.get(lbl, 0), support_map.get(lbl, 0.0))
                )
                row_labels.append(matches)
                row_primary.append(primary)
            else:
                row_labels.append(["__noise__"])
                row_primary.append("__noise__")

        df[f"{col}_itemsets"] = row_labels
        df[f"{col}_primary_itemset"] = row_primary

    return df


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

    closed = closed_df if closed_df is not None else pd.DataFrame()
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


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
# @timeit
def run_frequent_pattern_analysis(df, mod_cols, min_support):
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

    # closed_itemsets_by_col = mine_closed_itemsets_per_category(df, mod_cols, min_percent)
    # df_out = build_itemset_membership(df, mod_cols, closed_itemsets_by_col)

    closed_itemsets_by_col = mine_closed_itemsets_combined(df, mod_cols, min_support)  # for debugging / inspection only
    print(closed_itemsets_by_col)

    df_out = build_itemset_membership_combined(df, mod_cols, closed_itemsets_by_col)
    # an issue with rRNA 18S is that there would be many itemsets, but after finding membership many of those get depleted in

    # TODO: we need to figure out how to merge very low count itemsets into the most similar higher count itemset (e.g., by Jaccard overlap) so that we don't have a combinatorial explosion of itemsets with very low counts. This is a known issue with frequent pattern mining, and we need to implement a merging strategy to handle this.

    # Count occurrences of combined primary itemsets and drop any closed itemsets
    # whose primary label has a (negative) count below zero. This is defensive
    # (negative counts are unexpected) but mirrors the user's request.
    if "combined_primary_itemset" in df_out.columns:
        counts = df_out["combined_primary_itemset"].value_counts()
        print(counts)
        # find labels with count < 0 (rare/unexpected)
        to_remove = [lab for lab, c in counts.items() if c < min_support * len(df_out)]
        if to_remove:
            print("Removing closed itemsets with primary labels (count<min_support):", to_remove)
            # filter closed_itemsets_by_col to remove any matching labels
            def label_of_itemset(s):
                return "-".join(sorted(s))

            if closed_itemsets_by_col is not None and not closed_itemsets_by_col.empty:
                mask = ~closed_itemsets_by_col["itemset"].apply(lambda s: label_of_itemset(s) in set(to_remove))
                closed_itemsets_by_col = closed_itemsets_by_col[mask].reset_index(drop=True)

    print(closed_itemsets_by_col)


    return df_out, closed_itemsets_by_col


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
    MAX_CLUSTER_SIZE_IN_SAMPLE = 200

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
