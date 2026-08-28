import re
import sys

import pandas
import pysam
import numpy
from collections import Counter

from sklearn.preprocessing import MultiLabelBinarizer

from rqc_modules.constants import PYSAM_MOD_TUPLES
from rqc_modules.utils import process_input_files, process_annotation_file

from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import matplotlib

import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import seaborn as sns
from scipy.stats import fisher_exact, combine_pvalues
import itertools
from scipy.stats import spearmanr

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

def sanitize(tok):
    # keep names filesystem/column-friendly: no spaces, slashes, etc.
    return re.sub(r"[^0-9a-zA-Z_]+", "-", tok)

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


def find_mutually_exclusive_pairs(rows, allowed_tokens=None, target_prefixes=None,
                                    predictor_prefixes=None, min_individual_count=5):
    """
    Rank token pairs by mutual exclusivity, prioritizing pairs that are
    BOTH reasonably common individually AND never (or almost never)
    co-occur — rather than pairs that just happen to have few reads.

    Returns a dataframe sorted with the most exclusive pairs first.
    """
    if allowed_tokens is not None:
        rows = [[t for t in tokens if t in allowed_tokens] for tokens in rows]

    mlb_full = MultiLabelBinarizer()
    X_full = mlb_full.fit_transform(rows)
    tokens = numpy.array(mlb_full.classes_)
    n_reads = len(rows)

    target_idx = [i for i, t in enumerate(tokens) if target_prefixes is None or t.startswith(target_prefixes)]
    predictor_idx_all = [i for i, t in enumerate(tokens) if predictor_prefixes is None or t.startswith(predictor_prefixes)]

    results = []
    for ti in target_idx:
        for pi in predictor_idx_all:
            if ti == pi:
                continue

            a = X_full[:, ti]
            b = X_full[:, pi]

            n_a = int(a.sum())
            n_b = int(b.sum())
            n_both = int(((a == 1) & (b == 1)).sum())

            # both tokens need to individually be common enough to be meaningful
            if n_a < min_individual_count or n_b < min_individual_count:
                continue

            # expected co-occurrence under independence
            expected_both = (n_a * n_b) / n_reads

            # exclusivity score: how far actual co-occurrence falls short of
            # what you'd expect by chance, normalized by expectation.
            # 1.0 = perfectly exclusive (0 actual vs some expected)
            # 0.0 = exactly as expected by chance
            # negative = co-occurring MORE than expected
            exclusivity_score = (expected_both - n_both) / expected_both if expected_both > 0 else numpy.nan

            table = [[n_both, n_a - n_both], [n_b - n_both, n_reads - n_a - n_b + n_both]]
            _, p = fisher_exact(table)

            results.append({
                "token_a": tokens[ti],
                "token_b": tokens[pi],
                "n_a": n_a,
                "n_b": n_b,
                "n_both": n_both,
                "expected_both": round(expected_both, 2),
                "exclusivity_score": exclusivity_score,
                "p_value": p,
            })

    if not results:
        return pandas.DataFrame()

    res_df = pandas.DataFrame(results)
    # dedupe symmetric pairs (a,b) and (b,a) if target/predictor prefixes overlap
    res_df["pair_key"] = res_df.apply(lambda r: tuple(sorted([r["token_a"], r["token_b"]])), axis=1)
    res_df = res_df.drop_duplicates("pair_key").drop(columns="pair_key")

    res_df["p_adj"] = multipletests(res_df["p_value"], method="fdr_bh")[1]
    res_df = res_df.sort_values(["n_both", "exclusivity_score"], ascending=[True, False])

    return res_df

def collapse_redundant_predictor_combos(res_df, lift_col="lift_combo_absent",
                                          tolerance=0.05, group_col="intron"):
    """
    For each target (e.g. intron), collapse overlapping predictor
    combinations down to non-redundant ones. If a larger combo is a
    superset of an already-kept simpler one AND meaningfully improves
    on it, the simpler one is REMOVED (superseded) — only the larger,
    better combo is kept. If the larger combo does NOT improve
    meaningfully, IT is dropped instead, and the simpler one stays.
    """
    keep_rows = {}  # keyed by frozenset(tokens) -> (row, lift), acts as running "kept" set

    if group_col not in res_df.columns:
        return res_df

    for target, group in res_df.groupby(group_col):
        group = group.sort_values("n_predictor_tokens")  # simplest first
        kept = []  # list of dicts: {"tokens": set, "lift": float, "row": Series}

        for _, row in group.iterrows():
            this_tokens = set(row["predictor_combo"].split(" | "))
            this_lift = row[lift_col]

            # find any already-kept combo that this one is a superset of
            superseded = []
            is_redundant = False

            for k in kept:
                if k["tokens"].issubset(this_tokens) and k["tokens"] != this_tokens:
                    if this_lift > k["lift"] * (1 + tolerance):
                        # this combo meaningfully improves on the simpler kept one
                        # -> the simpler one is now superseded, mark for removal
                        superseded.append(k)
                    else:
                        # this combo does NOT meaningfully improve -> it's the redundant one
                        is_redundant = True
                        break

            if is_redundant:
                continue

            # remove any kept combos this one superseded
            for s in superseded:
                kept.remove(s)

            kept.append({"tokens": this_tokens, "lift": this_lift, "row": row})

        keep_rows[target] = [k["row"] for k in kept]

    all_rows = [row for rows in keep_rows.values() for row in rows]
    return pandas.DataFrame(all_rows)

def _get_m6a_position(token):
    """Extract the integer position from an 'm6A_<pos>' token."""
    return int(token.split("_", 1)[1])


def _get_intron_range(token):
    """Extract (start, end) ints from an 'introns_<start>-<end>' token."""
    coords = token.split("_", 1)[1]
    start_str, end_str = coords.split("-")
    return int(start_str), int(end_str)


def _within_distance_of_intron(predictor_token, intron_token, max_distance_nt):
    """
    True if predictor_token's position is within max_distance_nt of the
    intron's boundaries (0 if the predictor falls inside the intron span
    itself). Non-m6A / unparseable predictor tokens are always kept
    (distance filtering only applies to positional m6A predictors).
    """
    if max_distance_nt is None:
        return True
    if not predictor_token.startswith("m6A_"):
        return True

    try:
        pos = _get_m6a_position(predictor_token)
        start, end = _get_intron_range(intron_token)
    except (ValueError, IndexError):
        return True  # can't parse -> don't filter it out

    if pos < start:
        dist = start - pos
    elif pos > end:
        dist = pos - end
    else:
        dist = 0  # inside the intron span

    return dist <= max_distance_nt


def evaluate_cluster_tokens_as_intron_predictors(rows, cluster_important_tokens, target_prefix, min_count,
                                                    LIFT_THRESHOLD, max_distance_nt=None):
    """
    For each UNIQUE combination of non-intron predictor tokens across all
    clusters (deduplicated), test the union ("has any of these tokens")
    against EVERY intron token found across ALL clusters — but for each
    (predictor_combo, intron) pair, ONLY predictor tokens within
    max_distance_nt of that specific intron's boundaries are included in
    the union. Predictor tokens too far from a given intron are dropped
    for that intron's test (they may still be used for other, nearer
    introns).

    Parameters
    ----------
    max_distance_nt : int or None
        Maximum distance (nt) an m6A predictor token may be from an
        intron's start/end to be considered for that intron's test.
        None disables distance filtering entirely (original behavior).

    Returns
    -------
    pandas.DataFrame, one row per (predictor_combo, intron) pair, with
    the same columns as before, plus n_predictor_tokens now reflecting
    the DISTANCE-FILTERED predictor count actually used for that test.
    """
    mlb = MultiLabelBinarizer()
    X = mlb.fit_transform(rows)
    tokens = numpy.array(mlb.classes_)
    tok_to_idx = {t: i for i, t in enumerate(tokens)}
    n_reads = len(rows)

    all_intron_tokens = sorted({
        tok for imp_tokens in cluster_important_tokens.values()
        for tok in imp_tokens if tok.startswith(target_prefix)
    })

    if not all_intron_tokens:
        print("No intron tokens found among any cluster's important tokens.")
        return pandas.DataFrame()

    # --- deduplicate predictor sets across clusters (full set, before distance filtering) ---
    combo_to_clusters = {}
    for cid, imp_tokens in cluster_important_tokens.items():
        predictor_tokens = frozenset(t for t in imp_tokens if not t.startswith(target_prefix) and t in tok_to_idx)
        if not predictor_tokens:
            continue
        combo_to_clusters.setdefault(predictor_tokens, []).append(cid)

    results = []
    for predictor_set, cids in combo_to_clusters.items():
        for intron_tok in all_intron_tokens:
            if intron_tok not in tok_to_idx:
                continue

            # --- filter predictors to those near THIS intron ---
            filtered_predictor_tokens = sorted(
                t for t in predictor_set
                if _within_distance_of_intron(t, intron_tok, max_distance_nt)
            )
            if not filtered_predictor_tokens:
                continue  # nothing left close enough to this intron -> skip this pair

            pred_idx = [tok_to_idx[t] for t in filtered_predictor_tokens]
            p_combo = (X[:, pred_idx].sum(axis=1) >= 1).astype(int)

            n_present = int(p_combo.sum())
            n_absent = n_reads - n_present
            if n_present < min_count or n_absent < min_count:
                continue

            y = X[:, tok_to_idx[intron_tok]]
            p_intron_overall = y.mean()
            if p_intron_overall == 0:
                continue

            p_intron_given_present = y[p_combo == 1].mean()
            p_intron_given_absent = y[p_combo == 0].mean()

            lift_present = p_intron_given_present / p_intron_overall
            lift_absent = p_intron_given_absent / p_intron_overall

            n11 = int(((y == 1) & (p_combo == 1)).sum())
            n10 = int(((y == 1) & (p_combo == 0)).sum())
            n01 = int(((y == 0) & (p_combo == 1)).sum())
            n00 = int(((y == 0) & (p_combo == 0)).sum())
            _, p_value = fisher_exact([[n11, n10], [n01, n00]])

            results.append({
                "clusters": ", ".join(str(c) for c in sorted(cids, key=str)),
                "predictor_combo": " | ".join(filtered_predictor_tokens),
                "n_predictor_tokens": len(filtered_predictor_tokens),
                "intron": intron_tok,
                "p_intron_overall": p_intron_overall,
                "p_intron_given_combo_present": p_intron_given_present,
                "lift_combo_present": lift_present,
                "p_intron_given_combo_absent": p_intron_given_absent,
                "lift_combo_absent": lift_absent,
                "p_value": p_value,
            })

    res_df = pandas.DataFrame(results)
    if len(res_df) > 0:
        INV_LIFT_THRESHOLD = 1 / LIFT_THRESHOLD

        res_df["p_adj"] = multipletests(res_df["p_value"], method="fdr_bh")[1]

        res_df = res_df[
            (res_df["p_adj"] < 0.05) &
            (res_df["lift_combo_absent"] >= LIFT_THRESHOLD) &
            (res_df["lift_combo_present"] <= INV_LIFT_THRESHOLD)
        ].sort_values("lift_combo_absent", ascending=False)

        res_df = res_df.sort_values("intron")

    print(res_df)

    return res_df


# https://uc-r.github.io/hc_clustering
def run_pairwise_clustering(df, feature_cols, min_support, distance_threshold=0.1, min_feature_freq=0.01, exclusivity_score_threshold=0.5, show_dendrogram=False, LIFT_THRESHOLD=1.5, FEATURE_DISTANCE_THRESHOLD=100):
    rows_of_tokens = []

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
        rows_of_tokens.append(tokens)

    # --- Drop rare features: keep only tokens present in >= min_feature_freq of reads ---
    min_token_rows = int(numpy.ceil(min_feature_freq * len(df)))
    token_row_counts = Counter()
    for tokens in rows_of_tokens:
        for t in set(tokens):  # count once per row, even if repeated within a row
            token_row_counts[t] += 1

    keep_tokens = {t for t, c in token_row_counts.items() if c >= min_token_rows}
    rows = [[t for t in tokens if t in keep_tokens] for tokens in rows_of_tokens]

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

    min_cluster_freq = 0.5      # token must be present in at least 50% of this cluster's reads
    min_enrichment = 2.0        # AND at least 2x more frequent in this cluster than elsewhere

    cluster_names = {}
    cluster_summary = []
    important_tokens = set()

    n_total = len(rows)
    global_counts = Counter()
    for tokens in rows:
        global_counts.update(set(tokens))

    for cluster_id in sorted(set(labels)):
        idx = [i for i, lbl in enumerate(labels) if lbl == cluster_id]
        n_cluster = len(idx)

        counts = Counter()
        for i in idx:
            counts.update(set(rows[i]))

        # --- filter to tokens that are actually characteristic of this cluster ---
        candidates = []
        for tok, c_in in counts.items():
            freq_in_cluster = c_in / n_cluster

            c_total = global_counts[tok]
            c_out = c_total - c_in
            n_out = n_total - n_cluster
            freq_outside = c_out / n_out if n_out > 0 else 0.0

            # avoid divide-by-zero when a token is unique to this cluster (freq_outside == 0)
            enrichment = freq_in_cluster / freq_outside if freq_outside > 0 else numpy.inf

            if freq_in_cluster >= min_cluster_freq and enrichment >= min_enrichment:
                candidates.append((tok, c_in, freq_in_cluster, enrichment))

        # sort defining tokens by within-cluster frequency (most consistent first)
        candidates.sort(key=lambda x: x[2], reverse=True)

        top_n = 10
        top = candidates[:top_n]
        top_str = ", ".join(f"{tok} ({c}/{n_cluster}, {enr:.1f}x enriched)" for tok, c, freq, enr in top)

        top_tokens_clean = [sanitize(tok) for tok, _, _, _ in top]
        cluster_name = "_".join([f"cluster{cluster_id}"] + top_tokens_clean) if top_tokens_clean else f"cluster{cluster_id}"
        cluster_names[cluster_id] = cluster_name

        important_tokens.update(tok for tok, _, _, _ in top)

        cluster_summary.append({
            "cluster": cluster_name,
            "n_samples": n_cluster,
            "top_positions": top_str,
        })


    cluster_important_tokens = {}   # <-- NEW: cluster_id -> list of ALL its defining tokens (not just top_n)

    for cluster_id in sorted(set(labels)):
        idx = [i for i, lbl in enumerate(labels) if lbl == cluster_id]
        n_cluster = len(idx)

        counts = Counter()
        for i in idx:
            counts.update(set(rows[i]))

        candidates = []
        for tok, c_in in counts.items():
            freq_in_cluster = c_in / n_cluster
            c_total = global_counts[tok]
            c_out = c_total - c_in
            n_out = n_total - n_cluster
            freq_outside = c_out / n_out if n_out > 0 else 0.0
            enrichment = freq_in_cluster / freq_outside if freq_outside > 0 else numpy.inf

            if freq_in_cluster >= min_cluster_freq and enrichment >= min_enrichment:
                candidates.append((tok, c_in, freq_in_cluster, enrichment))

        candidates.sort(key=lambda x: x[2], reverse=True)
        cluster_important_tokens[cluster_id] = [tok for tok, _, _, _ in candidates]   # <-- store ALL, not just top_n


    df_out = df.copy()
    df_out["combined_primary_itemset"] = [cluster_names[int(lbl)] for lbl in labels]

    # --- optional: test for correlated / mutually exclusive feature pairs ---
    # Here is where we identify sets of features which are good predictors of each other, other co-occuring or mutually exclusive
    cluster_summary_df = pandas.DataFrame(cluster_summary).sort_values("n_samples", ascending=False).reset_index(drop=True)

    # for each cluster, grab the top tokens and test whether those tokens are better at predicting each intron than the background
    intron_predictor_results = evaluate_cluster_tokens_as_intron_predictors(
        rows, cluster_important_tokens, "introns", 5, LIFT_THRESHOLD, FEATURE_DISTANCE_THRESHOLD
    )

    if not intron_predictor_results.empty:
        collapsed_predictor_results = collapse_redundant_predictor_combos(
            intron_predictor_results, lift_col="lift_combo_absent", tolerance=0.05
        )
        # 2. Explode "clusters" so each row has ONE cluster_id
        collapsed_predictor_results["clusters"] = collapsed_predictor_results["clusters"].astype(str)
        predictor_exploded = collapsed_predictor_results.assign(
            cluster_id=collapsed_predictor_results["clusters"].str.split(",")
        ).explode("cluster_id")
        predictor_exploded["cluster_id"] = predictor_exploded["cluster_id"].str.strip().astype(int)
        predictor_exploded = predictor_exploded.drop(columns="clusters")

        # 3. Merge onto cluster_summary_df (this WILL produce multiple rows per cluster again — expected)
        cluster_summary_df["cluster_id"] = cluster_summary_df["cluster"].str.extract(r"^cluster(\d+)")[0].astype(int)
        merged_long = cluster_summary_df.merge(predictor_exploded, on="cluster_id", how="left")

        # 4. Collapse BACK to one row per cluster, aggregating predictor-result columns into lists
        list_cols = [c for c in predictor_exploded.columns if c != "cluster_id"]
        cluster_cols = [c for c in cluster_summary_df.columns if c != "cluster_id"]

        agg_dict = {c: (lambda s: list(s.dropna())) for c in list_cols}
        agg_dict.update({c: "first" for c in cluster_cols})

        cluster_summary_df = (
            merged_long
            .groupby("cluster_id", as_index=False)
            .agg(agg_dict)
        )



    # ------------------ plot dendrogram ---------------------
    if show_dendrogram:
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

            plt.xlabel("")
            plt.tight_layout()  # reserve right 15% of figure for the legend

            plt.show()

    return df_out, cluster_summary_df

def overlap_length(read, region_start, region_end):
    """Compute how many bp of a read's aligned blocks overlap [region_start, region_end)."""
    total = 0
    for block_start, block_end in read.get_blocks():  # aligned (ref-coordinate) blocks, excludes introns/deletions
        ov_start = max(block_start, region_start)
        ov_end = min(block_end, region_end)
        if ov_end > ov_start:
            total += ov_end - ov_start
    return total

def gather_read_entries_for_region(
    row,
    bam_labels,
    input_files,
    bam_handles,
    COVERAGE_PADDING,
    MIN_DELETION_LENGTH,
    MODS,
    PYSAM_MOD_THRESHOLD,
    MIN_GENE_OVERLAP
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

        # TODO, consider a % here
        if MIN_GENE_OVERLAP > 0:
            READS_IN_REGION = [
                read for read in READS_IN_REGION
                if overlap_length(read, row["start"], row["end"]) >= MIN_GENE_OVERLAP
            ]

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
    EXCLUSIVITY_THRESHOLD = args.exclusivity_threshold
    MIN_GENE_OVERLAP = args.min_gene_overlap
    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)
    LIFT_THRESHOLD = args.lift_threshold
    FEATURE_DISTANCE_THRESHOLD = args.feature_distance_threshold

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
    all_cluster_summaries = []

    try:
        # ------------------- MAIN READ AND FEATURE EXTRACTION LOOP FROM BAMFILES ------------------- #
        for _, row in matches.iterrows():
            if row["seq_id"] != "Pf3D7_12_v3":
                print("skipping {}...".format(row["ID"]))
                continue
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
                MIN_GENE_OVERLAP
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
                    # print("USING MIN_CLUSTER_PERC {}: min_cluster_size = {}".format(MIN_CLUSTER_PERC, min_cluster_size))
                else:
                    min_cluster_size = MIN_CLUSTER_SIZE_BULK
                    # print("USING MIN_CLUSTER_SIZE: min_cluster_size = {}".format(min_cluster_size))

                # So for a gene you can cluster by continuous variables (euclidean distance) or by categorical features (Jaccard distance) or by a combination of both (e.g. weighted sum of distances). The latter is experimental and may not work well, but it is possible to implement.
                df_clustered, cluster_summary = run_pairwise_clustering(
                    df, CLUSTER_COLS, min_cluster_size, DISTANCE_THRESHOLD, MIN_FEATURE_FREQ, EXCLUSIVITY_THRESHOLD,SHOW_DENDROGRAM, LIFT_THRESHOLD, FEATURE_DISTANCE_THRESHOLD
                )

                df_clustered["cluster"] = df_clustered["combined_primary_itemset"]

                # ---------- tag cluster_summary with the same row_key used for ct, and accumulate ----------
                cs = cluster_summary.copy()
                cs["row_key"] = str(row["ID"]) + "_" + cs["cluster"].astype(str)

                ct = (
                    df_clustered
                    .assign(row_key=lambda d: d["ID"].astype(str) + "_" + d["cluster"].astype(str))
                    .groupby(["row_key", "label"])
                    .size()
                    .unstack(fill_value=0)
                )

                # ---------- optional UMAP only for visualization ----------
                if len(matches) == 1:
                    UMAP_OUTPUT_FILE = "{}.umap".format(OUTPUT_FILE)
                    df_clustered.to_csv(UMAP_OUTPUT_FILE, sep="\t", index=False)
                    print("UMAP_OUTPUT_FILE")
                    print(f"Done. Wrote: {UMAP_OUTPUT_FILE}")

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

                cs = pandas.DataFrame([{
                    "cluster": "clusterNA",
                    "n_samples": 0 if (df is None or df.empty) else len(df),
                    "row_key": f"{row['ID']}_clusterNA",
                }])

            all_cluster_summaries.append(cs)
            all_count_tables.append(ct)
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

    # ---------- write cluster summaries, ordered the same as the count table ----------
    main_cluster_summary = (
        pandas.concat(all_cluster_summaries, axis=0)
        .set_index("row_key")
        .rename_axis("ID_cluster")
        .reindex(main_count_table.index)   # same order as main_count_table
        .reset_index()
    )

    CLUSTER_SUMMARY_OUTPUT_FILE = "{}.cluster_summary.tsv".format(OUTPUT_FILE)
    main_cluster_summary.to_csv(CLUSTER_SUMMARY_OUTPUT_FILE, sep="\t", index=False)
    print(f"Done. Wrote: {CLUSTER_SUMMARY_OUTPUT_FILE}")
