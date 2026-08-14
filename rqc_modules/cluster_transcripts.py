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

def combine_exclusivity_scores(exclusivity_df, target_prefixes, min_expected_both=1.0):
    """
    For each target token, combine exclusivity_score across all its
    predictors using a weighted average — weighted by STATISTICAL
    CONFIDENCE (-log10(p_value)) rather than expected_both, so a common
    but weakly-associated predictor doesn't dominate the average just
    because it has high expected co-occurrence.

    Also reports the single strongest individual predictor per side, so
    a diluted combined average doesn't hide a genuinely strong signal.
    """
    empty_columns = [
        "target",
        "n_exclusive_predictors", "exclusive_predictors", "combined_exclusivity_score",
        "top_exclusive_predictor", "top_exclusive_score",
        "n_cooccurring_predictors", "cooccurring_predictors", "combined_cooccurrence_score",
        "top_cooccurring_predictor", "top_cooccurring_score",
    ]

    if exclusivity_df is None or len(exclusivity_df) == 0:
        return pandas.DataFrame(columns=empty_columns)

    target_to_pairs = {}
    for _, r in exclusivity_df.iterrows():
        a, b = r["token_a"], r["token_b"]
        score = r["exclusivity_score"]
        expected_both = r["expected_both"]
        p_value = r["p_value"]

        if expected_both < min_expected_both or pandas.isna(score):
            continue

        a_is_target = a.startswith(target_prefixes)
        b_is_target = b.startswith(target_prefixes)

        if a_is_target and not b_is_target:
            target_tok, predictor_tok = a, b
        elif b_is_target and not a_is_target:
            target_tok, predictor_tok = b, a
        else:
            continue

        # confidence weight: -log10(p), clipped so p=0 doesn't blow up to inf
        p_clipped = max(p_value, 1e-300)
        weight = -numpy.log10(p_clipped)

        target_to_pairs.setdefault(target_tok, []).append((predictor_tok, score, weight))

    results = []
    for target_tok, entries in target_to_pairs.items():
        exclusive_entries = [e for e in entries if e[1] >= 0]
        cooccurring_entries = [e for e in entries if e[1] < 0]

        def weighted_avg(subset):
            if not subset:
                return numpy.nan
            w_sum = sum(e[2] for e in subset)
            if w_sum == 0:
                return numpy.nan
            return sum(s * w for _, s, w in subset) / w_sum

        def top_entry(subset, most_extreme_fn):
            if not subset:
                return None, numpy.nan
            best = most_extreme_fn(subset, key=lambda e: e[1])
            return best[0], best[1]

        excl_score = weighted_avg(exclusive_entries)
        coocc_score = weighted_avg(cooccurring_entries)

        # print(f"Target {target_tok}: {len(exclusive_entries)} exclusive, {len(cooccurring_entries)} co-occurring predictors. Combined exclusivity score: {excl_score:.3f}, combined co-occurrence score: {coocc_score:.3f}")
        # print(f"Exclusive entries: {exclusive_entries}")
        # print(f"Co-occurring entries: {cooccurring_entries}")

        top_excl_pred, top_excl_score = top_entry(exclusive_entries, max)          # highest (most exclusive)
        top_coocc_pred, top_coocc_score = top_entry(cooccurring_entries, min)      # most negative (most co-occurring)

        results.append({
            "target": target_tok,
            "n_exclusive_predictors": len(exclusive_entries),
            "exclusive_predictors": [e[0] for e in exclusive_entries],
            "combined_exclusivity_score": excl_score,
            "top_exclusive_predictor": top_excl_pred,
            "top_exclusive_score": top_excl_score,
            "n_cooccurring_predictors": len(cooccurring_entries),
            "cooccurring_predictors": [e[0] for e in cooccurring_entries],
            "combined_cooccurrence_score": coocc_score,
            "top_cooccurring_predictor": top_coocc_pred,
            "top_cooccurring_score": top_coocc_score,
        })

    return pandas.DataFrame(results, columns=empty_columns) if results else pandas.DataFrame(columns=empty_columns)

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


def combine_predictor_evidence(exclusivity_df, target_prefixes, method="fisher"):
    """
    For each target token, combine the p-values of ALL its associated
    predictors (regardless of individual significance) into a single
    combined p-value, using Fisher's method (or another scipy-supported
    method). This lets weak-but-consistent evidence across multiple
    predictors accumulate into stronger evidence for the target.

    CAVEAT: this assumes predictors are reasonably independent. If your
    predictors are correlated (e.g. nearby m6A sites that tend to be
    modified together), this will overstate significance. Consider
    checking predictor correlation before trusting results here.

    Parameters
    ----------
    exclusivity_df : pandas.DataFrame
        Output of find_mutually_exclusive_pairs — needs token_a, token_b,
        p_value (raw, not p_adj — combine raw p-values, then correct once
        at the end across targets).
    target_prefixes : str or tuple[str, ...]
        Prefix(es) identifying which token in each pair is the target.
    method : str
        Method passed to scipy.stats.combine_pvalues ("fisher", "pearson",
        "tippett", "stouffer", etc.).

    Returns
    -------
    pandas.DataFrame with columns:
        target, n_predictors, predictors, combined_statistic,
        combined_p_value, combined_p_adj
    One row per target token.
    """
    target_to_pvals = {}
    target_to_predictors = {}

    for _, r in exclusivity_df.iterrows():
        a, b, p = r["token_a"], r["token_b"], r["p_value"]

        a_is_target = a.startswith(target_prefixes)
        b_is_target = b.startswith(target_prefixes)

        if a_is_target and not b_is_target:
            target_tok, predictor_tok = a, b
        elif b_is_target and not a_is_target:
            target_tok, predictor_tok = b, a
        else:
            continue  # ambiguous — skip

        target_to_pvals.setdefault(target_tok, []).append(p)
        target_to_predictors.setdefault(target_tok, []).append(predictor_tok)

    results = []
    for target_tok, pvals in target_to_pvals.items():
        # combine_pvalues needs at least 1 p-value; guard tiny/degenerate p's
        pvals_clean = [max(p, 1e-300) for p in pvals]  # avoid log(0) issues in fisher's method
        stat, combined_p = combine_pvalues(pvals_clean, method=method)

        results.append({
            "target": target_tok,
            "n_predictors": len(pvals),
            "predictors": target_to_predictors[target_tok],
            "combined_statistic": stat,
            "combined_p_value": combined_p,
        })

    res_df = pandas.DataFrame(results)
    if len(res_df) > 0:
        res_df["combined_p_adj"] = multipletests(res_df["combined_p_value"], method="fdr_bh")[1]
        res_df = res_df.sort_values("combined_p_adj")

    return res_df

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

def annotate_clusters_with_combined_exclusivity(cluster_summary_df, labels, rows,
                                                   combined_scores_df,
                                                   exclusivity_score_threshold=0.7,
                                                   present_frac_threshold=0.5):
    """
    For each cluster (one row per cluster), find which target tokens it
    carries (>= present_frac_threshold of reads), and attach that target's
    pre-combined exclusivity/co-occurrence scores from combined_scores_df.

    Parameters
    ----------
    combined_scores_df : pandas.DataFrame
        Output of combine_exclusivity_scores — one row per target, with
        combined_exclusivity_score, combined_cooccurrence_score, and the
        associated predictor lists.
    exclusivity_score_threshold : float
        Only targets with combined_exclusivity_score >= this are matched.
        Set to None to skip this filter and match on target presence alone.
    present_frac_threshold : float
        Fraction of a cluster's reads that must carry a target token for
        it to count as "present" in that cluster.

    Returns
    -------
    pandas.DataFrame with columns:
        cluster, n_samples, target, combined_exclusivity_score,
        exclusive_predictors, combined_cooccurrence_score, cooccurring_predictors
    Each of target/combined_exclusivity_score/etc. holds a LIST — one
    entry per matched target token, aligned by position. Empty lists if
    no matches for that cluster.
    """

    required_cols = ["target", "combined_exclusivity_score", "exclusive_predictors", "combined_cooccurrence_score", "cooccurring_predictors"]
    if combined_scores_df is None or len(combined_scores_df) == 0 or "combined_exclusivity_score" not in combined_scores_df.columns:
        candidates = pandas.DataFrame(columns=required_cols)
    else:
        candidates = combined_scores_df
        if exclusivity_score_threshold is not None:
            candidates = candidates[candidates["combined_exclusivity_score"] >= exclusivity_score_threshold]

    candidates = combined_scores_df
    if exclusivity_score_threshold is not None:
        candidates = candidates[
            candidates["combined_exclusivity_score"] >= exclusivity_score_threshold
        ]

    n_samples_lookup = dict(zip(cluster_summary_df["cluster"], cluster_summary_df["n_samples"]))
    cluster_ids_sorted = sorted(set(labels))
    cluster_name_lookup = {cid: name for cid, name in zip(cluster_ids_sorted, cluster_summary_df["cluster"])}

    per_cluster_tokens = {}
    for cid in cluster_ids_sorted:
        idx = [i for i, lbl in enumerate(labels) if lbl == cid]
        counts = Counter()
        for i in idx:
            counts.update(set(rows[i]))
        n = len(idx)
        per_cluster_tokens[cid] = {tok: c / n for tok, c in counts.items()}

    result_rows = []
    for cid in cluster_ids_sorted:
        tok_freqs = per_cluster_tokens[cid]
        present_tokens = {t for t, f in tok_freqs.items() if f >= present_frac_threshold}
        cluster_name = cluster_name_lookup.get(cid, f"cluster_{cid}")
        n_samples = n_samples_lookup.get(cluster_name, None)

        target_list = []
        excl_score_list = []
        excl_predictors_list = []
        coocc_score_list = []
        coocc_predictors_list = []

        for _, r in candidates.iterrows():
            if r["target"] in present_tokens:
                target_list.append(r["target"])
                excl_score_list.append(r["combined_exclusivity_score"])
                excl_predictors_list.append(r["exclusive_predictors"])
                coocc_score_list.append(r["combined_cooccurrence_score"])
                coocc_predictors_list.append(r["cooccurring_predictors"])

        result_rows.append({
            "cluster": cluster_name,
            "n_samples": n_samples,
            "target": target_list,
            "combined_exclusivity_score": excl_score_list,
            "exclusive_predictors": excl_predictors_list,
            "combined_cooccurrence_score": coocc_score_list,
            "cooccurring_predictors": coocc_predictors_list,
        })

    return pandas.DataFrame(result_rows, columns=[
        "cluster", "n_samples", "target", "combined_exclusivity_score",
        "exclusive_predictors", "combined_cooccurrence_score", "cooccurring_predictors"
    ])

# https://uc-r.github.io/hc_clustering
def run_pairwise_clustering(df, feature_cols, min_support, distance_threshold=0.1, min_feature_freq=0.01, exclusivity_score_threshold=0.5, show_dendrogram=False):
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

    def sanitize(tok):
        # keep names filesystem/column-friendly: no spaces, slashes, etc.
        return re.sub(r"[^0-9a-zA-Z_]+", "-", tok)

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

    df_out = df.copy()
    df_out["combined_primary_itemset"] = [cluster_names[int(lbl)] for lbl in labels]

    # --- optional: test for correlated / mutually exclusive feature pairs ---
    exclusivity_df = find_mutually_exclusive_pairs(
        rows,
        allowed_tokens=important_tokens,
        target_prefixes="introns",
        predictor_prefixes="m6A",
        min_individual_count=5,
    )

    combined_scores = combine_exclusivity_scores(exclusivity_df, target_prefixes="introns")
    # print(combined_scores)
    # print(exclusivity_df)
    # print(important_tokens)
    cluster_summary_df = pandas.DataFrame(cluster_summary)
    cluster_summary_df = annotate_clusters_with_combined_exclusivity(
        cluster_summary_df,
        labels=labels,
        rows=rows,
        combined_scores_df=combined_scores,
        exclusivity_score_threshold=exclusivity_score_threshold,
        present_frac_threshold=0.5,
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
    EXCLUSIVITY_THRESHOLD = args.exclusivity_threshold
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
                    df, CLUSTER_COLS, min_cluster_size, DISTANCE_THRESHOLD, MIN_FEATURE_FREQ, EXCLUSIVITY_THRESHOLD,SHOW_DENDROGRAM
                )

                df_clustered["cluster"] = df_clustered["combined_primary_itemset"]
                # print(cluster_summary)

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
