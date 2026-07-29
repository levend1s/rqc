import pandas
import pysam
import numpy
import matplotlib.pyplot as plt
import umap
import hdbscan
import re

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances
from sklearn.feature_extraction import DictVectorizer
import scipy.sparse as sp
from sklearn.preprocessing import normalize
from sklearn.cluster import MiniBatchKMeans
from sklearn.mixture import GaussianMixture
from sklearn.feature_extraction.text import TfidfTransformer

from rqc_modules.constants import PYSAM_MOD_TUPLES
from rqc_modules.utils import process_input_files, process_annotation_file


import numpy as np
import pandas as pd
from scipy import sparse as sp

def explain_hdbscan_clusters(
    X,                      # matrix used for clustering (sparse or dense), shape (n_reads, n_features)
    labels,                 # cluster labels from HDBSCAN, shape (n_reads,)
    feature_names,          # feature names, len = n_features
    probabilities=None,     # optional: clusterer.probabilities_
    outlier_scores=None,    # optional: clusterer.outlier_scores_
    top_n=10,
    min_prevalence_in_cluster=0.10,
    include_noise=False
):
    """
    Returns:
      feature_table: DataFrame with top driving features per cluster
      cluster_table: DataFrame with per-cluster summary stats
      text_summary: dict cluster_id -> explanation string
    """

    labels = np.asarray(labels).astype(int)
    feature_names = np.asarray(feature_names, dtype=object)

    if sp.issparse(X):
        Xb = (X > 0).astype(np.int8)  # presence/absence
    else:
        Xb = (np.asarray(X) > 0).astype(np.int8)

    # Select rows to explain
    if include_noise:
        row_mask = np.ones(labels.shape[0], dtype=bool)
    else:
        row_mask = labels != -1

    Xb_use = Xb[row_mask]
    y = labels[row_mask]

    if y.size == 0:
        raise ValueError("No rows to explain (all points are noise).")

    # Cluster-level summary
    clusters = sorted(np.unique(y))
    cluster_rows = []
    for c in clusters:
        m = y == c
        n = int(m.sum())

        row = {
            "cluster": int(c),
            "n_reads": n,
            "fraction_of_explained_reads": n / len(y)
        }

        # Optional confidence/outlier summaries from full arrays
        if probabilities is not None:
            p = np.asarray(probabilities)[row_mask][m]
            row["mean_membership_probability"] = float(np.mean(p))
            row["median_membership_probability"] = float(np.median(p))
        if outlier_scores is not None:
            o = np.asarray(outlier_scores)[row_mask][m]
            row["mean_outlier_score"] = float(np.mean(o))
            row["median_outlier_score"] = float(np.median(o))

        cluster_rows.append(row)

    cluster_table = pd.DataFrame(cluster_rows).sort_values("cluster")

    # Feature enrichment per cluster vs rest
    feat_rows = []
    text_summary = {}

    # Convert mean(axis=0) safely for sparse/dense
    def col_mean(mat):
        if sp.issparse(mat):
            return np.asarray(mat.mean(axis=0)).ravel()
        return np.asarray(mat.mean(axis=0)).ravel()

    for c in clusters:
        in_c = (y == c)
        out_c = ~in_c

        # Need both sides to compare
        if in_c.sum() == 0 or out_c.sum() == 0:
            text_summary[int(c)] = f"Cluster {c}: cannot compare against rest."
            continue

        p_in = col_mean(Xb_use[in_c])     # prevalence in cluster
        p_out = col_mean(Xb_use[out_c])   # prevalence outside cluster

        diff = p_in - p_out
        l2fc = np.log2((p_in + 1e-6) / (p_out + 1e-6))

        # Keep enriched, sufficiently prevalent features
        candidate = np.where((diff > 0) & (p_in >= min_prevalence_in_cluster))[0]

        if candidate.size == 0:
            # show "best available" weak features instead of generic message
            # rank by prevalence diff (then log2FC), even if below thresholds
            top_weak_idx = np.lexsort((-l2fc, -diff))[::-1][:top_n]

            weak_parts = []
            for j in top_weak_idx:
                # only report features that are at least somewhat enriched
                if diff[j] <= 0:
                    continue
                weak_parts.append(
                    f"{feature_names[j]} ({p_in[j]*100:.1f}% in vs {p_out[j]*100:.1f}% out, "
                    f"Δ {diff[j]*100:.1f} pp, log2FC {l2fc[j]:.2f})"
                )

            cluster_n = int(in_c.sum())
            cluster_frac = 100.0 * cluster_n / len(y)

            if weak_parts:
                text_summary[int(c)] = (
                    f"Cluster {c} (n={cluster_n}, {cluster_frac:.1f}%): "
                    f"no strong single marker passed thresholds "
                    f"(min_prevalence_in_cluster={min_prevalence_in_cluster:.2f}). "
                    f"Top weakly enriched features: " + " ; ".join(weak_parts)
                )
            else:
                text_summary[int(c)] = (
                    f"Cluster {c} (n={cluster_n}, {cluster_frac:.1f}%): "
                    f"no enriched features over background (cluster likely defined by multifeature pattern or low signal)."
                )
            continue

        # Rank by prevalence gain, then log2 fold-change
        order = np.lexsort((-l2fc[candidate], -diff[candidate]))[::-1]
        top_idx = candidate[order[:top_n]]

        parts = []
        for j in top_idx:
            feat_rows.append({
                "cluster": int(c),
                "feature": str(feature_names[j]),
                "prevalence_in_cluster": float(p_in[j]),
                "prevalence_outside_cluster": float(p_out[j]),
                "prevalence_diff": float(diff[j]),
                "log2_fc": float(l2fc[j]),
            })
            parts.append(
                f"{feature_names[j]} ({p_in[j]*100:.1f}% in vs {p_out[j]*100:.1f}% out, "
                f"Δ {diff[j]*100:.1f} pp)"
            )

        text_summary[int(c)] = " ; ".join(parts)

    feature_table = pd.DataFrame(feat_rows)
    if not feature_table.empty:
        feature_table = feature_table.sort_values(
            ["cluster", "prevalence_diff", "log2_fc"],
            ascending=[True, False, False]
        )

    return feature_table, cluster_table, text_summary

import numpy as np
from sklearn.metrics import pairwise_distances

def diagnose_cluster_geometry(X, labels, clusterer, df, feature_names, c, top_n=15):
    y = np.asarray(labels).astype(int)
    in_c = (y == c)
    out_c = ~in_c
    if in_c.sum() == 0 or out_c.sum() == 0:
        print(f"Cluster {c}: empty in/out.")
        return

    # 1) core strength
    probs = getattr(clusterer, "probabilities_", None)
    if probs is not None:
        p_in = probs[in_c]
        print(f"Cluster {c}: n={in_c.sum()} mean_prob={p_in.mean():.3f} median_prob={np.median(p_in):.3f}")

    # 2) medoid of cluster
    Xi = X[in_c]
    D_in = pairwise_distances(Xi, metric="cosine")
    medoid_local = np.argmin(D_in.mean(axis=1))
    medoid_global = np.where(in_c)[0][medoid_local]

    # 3) nearest external reads to medoid
    D_out = pairwise_distances(X[medoid_global], X[out_c], metric="cosine").ravel()
    out_idx = np.where(out_c)[0]
    nn_out = out_idx[np.argsort(D_out)[:20]]

    # 4) compare centroid vs nearest-out centroid in feature space
    ci = np.asarray((X[in_c] > 0).mean(axis=0)).ravel()
    co = np.asarray((X[nn_out] > 0).mean(axis=0)).ravel()
    diff = ci - co

    top_pos = np.argsort(-diff)[:top_n]   # more present in cluster
    top_neg = np.argsort(diff)[:top_n]    # less present in cluster

    print("\nTop +features vs nearest outside:")
    for j in top_pos:
        if diff[j] <= 0: break
        print(f"  {feature_names[j]}: in={ci[j]:.3f}, near_out={co[j]:.3f}, Δ={diff[j]:.3f}")

    print("\nTop -features vs nearest outside:")
    for j in top_neg:
        if diff[j] >= 0: break
        print(f"  {feature_names[j]}: in={ci[j]:.3f}, near_out={co[j]:.3f}, Δ={diff[j]:.3f}")

    # 5) sample composition + technical covariates
    print("\nSample composition in cluster:")
    print(df.loc[in_c, "label"].value_counts(normalize=True))
    for col in ["read_length","average_quality","poly_a_length"]:
        if col in df.columns:
            print(f"{col}: in={df.loc[in_c,col].median():.2f}, out={df.loc[out_c,col].median():.2f}")

# scan region and pileup mods (list of genomic positions and their mod / unmod ratios)
# create a table where read_ids are row, columns are mods (m6A, m5C, pseU, m6A_inosine) with a list of mod positions (genomic space)
# other columns include read_start, read_end, read_length, read_strand, poly_A length, average_read_quality (could we cluster by individual base quality?)
# Perform dimensionality reduction (PCA, tSNE, UMAP) on this table and cluster reads based on their mod positions and other features
# create cartoon representation of read type that each cluster represents (e.g. m6A at position 100, m5C at position 150, etc.)


# Ok I'm considering how this will scale up. This could one day completely ignore annotations.
# 
# But I think the first solution
# - will loop through all gene annotations
# - generate tsv
# - umap
# - cluster (HBDBSCAN)
# - record gene_id_cluster_id and each read id associated with it
# output will be a featureCounts like file which can be passed to edgeR or salmon or w/e

# - If analysing a single gene:
# - segregate BAM files?
# - IGV screenshots

def build_mod_matrix(df, mod_cols):
    """Return sparse matrix from selected mod columns (list-like positions per row)."""
    mats = []
    names = []

    for col in mod_cols:
        row_dicts = []
        for vals in df[col]:
            d = {}
            if isinstance(vals, list):
                for v in vals:
                    k = f"{col}_{v}"          # keeps feature namespace per mod type
                    d[k] = d.get(k, 0) + 1    # counts per read
            row_dicts.append(d)

        vec = DictVectorizer(sparse=True)
        Xm = vec.fit_transform(row_dicts).tocsr()
        mats.append(Xm)
        names.extend(vec.get_feature_names_out().tolist())

    if not mats:
        return sp.csr_matrix((len(df), 0), dtype=numpy.float32), numpy.array([], dtype=object)

    return sp.hstack(mats, format="csr"), numpy.array(names, dtype=object)


def build_introns_matrix(df, intron_col="introns_binary"):
    """
    introns_binary can be:
    - scipy sparse matrix already
    - list/set of active intron ids per read (binary presence)
    """
    row_dicts = []
    for vals in df["introns"]:
        d = {}
        if isinstance(vals, list):
            for intr in vals:
                if isinstance(intr, tuple) and len(intr) == 2:
                    s, e = intr
                    k = f"intron_{s}_{e}"
                    d[k] = 1   # binary presence per read
        row_dicts.append(d)

    vec = DictVectorizer(sparse=True)
    Xi = vec.fit_transform(row_dicts).tocsr()
    return Xi, numpy.array(vec.get_feature_names_out(), dtype=object)


def run_clustering(
    df,
    min_cluster_size,
    use_mods=True,
    mod_cols=("m6A_positions",),
    use_introns=False,
    intron_col="introns",
    min_feature_reads=None,          # e.g. 0.01*n_reads or min_cluster_size
    weight_common=True,
    mod_weight=1.0,
    intron_weight=1.0
):
    blocks = []
    feat_names = []

    if use_mods:
        Xm, fnm = build_mod_matrix(df, mod_cols)
        if Xm.shape[1] > 0:
            blocks.append(Xm * numpy.sqrt(mod_weight))
            feat_names.append(fnm)

    if use_introns:
        Xi, fni = build_introns_matrix(df, intron_col=intron_col)
        if Xi.shape[1] > 0:
            blocks.append(Xi * numpy.sqrt(intron_weight))
            feat_names.append(fni)

    if not blocks:
        raise ValueError("No feature blocks selected or all selected blocks are empty.")

    X = sp.hstack(blocks, format="csr").astype(numpy.float64)
    feature_names = numpy.concatenate(feat_names) if feat_names else numpy.array([], dtype=object)

    # feature frequency filter
    if min_feature_reads is None:
        min_feature_reads = min_cluster_size  # or int(numpy.ceil(0.01 * n_reads))
    feat_presence = numpy.asarray((X > 0).sum(axis=0)).ravel()
    keep = feat_presence >= int(min_feature_reads)
    X = X[:, keep]
    feat_presence = feat_presence[keep]
    feature_names = feature_names[keep]

    if X.shape[1] == 0:
        raise ValueError(f"No features remain after filtering with min_feature_reads={min_feature_reads}")

    # TF-IDF
    tfidf = TfidfTransformer(
        norm="l2",          # usually best for cosine-based clustering
        use_idf=False,
        smooth_idf=True,
        sublinear_tf=False  # set True if you want 1+log(tf)
    )
    X = tfidf.fit_transform(X).astype(np.float64)


    if X.shape[1] == 1:
        x = X.toarray().ravel() if hasattr(X, "toarray") else numpy.asarray(X).ravel()
        # Option A: GMM 2-way split
        gm = GaussianMixture(n_components=2, random_state=0)
        labels = gm.fit_predict(x.reshape(-1, 1))

        # optional: mark uncertain points as noise
        probs = gm.predict_proba(x.reshape(-1, 1)).max(axis=1)
        labels = numpy.where(probs < 0.6, -1, labels)

    else:
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=max(5, min_cluster_size // 2),
            metric="cosine",
            cluster_selection_method="eom",
            allow_single_cluster=True
        )
        labels = clusterer.fit_predict(X)

        feature_table, cluster_table, text_summary = explain_hdbscan_clusters(
            X=X,
            labels=labels,
            feature_names=feature_names,
            probabilities=getattr(clusterer, "probabilities_", None),
            outlier_scores=getattr(clusterer, "outlier_scores_", None),
            top_n=8,
            min_prevalence_in_cluster=0.10,
            include_noise=False
        )

        # FIXME PF3D7_0210000 is hallucinating clusters, which actually seem real, but how is this happening from just introns?
        print(feature_names)
        print(cluster_table)
        print(feature_table.head(50))
        for c in sorted(text_summary):
            diagnose_cluster_geometry(X, labels, clusterer, df.reset_index(drop=True), feature_names, c=c, top_n=20)

    out = df.copy()
    out["cluster"] = labels.astype(int).astype(str)
    return out, X, feature_names

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

    USE_MOD_INFO = True
    USE_INTRON_INFO = True
    USE_CONTINUOUS_INFO = False

    MOD_WEIGHT = 0.3
    INTRON_WEIGHT = 0.3
    CONTINUOUS_WEIGHT = 0.4

    NUM_SAMPLES = 4
    MIN_CLUSTER_SIZE = 10 * NUM_SAMPLES
    MAX_CLUSTER_SIZE = 200 * NUM_SAMPLES
    MIN_CLUSTER_PERC = 0.01

    MINIMUM_READS_TO_PROCESS = 40

    UMAP_MIN_DIST = 0.3
    RANDOM_SEED = 42
    numpy.random.seed(RANDOM_SEED)

    # ---------- weights ----------
    weights = {
        "mod": MOD_WEIGHT if USE_MOD_INFO else 0.0,
        "intron": INTRON_WEIGHT if USE_INTRON_INFO else 0.0,
        "cont": CONTINUOUS_WEIGHT if USE_CONTINUOUS_INFO else 0.0,
    }
    s = sum(weights.values())
    if s == 0:
        raise ValueError("At least one feature type must be enabled.")
    weights = {k: v / s for k, v in weights.items()}

    w_mod = weights["mod"]
    w_intron = weights["intron"]
    w_cont = weights["cont"]

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

    MODS = ['m6A', 'm5C', 'pseU', 'm6A_inosine']

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

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)
    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']

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
            read_table = pandas.DataFrame(columns=read_table_header)
            read_table_index = 0

            for label in bam_labels:
                samfile_path = input_files[label]['path']
                # print("PROCESSING BAM: {}".format(samfile_path))

                samfile = bam_handles[label]

                READS_IN_REGION = list(samfile.fetch(
                                    contig=row['seq_id'], 
                                    start=row['start']-COVERAGE_PADDING, 
                                    stop=row['end']+COVERAGE_PADDING
                                ))
                # filter reads
                for i in range(len(READS_IN_REGION)):
                    r = READS_IN_REGION[i]
                    if r.is_secondary or r.is_supplementary:
                        continue
                    if (row["strand"] == "+" and r.is_reverse) or (row["strand"] == "-" and r.is_forward):
                        continue

                    mod_positions = {}

                    for mod in MODS:
                        mod_positions["{}_positions".format(mod)] = []
                        if r.is_forward:
                            pysam_mod_tuple_code = '{}_for'.format(mod)
                        else:
                            pysam_mod_tuple_code = '{}_rev'.format(mod)

                        if pysam_mod_tuple_code in PYSAM_MOD_TUPLES:

                            ref_pos = r.get_reference_positions(full_length=True)

                            mb = r.modified_bases
                            mods_probs = mb.get(PYSAM_MOD_TUPLES[pysam_mod_tuple_code]) if mb else None

                            if mods_probs:
                                # keep only mod positions which are above mod prob threshold
                                read_mod_positions = [x[0] for x in mods_probs if x[1] >= PYSAM_MOD_THRESHOLD]
                                genomic_mod_positions = [ref_pos[mod] for mod in read_mod_positions if ref_pos[mod] is not None]

                                mod_positions["{}_positions".format(mod)] = genomic_mod_positions

                    # introns
                    introns = []
                    ref_pos = r.reference_start

                    for op, length in r.cigartuples:
                        if op == 2:  # D = 
                            if length >= MIN_DELETION_LENGTH:
                                introns.append((ref_pos, ref_pos + length))

                        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
                            ref_pos += length

                    # phred quality
                    avg_quality = (
                        sum(r.query_qualities) / len(r.query_qualities)
                        if r.query_qualities
                        else 0
                    )                

                    poly_a_length = 0

                    if r.has_tag('pt:i'):
                        poly_a_length = r.get_tag('pt:i')

                    read_strand = '+' if r.is_forward else '-'

                    read_entry = [
                        r.query_name,
                        label,
                        MOD_PROB_THRESHOLD,
                        samfile_path,
                        row['seq_id'],
                        row['ID'],
                        r.reference_start,
                        r.reference_end,
                        read_strand,
                        r.query_length,
                        poly_a_length,
                        avg_quality,
                        introns,
                    ]

                    for mod in MODS:
                        read_entry.append(mod_positions["{}_positions".format(mod)])

                    read_table.loc[read_table_index] = read_entry
                    read_table_index += 1

            # ------------------- CLUSTER PROCESSING (OPTIMIZED) ------------------- #
            read_table = read_table.drop(columns=["poly_a_length", "average_quality", "read_start", "read_end", "read_length", "read_strand"])

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

                
                df_clustered, X_final, feat_names = run_clustering(
                    df=df,
                    min_cluster_size=min_cluster_size,
                    use_mods=CLUSTER_MODS,
                    mod_cols=mod_col_names,  # pick any mod columns
                    use_introns=CLUSTER_INTRONS,
                    intron_col="introns",
                    min_feature_reads=min_cluster_size,            # or int(np.ceil(0.01 * len(df)))
                    weight_common=True,
                    mod_weight=1.0,
                    intron_weight=1.0
                )

                # print(feat_names)
        
                # ---------- optional UMAP only for visualization ----------
                if len(matches) == 1:
                    print("generating UMAP embedding (viz only)...")
                    n_neighbors = min(50, max(15, min_cluster_size // 2))
                    emb = umap.UMAP(
                        n_neighbors=n_neighbors,
                        min_dist=0.3,
                        metric="cosine",
                        random_state=42
                    ).fit_transform(X_final)
            
                    df_clustered["umap_x"] = emb[:, 0]
                    df_clustered["umap_y"] = emb[:, 1]
                    print("DONE")
        
                    # ---------- write ----------
                    UMAP_OUTPUT_FILE = "{}.umap".format(OUTPUT_FILE)
                    df_clustered.to_csv(UMAP_OUTPUT_FILE, sep="\t", index=False)
                    print(f"Done. Wrote: {UMAP_OUTPUT_FILE}")
                    # print(df_clustered[["read_id", "cluster"]].head())

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
                    ct = pd.DataFrame(index=[f"{row['ID']}_clusterNA"])
                else:
                    ct = (
                        df.groupby("label")
                        .size()
                        .to_frame(name="n")
                        .T
                    )
                    ct.index = [f"{row['ID']}_clusterNA"]
            all_count_tables.append(ct)

            # HACK remove
            if row['seq_id'] == "Pf3D7_03_v3":
                break

            
    finally:
        # Always close all handles
        for fh in bam_handles.values():
            fh.close()

    # after loop: main combined table
    main_count_table = (
        pd.concat(all_count_tables, axis=0)
        .fillna(0)
        .astype(int)
        .sort_index()
    )

    # ---------- write ----------
    out_df = main_count_table.rename_axis("ID_cluster").reset_index()
    out_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"Done. Wrote: {OUTPUT_FILE}")
