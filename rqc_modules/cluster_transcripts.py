import pandas
import pysam
import numpy
import umap

from sklearn.feature_extraction import DictVectorizer
from sklearn.preprocessing import normalize
from sklearn.cluster import HDBSCAN

import scipy.sparse as sp

from rqc_modules.constants import PYSAM_MOD_TUPLES
from rqc_modules.utils import process_input_files, process_annotation_file

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

def shorten_feature(f):
    # intron_409522_409679 -> I409522_409679
    # m6A_positions_12345 -> m6A12345

    if f.startswith("intron_"):
        return "I" + f.replace("intron_", "")

    if "_positions_" in f:
        mod, pos = f.split("_positions_")
        return mod + pos

    return f

def make_isoform_name(row, gene_id):
    """
    Generate an isoform-style cluster name from a binary feature vector.

    Example:
        PF3D7_0210000_I409522_409679_m6A410000_m6A410200
    """

    present = [
        shorten_feature(f)
        for f in row.index[row]
    ]

    if len(present) == 0:
        return f"{gene_id}_no_features"

    return f"{gene_id}_" + "_".join(present)


def assign_mod_cluster_ids(df, mod_cols):
    """
    Append a ``mod_cluster_id`` column to ``df`` based on frequent mod-site combinations.

    For each row, build a tuple of ``(mod_column, position)`` pairs from the selected
    mod columns. Count combinations across reads, rank them by frequency, and keep
    combinations that occur in more than 1% of reads. Retained combinations are
    named by joining their mod positions (e.g. "col1:100_col2:205"); rows with
    non-retained combinations are labeled "noise".
    """
    if not isinstance(mod_cols, (list, tuple)):
        mod_cols = [mod_cols]

    mod_cols = [col for col in mod_cols if col in df.columns]

    if not mod_cols:
        df["mod_cluster_id"] = pandas.Series("noise", index=df.index, dtype="object")
        return df

    def _build_combo(row):
        sites = []
        for col in mod_cols:
            vals = row[col]
            if isinstance(vals, (list, tuple, set)):
                for v in vals:
                    if pandas.notna(v):
                        sites.append((col, v))
            elif pandas.notna(vals):
                sites.append((col, vals))

        if not sites:
            return pandas.NA

        return tuple(sorted(sites))

    def _combo_to_name(combo):
        """Turn a tuple of (col, pos) pairs into a readable cluster name."""
        return "_".join(f"{col}:{pos}" for col, pos in combo)

    combo_series = df.apply(_build_combo, axis=1)
    counts = combo_series.value_counts(dropna=False)
    threshold = max(1, int(numpy.ceil(0.01 * len(df))))
    kept = counts[counts > threshold]

    if kept.empty:
        df["mod_cluster_id"] = pandas.Series("noise", index=df.index, dtype="object")
        return df

    combo_to_cluster_id = {
        combo: _combo_to_name(combo)
        for combo in kept.sort_values(ascending=False).index
        if pandas.notna(combo)
    }

    assigned_ids = [
        combo_to_cluster_id.get(combo, "noise")
        for combo in combo_series
    ]
    df["mod_cluster_id"] = pandas.Series(assigned_ids, index=df.index, dtype="object")
    return df


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


def compute_freq_weighted(X):
    """
    Equivalent to R compute_freq_weighted()
    X: sparse or dense count matrix
    """

    if sp.issparse(X):
        row_sums = numpy.asarray(X.sum(axis=1)).ravel()
        row_sums[row_sums == 0] = 1

        # row normalisation
        TF = X.multiply(1 / row_sums[:, None])

        # column frequency
        doc_freq = numpy.asarray((X > 0).sum(axis=0)).ravel()

    else:
        row_sums = X.sum(axis=1)
        row_sums[row_sums == 0] = 1

        TF = X / row_sums[:, None]

        doc_freq = (X > 0).sum(axis=0)

    freq_weight = numpy.log(1 + doc_freq)

    return TF.multiply(freq_weight) if sp.issparse(TF) else TF * freq_weight

def run_clustering(
    df,
    min_cluster_size,
    use_mods=True,
    mod_cols=("m6A_positions",),
    use_introns=False,
    intron_col="introns",
    min_feature_reads=None,          # e.g. 0.01*n_reads or min_cluster_size
    mod_weight=1.0,
    intron_weight=1.0
):
    blocks = []
    feat_names = []

    df = assign_mod_cluster_ids(df, mod_cols)
    print(df)

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

    X = sp.hstack(blocks, format="csr").astype(numpy.float32)
    feature_names = numpy.concatenate(feat_names) if feat_names else numpy.array([], dtype=object)

    # feature frequency filter
    min_feature_reads = int(numpy.ceil(0.01 * len(df)))
    feat_presence = numpy.asarray((X > 0).sum(axis=0)).ravel()

    keep = feat_presence >= int(min_feature_reads)
    X = X[:, keep]
    feature_names = feature_names[keep]
    feat_presence = feat_presence[keep]

    if X.shape[1] == 0:
        raise ValueError(f"No features remain after filtering with min_feature_reads={min_feature_reads}")

    X = compute_freq_weighted(X)

    X_dense = X.toarray().astype(numpy.float32)
    X_norm = normalize(X_dense, norm="l2")

    # new idea: take the most common features (>1% of reads)
    # bin reads by each permutation
    # assign a mod binary feature cluster id as a new continuous variable
    # use HDBSCAN on all continuous variables

    clusterer = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=20,
        metric="euclidean",
        allow_single_cluster=True
    )
    labels = clusterer.fit_predict(X_norm)

    cluster_names = {}

    for c in sorted(numpy.unique(labels)):
        if c == -1:
            cluster_names[c] = f"_noise"
            continue

        # Use original binary feature matrix
        cluster_X = X.tocsr()[labels == c]
        prevalence = numpy.asarray(cluster_X.mean(axis=0)).ravel()
        # Include features seen in at least 10% of cluster reads
        present = feature_names[prevalence >= 0.10]
        short = [shorten_feature(f) for f in present]

        if len(short) == 0:
            # fallback: show the most common features even if rare
            top = numpy.argsort(prevalence)[::-1][:5]

            short = [
                shorten_feature(feature_names[i])
                for i in top
                if prevalence[i] > 0
            ]

        if len(short) == 0:
            cluster_names[c] = f"_empty"
        else:
            cluster_names[c] = f"_".join(short)

    df["cluster"] = [cluster_names[c] for c in labels]

    out = df.copy()
    return out, X_norm, feature_names

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

                df_clustered, X_final, feat_names = run_clustering(
                    df=df,
                    min_cluster_size=min_cluster_size,
                    use_mods=CLUSTER_MODS,
                    mod_cols=mod_col_names,  # pick any mod columns
                    use_introns=CLUSTER_INTRONS,
                    intron_col="introns",
                    min_feature_reads=min_cluster_size,            # or int(numpy.ceil(0.01 * len(df)))
                    mod_weight=1.0,
                    intron_weight=1.0
                )
        
                # ---------- optional UMAP only for visualization ----------
                if len(matches) == 1:
                    print("generating UMAP embedding (viz only)...")
                    emb = umap.UMAP(
                        n_neighbors=20,
                        min_dist=0.3,
                        metric="euclidean"
                    ).fit_transform(X_final)
            
                    df_clustered["umap_x"] = emb[:, 0]
                    df_clustered["umap_y"] = emb[:, 1]
                    print("DONE")
        
                    # ---------- write ----------
                    UMAP_OUTPUT_FILE = "{}.umap".format(OUTPUT_FILE)
                    df_clustered.to_csv(UMAP_OUTPUT_FILE, sep="\t", index=False)
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
