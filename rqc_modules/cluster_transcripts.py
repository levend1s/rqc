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


def cluster_transcripts(args):
    INPUT = args.input
    ANNOTATION_FILE = args.annotation
    IDS = args.ids
    COVERAGE_PADDING = args.padding
    MOD_PROB_THRESHOLD = args.mod_prob_threshold
    OUTPUT_FILE = args.outfile
    MIN_DELETION_LENGTH = args.min_deletion_length

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)

    ALL_CONTINUOUS_COLS = ["poly_a_length"]

    USE_MOD_INFO = True
    USE_INTRON_INFO = True
    USE_CONTINUOUS_INFO = False

    MOD_TYPES = ["m6A"]  # e.g. ["m6A", "m5C", "pseU", "m6A_inosine"]
    if not USE_MOD_INFO:
        MOD_TYPES = []

    MOD_WEIGHT = 0.3
    INTRON_WEIGHT = 0.3
    CONTINUOUS_WEIGHT = 0.4

    NUM_SAMPLES = 4
    MIN_CLUSTER_SIZE = 10 * NUM_SAMPLES
    MAX_CLUSTER_SIZE = 200 * NUM_SAMPLES
    MIN_CLUSTER_PERC = 0.01

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

    matches = gff_df[gff_df['ID'].isin(IDS)]

    if matches.empty:
        print("ERROR: no matches found for ids {}".format(IDS))

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

    read_table = pandas.DataFrame(columns=read_table_header)
    read_table_index = 0

    input_files = process_input_files(INPUT)

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)
    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']

    # Open once: label -> handle
    bam_handles = {
        label: pysam.AlignmentFile(input_files[label]['path'], "rb")
        for label in bam_labels
    }

    try:
        # ------------------- MAIN READ AND FEATURE EXTRACTION LOOP FROM BAMFILES ------------------- #
        for _, row in matches.iterrows():
            print("processing {}...".format(row["ID"]))
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
                            mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES[pysam_mod_tuple_code])

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
            df = read_table.reset_index(drop=True).copy()
            num_total_reads = len(df)
    
            min_cluster_size = min(
                MAX_CLUSTER_SIZE,
                max(MIN_CLUSTER_SIZE, int(numpy.ceil(MIN_CLUSTER_PERC * num_total_reads / NUM_SAMPLES)))
            )
    
            # ---------- sparse feature matrix ----------
            print("generating mod sparse matrix...")
            row_dicts = []
            for vals in df["m6A_positions"]:
                d = {}
                if isinstance(vals, list):
                    for v in vals:
                        k = f"m6A_{v}"
                        d[k] = d.get(k, 0) + 1
                row_dicts.append(d)
    
            X = DictVectorizer(sparse=True).fit_transform(row_dicts)  # sparse count matrix
            # 1) Drop features present in <1% of reads
            feat_presence = numpy.asarray((X > 0).sum(axis=0)).ravel()   # per-feature document frequency
            keep = feat_presence >= min_cluster_size
            X = X[:, keep]
            feat_presence = feat_presence[keep]
            # 2) Row-normalize counts (TF)
            X = normalize(X, norm="l1", axis=1)
            # 3) Up-weight common features: weight_j = log(1 + df_j)
            weights = numpy.log1p(feat_presence).astype(numpy.float32)      # bigger df => bigger weight
            X = X @ sp.diags(weights, format="csr")
            print(X)
            print("DONE")
    
            # ---------- CLUSTER ON FEATURE SPACE (NOT UMAP) ----------
            print("clustering on feature matrix...")
            clusterer = hdbscan.HDBSCAN(
                min_cluster_size=min_cluster_size,
                min_samples=max(5, min_cluster_size // 2),
                metric="cosine"   # if using sparse and this errors, switch to metric="cosine" with dense/reduced data
            )
            clusters = clusterer.fit_predict(X)
            df["cluster"] = clusters.astype(int).astype(str)
            print("DONE")
    
            # ---------- optional UMAP only for visualization ----------
            if len(matches) == 1:
                print("generating UMAP embedding (viz only)...")
                n_neighbors = min(50, max(15, min_cluster_size // 2))
                emb = umap.UMAP(
                    n_neighbors=n_neighbors,
                    min_dist=0.3,
                    metric="cosine",
                    random_state=42
                ).fit_transform(X)
        
                df["umap_x"] = emb[:, 0]
                df["umap_y"] = emb[:, 1]
                print("DONE")
    
                # ---------- write ----------
                df.to_csv(OUTPUT_FILE, sep="\t", index=False)
                print(f"Done. Wrote: {OUTPUT_FILE}")
                print(df[["read_id", "cluster"]].head())        
    finally:
        # Always close all handles
        for fh in bam_handles.values():
            fh.close()
