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

def run_clustering(
    df,
    mod_cols=[],
):
    # TODO implement tolerance
    min_percent = 1.0
    n_reads = len(df)
    min_reads_cluster = (min_percent / 100.0) * n_reads  # strict ">" threshold
    _, frequent_features, _ = append_features_present_in_gt_percent(df, mod_cols, min_reads_cluster)
    print(frequent_features)
    print(n_reads)
    print(min_reads_cluster)
    genomic_feature_clusters = {}

    # build a dictionary of the most frequent unique clusters
    for _, row in df.iterrows():
        cluster_id = ""

        for col in mod_cols:
            row_features = set(row.get(col, []))
            frequent_features_col = frequent_features[col]
            intersect_features = sorted(row_features & frequent_features_col)
            something_to_add = "-".join([f"{col}_{f}" for f in intersect_features])
            if cluster_id and something_to_add:
                cluster_id += "-"
            cluster_id += something_to_add

        genomic_feature_clusters[cluster_id] = 1 if cluster_id not in genomic_feature_clusters else genomic_feature_clusters.get(cluster_id, 0) + 1

    # keep only clusters which explain >1% of reads
    filtered_genomic_feature_clusters = {k: v for k, v in genomic_feature_clusters.items() if v > min_reads_cluster}
    print(genomic_feature_clusters)
    print(filtered_genomic_feature_clusters)


    # TODO maybe remove
    # frequent_features_set = []
    # for k, v in frequent_features.items():
    #     frequent_features_set += [f"{k}_{e}" for e in v]

    # frequent_features_set = set(frequent_features_set)

    frequent_combination_features_set = []
    # convert back to sets for quick comparison
    for k, v in filtered_genomic_feature_clusters.items():
        features = [f for f in k.split("-") if f]
        frequent_combination_features_set.append( set(features))


    print(frequent_combination_features_set)

    cluster_ids = []
    # go through reads again, and assign them to one of the unique_features
    # if not assigned to a filtered cluster, try assign to frequent_features
    # otherwise, assign to noise group
    for _, row in df.iterrows():
        cluster_id = ""

        row_features = []

        for col in mod_cols:
            row_mod_features = row.get(col, [])
            row_features += [f"{col}_{f}" for f in row_mod_features]

        row_features = set(row_features)

        # find intersection in top combos
        if not row_features:
            matched_cluster = "empty"          # no features at all
        else:
            matched_cluster = "noise"     # has features, but no matching combo yet

            for combo in sorted(frequent_combination_features_set, key=len, reverse=True):
                if not isinstance(combo, set) or not combo:
                    continue

                intersect_features = row_features & combo
                if intersect_features:
                    matched_cluster = "-".join(sorted(intersect_features))
                    break

        cluster_ids.append(matched_cluster)

    df["cluster"] = cluster_ids
    return df

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

                df_clustered = run_clustering(
                    df=df,
                    mod_cols=mod_col_names,  # pick any mod columns
                )
        
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
