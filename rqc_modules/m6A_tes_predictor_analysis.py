import matplotlib.pyplot as plt
import matplotlib.colors
import scipy.stats
import os
import numpy
import pandas

from rqc_modules.utils import process_input_files, process_annotation_file, filter_gff_for_target_features, find_canonical_mods, get_filtered_reads_ids


def m6A_tes_predictor_analysis(args):
    ANNOTATION_FILE = args.annotation
    MOD_PROP_THRESHOLD = args.mod_prob_threshold
    MOD_RATIO = args.mod_ratio
    COVERAGE_PADDING = args.coverage_padding
    INPUT = args.input
    IDS = args.ids
    OUTPUT = args.output
    READ_DEPTH_THRESHOLD = args.read_depth
    OFFSET_PADDING = args.offset_padding
    TES_DISTANCE_THRESHOLD = args.tes_distance_threshold
    FEATURE_TYPE = args.type
    EXCLUDE_CONTIGS = args.exclude_contigs


    print("m6A_tes_predictor_analysis...")

    
    if OUTPUT:
        output_dir = os.path.dirname(OUTPUT)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        OUTPUT_COMPARE = OUTPUT.replace(".tsv", "_compare.tsv")
        OUTPUT_RAW = OUTPUT.replace(".tsv", "_raw.tsv")

    # process input file. Each line contains a label, the type of file, and the filepath
    input_files = process_input_files(INPUT)

    # load annotation file and find indexes for all parent children
    gff_df = process_annotation_file(ANNOTATION_FILE)

    if FEATURE_TYPE:
        matches = gff_df[gff_df['type'] == FEATURE_TYPE]

        if EXCLUDE_CONTIGS:
            matches = matches[~matches['seq_id'].isin(EXCLUDE_CONTIGS)]
    else:
        matches = gff_df[gff_df['ID'].isin(IDS)]

    if matches.empty:
        print("ERROR: no matches found for type {} and ids {}".format(FEATURE_TYPE, IDS))

    if matches.empty:
        print("ERROR: no matches found for ids {}".format(IDS))

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']

    canonical_mods = find_canonical_mods(
        matches, 
        input_files, 
        bam_labels, 
        MOD_RATIO,
        READ_DEPTH_THRESHOLD, 
        COVERAGE_PADDING,
        return_percent=True
    )
    # approximated_tes = approximate_tes(target_features, input_files, bam_labels)
    approximated_tes = None
    print(canonical_mods)

    canonical_mod_start_positions = {
        key: list(value.keys())
        for key, value in canonical_mods.items()
    }

    print(canonical_mod_start_positions)

    filtered_read_ids = get_filtered_reads_ids(
        matches, 
        input_files,
        bam_labels, 
        MOD_PROP_THRESHOLD,
        READ_DEPTH_THRESHOLD,
        COVERAGE_PADDING,
        canonical_mod_start_positions
        # approximated_tes = approximated_tes
    )

    # HACK manually remove total coverage for all mods
    for _, row in matches.iterrows():
        for label in bam_labels:
            filtered_read_ids[label][row.ID]['tx_end_sites'].pop('None', None)

    d_mod_tes_probs = {}

    # TODO: make more readable, refactor
    for _, row in matches.iterrows():
        # create hist / kde for all tes, and separate ones for reads containing each mod
        d_mod_tes_probs[row.ID] = {}

        for label in bam_labels:
            d_mod_tes_probs[row.ID][label] = {}
            for key, tx_end_sites in filtered_read_ids[label][row.ID]['tx_end_sites'].items():

                if key != 'all':
                    if row.strand == "-":
                        m6A_specific_tes_threshold = key - TES_DISTANCE_THRESHOLD
                    else:
                        m6A_specific_tes_threshold = key + TES_DISTANCE_THRESHOLD

                    tx_within = [x for x in tx_end_sites if (x >= m6A_specific_tes_threshold if row.strand == "-" else x <= m6A_specific_tes_threshold)]
                    num_tx_within = len(tx_within)

                    # remove those sites with tx_end_sites but the caonical m6A is present,
                    # this likely means the transcripts related to that mod are related to a different gene 

                    # HACK
                    if len(tx_end_sites) > READ_DEPTH_THRESHOLD:
                        prob = num_tx_within / len(tx_end_sites)
                        # if len(tx_end_sites) == 0 and canonical_mods[row.ID][key] > 0:
                        #     continue
                        d_mod_tes_probs[row.ID][label][key] = {
                            'probability': prob,
                            'percent_mod': canonical_mods[row.ID][key] / 100,
                            'num_tx_within': num_tx_within,
                            'total_tx': len(tx_end_sites)
                        }

    # ID, mod start, percent_mod, probability, num tx within total tx
    rows = []

    for gene_id, samples in d_mod_tes_probs.items():
        for sample, positions in samples.items():
            for mod_start, metrics in positions.items():
                rows.append({
                    "ID": gene_id,
                    "mod_start": mod_start,
                    "percent_mod": metrics["percent_mod"],
                    "probability": metrics["probability"],
                    "num_tx_within": metrics["num_tx_within"],
                    "total_tx": metrics["total_tx"],
                })

    summary_df = pandas.DataFrame(rows)

    if OUTPUT:
        summary_df.to_csv(OUTPUT, sep='\t', index=False)
