import pandas
import pysam
import numpy
import math
from collections import Counter

from rqc_modules.utils import getSubfeatures, process_input_files, process_annotation_file
from rqc_modules.constants import PYSAM_MOD_TUPLES, MODKIT_BEDMETHYL_HEADER, FEATURECOUNTS_HEADER, BAM_PILEUP_DEFAULT_FLAGS, PYSAM_PILEUP_MAX_DEPTH, BAM_REVERSE_STRAND

from rqc_modules.utils import START_CLOCK, STOP_CLOCK

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from functools import partial

import statsmodels.formula.api as smf
import statsmodels.api as sm

import multiprocessing
multiprocessing.set_start_method('fork')  # Use 'fork' instead of 'spawn'
numpy.seterr(divide='ignore', invalid='ignore')

# this is fine as long as reads_in_region is passed by reference cause its big
def process_row(row, label, samfile_path, coverage_padding, pysam_mod_threshold, read_depth_threshold, canonical_mod_prop_threshold):
    gene_length = row['end'] - row['start']

    samfile = pysam.AlignmentFile(samfile_path, 'rb')
    READS_IN_REGION = list(samfile.fetch(
        contig=row['seq_id'], 
        start=row['start'] - coverage_padding, 
        stop=row['end'] + coverage_padding
    ))

    read_indexes_to_process = []

    # generate coverage for this gene
    row_flag_filters = BAM_PILEUP_DEFAULT_FLAGS
    row_flags_requires = 0

    if row['strand'] == '+':
        row_flag_filters = BAM_PILEUP_DEFAULT_FLAGS | BAM_REVERSE_STRAND
        a_code = 'A'
        pysam_mod_tuple_code = 'm6A_for'
    else:
        row_flags_requires = BAM_REVERSE_STRAND
        a_code = 't'
        pysam_mod_tuple_code = 'm6A_rev'

    gene_length = row['end'] - row['start'] + 1
    d_coverage = {
        "total_depth": [0] * gene_length,
        "count_a": [0] * gene_length,
        "count_t": [0] * gene_length,
        "m6A": [0] * gene_length
    }

    for column in samfile.pileup(
        contig=row['seq_id'], 
        start=row['start'] - 1, 
        stop=row['end'],
        # min_mapping_quality=MIN_MAPQ,
        max_depth=PYSAM_PILEUP_MAX_DEPTH,
        flag_require=row_flags_requires,
        flag_filter=row_flag_filters,
        truncate = True,
        min_base_quality=0
    ):
        # take the number of aligned reads at this column position (column.n) minus the number of aligned reads which have either a skip or delete base at this column position (r.query_position)
        column_bases_read = list(column.get_query_sequences(add_indels=True))
        
        d_coverage['total_depth'][column.reference_pos - row['start']] = len(column_bases_read)
        d_coverage['count_a'][column.reference_pos - row['start']] = column_bases_read.count(a_code)

    max_read_depth = max(d_coverage['total_depth'])
    if max_read_depth < read_depth_threshold:
        row_summary = [label, row.ID, row.type, row.strand, max_read_depth, [], 0, 0, 0]
        print("\t".join(map(str, row_summary)))

        return row_summary


    read_indexes_to_process = []
    # filter reads
    for i in range(len(READS_IN_REGION)):
        r = READS_IN_REGION[i]

        # keep only reads in the same direction as this strand
        if (row['strand'] == "+" and r.is_forward) or (row['strand'] == "-" and r.is_reverse):
            if row['strand'] == "-":
                read_3p_end = r.reference_start

                if read_3p_end <= row.end:
                    read_indexes_to_process.append(i)
            else:
                read_3p_end = r.reference_end

                if read_3p_end >= row.start:
                    read_indexes_to_process.append(i)
                    

    filtered_reads_in_region = [READS_IN_REGION[i] for i in read_indexes_to_process]

    all_genomic_mod_positions = []

    for r in filtered_reads_in_region:
        ref_pos = r.get_reference_positions(full_length=True)
        mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES[pysam_mod_tuple_code])

        if mods_probs:
            # keep only mod positions which are above mod prob threshold
            read_mod_positions = [x[0] for x in mods_probs if x[1] >= pysam_mod_threshold]
            genomic_mod_positions = [ref_pos[mod] for mod in read_mod_positions if ref_pos[mod] is not None]

            all_genomic_mod_positions.append(genomic_mod_positions)

    flattened_genomic_mod_positions = [item for sublist in all_genomic_mod_positions for item in sublist]
    genomic_mod_counts = Counter(flattened_genomic_mod_positions)
    for pos, count in genomic_mod_counts.items():

        # this ignores mods in a read that are outside of the gene annotation
        if pos <= row.end and pos >= row.start:
            d_coverage['m6A'][pos - row['start']] = count

    # mod_ratio = numpy.array(d_coverage['m6A']) / numpy.array(d_coverage['total_depth'])

    mod_ratio = [(a / b) if a > 0 and b > 0 else 0 for a, b in zip(d_coverage['m6A'], d_coverage['total_depth'])]
    genomic_coords = [i + row['start'] for i in range(gene_length)]
    mod_ratio_read_depth_tuples = list(zip(mod_ratio, d_coverage['total_depth'], d_coverage['m6A'], genomic_coords))


    # now I have coverages for m6A, total depth and Adenosines per gene
    # Now i'd like to calculate average methylation for the gene. I want to do this for only "canonical" m6A positions (defined by a m6A/A ratio), for non-canonical m6A positions, and a total m6A/A ratio.
    # first, create a new list of m6A/A ratios for each position
    # Then create a list of tuples of (m6A/A ratio, read depth) for each position, this will allow me to calculate a weighted average for the m6A/A ratio
    # Split this into canonical and non-canonical mod positions, where canonical mod positions are defined by the m6A/A ratio being above the CANNONICAL_MOD_PROP_THRESHOLD
    canonical_mods = [(a, b, c, d) for a, b, c, d in mod_ratio_read_depth_tuples if a >= canonical_mod_prop_threshold and b >= read_depth_threshold]
    non_canonical_mods = [(a, b, c, d) for a, b, c, d in mod_ratio_read_depth_tuples if a < canonical_mod_prop_threshold or b < read_depth_threshold]

    if len(canonical_mods) > 0:
        weighted_average_canonical_mod = sum(a * b for a, b, c, d in canonical_mods) / sum([b for a, b, c, d in canonical_mods])
        if math.isnan(weighted_average_canonical_mod):
            weighted_average_canonical_mod = 0
    else:
        weighted_average_canonical_mod = 0

    if len(non_canonical_mods) > 0:
        weighted_average_non_canonical_mod = sum(a * b for a, b, c, d in non_canonical_mods) / sum([b for a, b, c, d in non_canonical_mods])
        if math.isnan(weighted_average_non_canonical_mod):
            weighted_average_non_canonical_mod = 0
    else:
        weighted_average_non_canonical_mod = 0

    sum_a = sum(d_coverage['count_a'])
    if sum_a > 0:
        total_mod_to_unmodified_ratio = sum(d_coverage['m6A']) / sum_a
    else:
        total_mod_to_unmodified_ratio = 0

    row_summary = [label, row.ID, row.type, row.strand, max_read_depth, [d for a, b, c, d in canonical_mods], weighted_average_canonical_mod, weighted_average_non_canonical_mod, total_mod_to_unmodified_ratio]
    print("\t".join(map(str, row_summary)))

    return row_summary


# NOTE i do not think this is base and read accurate
# for each mod, it get close to the methylation ratio, but not exact (when inspected in JBrowse)
def gene_methylation_analysis(args):
    INPUT = args.input
    ANNOTATION_FILE = args.annotation
    POLY_A_FILTER = args.poly_a_filter
    COVERAGE_PADDING = args.padding
    CANNONICAL_MOD_PROP_THRESHOLD = args.mod_ratio
    READ_DEPTH_THRESHOLD = args.read_depth
    FEATURE_TYPE = args.type
    VERBOSE = args.verbose
    IDS = args.ids
    GENOME = args.genome
    MOD_PROB_THRESHOLD = args.mod_prob_threshold
    OUTPUT = args.output
    COMPARE_METHYLATION_BETWEEN_TREATMENTS = args.compare_methylation_between_treatments

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)

    print("LOG: starting gene methylation analysis")
    print("LOG: input files: {}".format(INPUT))
    print("LOG: annotation file: {}".format(ANNOTATION_FILE))
    print("LOG: poly A filter: {}".format(POLY_A_FILTER))
    print("LOG: coverage padding: {}".format(COVERAGE_PADDING))
    print("LOG: canonical mod prop threshold: {}".format(CANNONICAL_MOD_PROP_THRESHOLD))
    print("LOG: read depth threshold: {}".format(READ_DEPTH_THRESHOLD))
    print("LOG: feature type: {}".format(FEATURE_TYPE))
    print("LOG: ids: {}".format(IDS))
    
    input_files = process_input_files(INPUT)
    gff_df = process_annotation_file(ANNOTATION_FILE)

    if FEATURE_TYPE:
        matches = gff_df[gff_df['type'] == FEATURE_TYPE]
    else:
        matches = gff_df[gff_df['ID'].isin(IDS)]
    if matches.empty:
        print("ERROR: no matches found for type {} and ids {}".format(FEATURE_TYPE, IDS))

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']

    raw_summary_header = [
        "label",
        "ID",
        "type",
        "strand",
        "max_depth",
        "canonical_mods",
        "canonical_wam",
        "non_canonical_wam",
        "total_wam"
    ]
    print("\t".join(map(str, raw_summary_header)))

    raw_results = []

    for label in bam_labels:
        with ThreadPoolExecutor() as executor:
            process_row_partial = partial(
                process_row,
                label=label,
                samfile_path=input_files[label]['path'],
                coverage_padding=COVERAGE_PADDING,
                pysam_mod_threshold=PYSAM_MOD_THRESHOLD,
                read_depth_threshold=READ_DEPTH_THRESHOLD,
                canonical_mod_prop_threshold=CANNONICAL_MOD_PROP_THRESHOLD
            )

            label_raw_results = list(executor.map(process_row_partial, [row for _, row in matches.iterrows()]))

        raw_results += label_raw_results

    raw_results_df = pandas.DataFrame(raw_results, columns=raw_summary_header)

    if OUTPUT:
        raw_results_df.to_csv("raw_{}".format(OUTPUT), sep='\t', index=False)

    if COMPARE_METHYLATION_BETWEEN_TREATMENTS:
        summary_header = [
            "ID",
            "type",
            "strand",
            "canonical_mods",
            "average_depth_g1",
            "average_depth_g2",
            "canonical_wam_g1",
            "canonical_wam_g2",
            "non_canonical_wam_g1",
            "non_canonical_wam_g2",
            "total_wam_g1",
            "total_wam_g2",
            "test_stat",
            "p_val"
        ]
        print("\t".join(map(str, summary_header)))
        results = []
        group_1_prefix = COMPARE_METHYLATION_BETWEEN_TREATMENTS[0]
        group_2_prefix = COMPARE_METHYLATION_BETWEEN_TREATMENTS[1]

        group1_rows = raw_results_df[raw_results_df['label'].str.startswith(group_1_prefix)]
        group2_rows = raw_results_df[raw_results_df['label'].str.startswith(group_2_prefix)]

        group1_bam_labels = group1_rows.label.to_list()
        group2_bam_labels = group2_rows.label.to_list()

        # go through each row, determine quantify change in 
        for row_index, row in matches.iterrows():
            these_rows_g1 = group1_rows[group1_rows['ID'] == row.ID]
            these_rows_g2 = group2_rows[group2_rows['ID'] == row.ID]

            canonical_mods = []

            for _, row_result in raw_results_df[raw_results_df.ID == row.ID].iterrows():
                canonical_mods.append(row_result.canonical_mods)

            canonical_mods = list(set([item for sublist in canonical_mods for item in sublist]))

            average_depth_g1 = int(sum(these_rows_g1.max_depth.to_list()) / len(these_rows_g1.max_depth.to_list()))
            average_depth_g2 = int(sum(these_rows_g2.max_depth.to_list()) / len(these_rows_g2.max_depth.to_list()))

            canonical_wam_g1 = sum(v * w for v, w in zip(these_rows_g1.canonical_wam.to_list(), these_rows_g1.max_depth.to_list())) / sum(these_rows_g1.max_depth.to_list())
            canonical_wam_g2 = sum(v * w for v, w in zip(these_rows_g2.canonical_wam.to_list(), these_rows_g2.max_depth.to_list())) / sum(these_rows_g2.max_depth.to_list())

            non_canonical_wam_g1 = sum(v * w for v, w in zip(these_rows_g1.non_canonical_wam.to_list(), these_rows_g1.max_depth.to_list())) / sum(these_rows_g1.max_depth.to_list())
            non_canonical_wam_g2 = sum(v * w for v, w in zip(these_rows_g2.non_canonical_wam.to_list(), these_rows_g2.max_depth.to_list())) / sum(these_rows_g2.max_depth.to_list())

            total_wam_g1 = sum(v * w for v, w in zip(these_rows_g1.total_wam.to_list(), these_rows_g1.max_depth.to_list())) / sum(these_rows_g1.max_depth.to_list())
            total_wam_g2 = sum(v * w for v, w in zip(these_rows_g2.total_wam.to_list(), these_rows_g2.max_depth.to_list())) / sum(these_rows_g2.max_depth.to_list())

            group_depth = these_rows_g1.max_depth.to_list() + these_rows_g2.max_depth.to_list()
            group_canonical_wam = these_rows_g1.canonical_wam.to_list() + these_rows_g2.canonical_wam.to_list()
            approx_modified = [v * w for v, w in zip(group_depth, group_canonical_wam)]

            # from chatGPT
            data = pandas.DataFrame({
                'group': [group_1_prefix] * len(group1_bam_labels) + [group_2_prefix] * len(group2_bam_labels),
                'approx_modified': approx_modified,
                'depth': group_depth
            })

            data['group_binary'] = (data['group'] == group_2_prefix).astype(int)

            # Failures column
            data['approx_unmodified'] = data['depth'] - data['approx_modified']
            response = data[['approx_modified', 'approx_unmodified']]

            # Fit Binomial GLM
            model = sm.GLM(response, sm.add_constant(data['group_binary']), family=sm.families.Binomial())
            result = model.fit()
            # print(result.summary())
            test_stat = result.params.group_binary
            pval = result.pvalues.group_binary

            row_summary = [row.ID, row.type, row.strand, canonical_mods, average_depth_g1, average_depth_g2, canonical_wam_g1, canonical_wam_g2, non_canonical_wam_g1, non_canonical_wam_g2, total_wam_g1, total_wam_g2, test_stat, pval]
            results.append(row_summary)
            print("\t".join(map(str, row_summary)))

        summary_df = pandas.DataFrame(results, columns=summary_header)

        if OUTPUT:
            summary_df.to_csv(OUTPUT, sep='\t', index=False)



            