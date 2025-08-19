import pandas
import pysam
import numpy
from rqc_modules.utils import getSubfeatures, process_input_files, process_annotation_file
from rqc_modules.constants import PYSAM_MOD_TUPLES, MODKIT_BEDMETHYL_HEADER, FEATURECOUNTS_HEADER, BAM_PILEUP_DEFAULT_FLAGS, PYSAM_PILEUP_MAX_DEPTH, BAM_REVERSE_STRAND

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
    

    input_files = process_input_files(INPUT)
    gff_df = process_annotation_file(ANNOTATION_FILE)

    if FEATURE_TYPE:
        matches = gff_df[gff_df['type'] == FEATURE_TYPE]
    else:
        matches = gff_df[gff_df['ID'].isin(IDS)]
    if matches.empty:
        print("ERROR: no matches found for type {} and ids {}".format(FEATURE_TYPE, IDS))

    # END SETUP

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']
    bedmethyl_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bedmethyl']


    # START CANNONICAL MOD IDENTIFICATION
    print("LOG: identifying canonical mod positions in bedmethyl files")
    cannonical_mods_start_pos = {}
    for label in bam_labels:
        prefix = label.split("_")[0]

        # HACK fixme to determine which bedmethyl is associated with this bam file
        mod_label = "{}_m6A_0.95".format(prefix)

        mods_file_df = pandas.read_csv(input_files[mod_label]['path'], sep='\t', names=MODKIT_BEDMETHYL_HEADER)
        mods_file_df['strand'] = mods_file_df['strand'].astype('category')
        mods_file_df['contig'] = mods_file_df['contig'].astype('category')
        
        # TODO maybe don't need to filter this for read depth, just filter the gene for read depth
        mods_file_df = mods_file_df[
            (mods_file_df.percent_mod >= (CANNONICAL_MOD_PROP_THRESHOLD * 100)) & 
            (mods_file_df.valid_cov >= READ_DEPTH_THRESHOLD)
        ]

        # TODO update this to use my method described in labarchives
        for row_index, row in matches.iterrows():

            if row['ID'] not in cannonical_mods_start_pos:
                cannonical_mods_start_pos[row['ID']] = []

            row_mods = mods_file_df[
                (mods_file_df.start >= (row['start'] - COVERAGE_PADDING)) &
                (mods_file_df.end <= (row['end'] + COVERAGE_PADDING)) &
                (mods_file_df.strand == row['strand']) &
                (mods_file_df.contig == row['seq_id'])
            ]

            for mod_index, mod in row_mods.iterrows():
                if mod['start'] not in cannonical_mods_start_pos[row['ID']]:
                    cannonical_mods_start_pos[row['ID']].append(mod['start'])



    # END CANNONICAL MOD IDENTIFICATION
    print("label\tgene id\treads used\treads in region\tfiltered (strand)\tfiltered (fc)\tfiltered (3p)\tfiltered (mod)")

    d_mod_info = {}

    for label in bam_labels:
        samfile = pysam.AlignmentFile(input_files[label]['path'], 'rb')
        d_mod_info[label] = {}
        
        # attempt to find the relevent featureCounts file in input_files
        # TODO load cannonical mod positions into array and convert to tx space
        # feature_counts_sample_label = label.split("_")[0] + "_featureCounts"
        # feature_counts_df = pandas.read_csv(input_files[feature_counts_sample_label]['path'], sep='\t', names=FEATURECOUNTS_HEADER)
        # feature_counts_df['targets'] = feature_counts_df['targets'].astype('category')

        # # TODO load cannonical mod positions into array and convert to tx space
        # prefix = label.split("_")[0]
        # mod_label = "{}_m6A_0.95".format(prefix)

        # mods_file_df = pandas.read_csv(input_files[mod_label]['path'], sep='\t', names=MODKIT_BEDMETHYL_HEADER)
        # mods_file_df['strand'] = mods_file_df['strand'].astype('category')
        # mods_file_df['contig'] = mods_file_df['contig'].astype('category')


        # generate coverage for all matches in this bam file
        for row_index, row in matches.iterrows():
            summary_df_index = 0
            read_on_different_strand = 0
            gene_length = row['end'] - row['start']
            row_name = row['ID']
            # gene_reads = feature_counts_df[feature_counts_df.targets == row['ID'].split(".")[0]]

            reads_in_region = samfile.fetch(
                contig=row['seq_id'], 
                start=row['start'] - COVERAGE_PADDING, 
                stop=row['end'] + COVERAGE_PADDING
            )
            reads_in_region = list(reads_in_region)

            # STOP_CLOCK("fetch", "stop")

            # 0.12086892127990723s
            print(cannonical_mods_start_pos[row['ID']])

            # !!!!! START NANOPORE SPECIFIC !!!!!
            # filter out reads where the 3' end is not in or beyond the last feature (3'UTR or last exon) of the target gene
            row_subfeatures = getSubfeatures(gff_df, row['ID'], "subfeature", 0)

            read_indexes_to_process = []

            read_outside_3p_end = []

            MOD_PROB_THRESHOLD = 0.95
            PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD) 

            # example: PF3D7_0709050.1
            if len(row_subfeatures) == 0:
                most_3p_subfeature = row
            else:
                if row['strand'] == "-":
                    most_3p_subfeature = row_subfeatures.iloc[0]
                else:
                    most_3p_subfeature = row_subfeatures.iloc[-1]

            num_a = 0
            num_mod = 0
            d_this_row_mod_counts = {}

            # generate coverage for this gene
            row_flag_filters = BAM_PILEUP_DEFAULT_FLAGS
            row_flags_requires = 0

            if row['strand'] == '+':
                row_flag_filters = BAM_PILEUP_DEFAULT_FLAGS | BAM_REVERSE_STRAND
            else:
                row_flags_requires = BAM_REVERSE_STRAND

            gene_length = row['end'] - row['start'] + 1
            d_coverage = {
                "total_depth": [0] * gene_length,
                "count_a": [0] * gene_length,
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
                truncate = True
            ):
                # take the number of aligned reads at this column position (column.n) minus the number of aligned reads which have either a skip or delete base at this column position (r.query_position)
                column_bases_read = list(filter(None, column.get_query_sequences()))
                
                d_coverage['total_depth'][column.reference_pos - row['start'] + 1] = len(column_bases_read)
                d_coverage['count_a'][column.reference_pos - row['start'] + 1] = column_bases_read.count('A')

            # NOTE: now the longest function in the TES analysis
            for i in range(len(reads_in_region)):
                r = reads_in_region[i]

                # keep only reads in the same direction as this strand
                if (row['strand'] == "+" and r.is_forward) or (row['strand'] == "-" and r.is_reverse):
                    if row['strand'] == "-":
                        read_3p_end = r.reference_start

                        if read_3p_end <= most_3p_subfeature.end:
                            # this is [(read index, 256 * mod_prob)...]
                            mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES['m6A_rev'])
                            if mods_probs:
                                # keep only mod positions which are above mod prob threshold
                                ref_pos = numpy.array(r.get_reference_positions(full_length=True))
                                read_mod_positions = [x[0] for x in mods_probs if x[1] >= PYSAM_MOD_THRESHOLD]
                                # print(ref_pos[read_mod_positions])

                            else:
                                read_indexes_to_process.append(i)
                        else:
                            read_outside_3p_end.append(r.query_name)

                    else:
                        read_3p_end = r.reference_end

                        if read_3p_end >= most_3p_subfeature.start:
                            num_a += r.query_sequence.count('A')

                            ref_pos = numpy.array(r.get_reference_positions(full_length=True))

                            mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES['m6A_for'])
                            if mods_probs:
                                # keep only mod positions which are above mod prob threshold
                                read_mod_positions = [x[0] for x in mods_probs if x[1] >= PYSAM_MOD_THRESHOLD]
                                genomic_mod_positions = ref_pos[read_mod_positions]

                                for mod in genomic_mod_positions:
                                    if mod == None:
                                        print(ref_pos)
                                        print(read_mod_positions)

                                    d_coverage['m6A'][mod - row['start'] + 1] += 1

                                # TODO keep only mod positions within the gene
                                num_mod += len(genomic_mod_positions)
                        else:
                            read_outside_3p_end.append(r.query_name)
                else:
                    read_on_different_strand += 1

            # convert d_coverage to genomic space
            # d_coverage['total_depth'] = [x + row.start for x in d_coverage['total_depth']]

            print("num_a: {}, num_mod: {}".format(num_a, num_mod))
            print(d_this_row_mod_counts)
        samfile.close()

    # we now have gone through each gene and calculated canonical mod depth, and non canonical mod depth, and total mod/unmodified ratio
    # gene id | type | strand | sample1 - list of canonical mod positions (genomic) | sample1 - total m6A/A | sample2 - list of canonical mod positions (genomic) | sample2 - total m6A/A | cMod average methylation (weighted by cMod read depth)) | non-cMod average methylation

    # canonical mods are determined by the modkit bedmethyl file, which is generated by modkit


            