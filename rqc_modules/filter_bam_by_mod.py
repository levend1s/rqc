def filter_bam_by_mod(args):

    if len(cannonical_mods_genome_space_filter):
        FILTER_READS_FOR_CANNONICAL_MODS = True
    else:
        FILTER_READS_FOR_CANNONICAL_MODS = False

    # START CANNONICAL MOD IDENTIFICATION
    cannonical_mods_start_pos = {}
    for label in bam_labels:
        prefix = label.split("_")[0]
        mod_label = "{}_m6A_0.95".format(prefix)

        if DEBUG:
            print("CANNONICAL_MOD_PROP_THRESHOLD: {}".format(CANNONICAL_MOD_PROP_THRESHOLD))
            print("READ_DEPTH_THRESHOLD: {}".format(READ_DEPTH_THRESHOLD))

        mods_file_df = pandas.read_csv(input_files[mod_label]['path'], sep='\t', names=MODKIT_BEDMETHYL_HEADER)
        mods_file_df['strand'] = mods_file_df['strand'].astype('category')
        mods_file_df['contig'] = mods_file_df['contig'].astype('category')
        
        # TODO maybe don't need to filter this for read depth, just filter the gene for read depth
        mods_file_df = mods_file_df[
            (mods_file_df.percent_mod >= (CANNONICAL_MOD_PROP_THRESHOLD * 100)) & 
            (mods_file_df.valid_cov >= READ_DEPTH_THRESHOLD)
        ]

        for row_index, row in matches.iterrows():

            if row['ID'] not in cannonical_mods_start_pos:
                cannonical_mods_start_pos[row['ID']] = []

            row_mods = mods_file_df[
                (mods_file_df.start >= (row['start'] - COVERAGE_PADDING)) &
                (mods_file_df.end <= (row['end'] + COVERAGE_PADDING)) &
                (mods_file_df.strand == row['strand']) &
                (mods_file_df.contig == row['seq_id'])
            ]

            if DEBUG:
                print(row_mods)

            for mod_index, mod in row_mods.iterrows():
                if mod['start'] not in cannonical_mods_start_pos[row['ID']]:
                    cannonical_mods_start_pos[row['ID']].append(mod['start'])

    if DEBUG:
        print("cannonical_mods_start_pos: {}".format(cannonical_mods_start_pos[row['ID']]))

    # END CANNONICAL MOD IDENTIFICATION
    for label in bam_labels:
        samfile = pysam.AlignmentFile(input_files[label]['path'], 'rb')
        num_bams += 1
        d_poly_a_lengths[label] = {}
        d_tts[label] = {}
        d_mod_info[label] = {}
        
        d_not_beyond_3p[label] = {}
        d_not_in_feature_counts[label] = {}

        # attempt to find the relevent featureCounts file in input_files
        feature_counts_sample_label = label.split("_")[0] + "_featureCounts"
        feature_counts_df = pandas.read_csv(input_files[feature_counts_sample_label]['path'], sep='\t', names=FEATURECOUNTS_HEADER)
        feature_counts_df['targets'] = feature_counts_df['targets'].astype('category')

        # TODO load cannonical mod positions into array and convert to tx space
        prefix = label.split("_")[0]
        mod_label = "{}_m6A_0.95".format(prefix)

        mods_file_df = pandas.read_csv(input_files[mod_label]['path'], sep='\t', names=MODKIT_BEDMETHYL_HEADER)
        mods_file_df['strand'] = mods_file_df['strand'].astype('category')
        mods_file_df['contig'] = mods_file_df['contig'].astype('category')

        # generate coverage for all matches in this bam file
        for row_index, row in matches.iterrows():
            # find subfeatures
            # 0.8235001564025879s - 1.4s
            # START_CLOCK("row_start")

            summary_df_index = 0
            read_on_different_strand = 0

            # 0.30355286598205566s
            gene_reads = feature_counts_df[feature_counts_df.targets == row['ID'].split(".")[0]]
            # START_CLOCK("fetch")

            # 0.0003178119659423828s
            reads_in_region = samfile.fetch(
                contig=row['seq_id'], 
                start=row['start'] - COVERAGE_PADDING, 
                stop=row['end'] + COVERAGE_PADDING
            )
            reads_in_region = list(reads_in_region)

            gene_length = row['end'] - row['start']
            row_name = row['ID']
            # STOP_CLOCK("fetch", "stop")

            # 0.12086892127990723s
            row_mods = mods_file_df[
                (mods_file_df['start'].isin(cannonical_mods_start_pos[row['ID']])) &
                (mods_file_df.strand == row['strand']) &
                (mods_file_df.contig == row['seq_id'])
            ]

            d_mod_info[label][row['ID']] = {}
            d_mod_info[label][row['ID']]['valid_cov'] = {}
            d_mod_info[label][row['ID']]['num_mod'] = {}

            for mod_index, mod in row_mods.iterrows():
                d_mod_info[label][row['ID']]['valid_cov'][mod['start']] = mod['valid_cov']
                d_mod_info[label][row['ID']]['num_mod'][mod['start']] = mod['num_mod']

            # make sure all cannonical mod locations exist, so add them as zero if not in bedmethyl
            for mod in cannonical_mods_start_pos[row['ID']]:
                if mod not in d_mod_info[label][row['ID']]['valid_cov']:
                    d_mod_info[label][row['ID']]['valid_cov'][mod] = 0
                    d_mod_info[label][row['ID']]['num_mod'][mod] = 0

            # !!!!! START NANOPORE SPECIFIC !!!!!
            # filter out reads where the 3' end is not in or beyond the last feature (3'UTR or last exon) of the target gene
            row_subfeatures = getSubfeatures(row['ID'], "subfeature", 0)

            read_indexes_to_process = []

            mod_of_interest = 'm6A'

            missing_cannonical_mods = []
            read_outside_3p_end = []

            MOD_PROB_THRESHOLD = 0.95
            pysam_mod_threshold = int(256 * MOD_PROB_THRESHOLD) 

            # example: PF3D7_0709050.1
            if len(row_subfeatures) == 0:
                most_3p_subfeature = row
            else:
                if row['strand'] == "-":
                    most_3p_subfeature = row_subfeatures.iloc[0]
                else:
                    most_3p_subfeature = row_subfeatures.iloc[-1]

            # NOTE: now the longest function in the TES analysis
            for i in range(len(reads_in_region)):
                r = reads_in_region[i]

                # keep only reads in the same direction as this strand
                if (row['strand'] == "+" and r.is_forward) or (row['strand'] == "-" and r.is_reverse):
                    if row['strand'] == "-":
                        read_3p_end = r.reference_start

                        if read_3p_end <= most_3p_subfeature.end:
                            # read_indexes_to_process.append(this_index)

                            if FILTER_READS_FOR_CANNONICAL_MODS:
                                # this is [(read index, 256 * mod_prob)...]
                                mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES['m6A_rev'])
                                if mods_probs:
                                    # keep only mod positions which are above mod prob threshold
                                    ref_pos = numpy.array(r.get_reference_positions(full_length=True))
                                    read_mod_positions = [x[0] for x in mods_probs if x[1] >= pysam_mod_threshold]
                                    
                                    # read mod positions is the position from the start of the read
                                    # aligned reads mayu contain indels, so we need to get reference index from get_reference_positions
                                    # print(ref_pos[read_mod_positions])
                                    if set(cannonical_mods_genome_space_filter).issubset(ref_pos[read_mod_positions]):
                                        # print("READ HAD ALL CANNONICAL MODS\n")
                                        if FILTER_FOR_M6A != "[]":
                                            read_indexes_to_process.append(i)
                                        if FILTER_OUT_M6A != "[]":
                                            missing_cannonical_mods.append(r.query_name)
                                    else:
                                        # print("READ DID NOT HAVE ALL CANNONICAL MODS: id={}".format(r.query_name))
                                        if FILTER_FOR_M6A != "[]":
                                            missing_cannonical_mods.append(r.query_name)
                                        if FILTER_OUT_M6A != "[]":
                                            read_indexes_to_process.append(i)
                            else:
                                read_indexes_to_process.append(i)
                        else:
                            read_outside_3p_end.append(r.query_name)

                    else:
                        read_3p_end = r.reference_end

                        if read_3p_end >= most_3p_subfeature.start:
                            # read_indexes_to_process.append(this_index)
                            if FILTER_READS_FOR_CANNONICAL_MODS:
                                # this is [(read index, 256 * mod_prob)...]
                                mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES['m6A_for'])
                                if mods_probs:
                                    # keep only mod positions which are above mod prob threshold
                                    ref_pos = numpy.array(r.get_reference_positions(full_length=True))
                                    read_mod_positions = [x[0] for x in mods_probs if x[1] >= pysam_mod_threshold]
                                    
                                    # read mod positions is the position from the start of the read
                                    # aligned reads mayu contain indels, so we need to get reference index from get_reference_positions
                                    # print(ref_pos[read_mod_positions])
                                    if set(cannonical_mods_genome_space_filter).issubset(ref_pos[read_mod_positions]):
                                        # print("READ HAD ALL CANNONICAL MODS\n")
                                        if FILTER_FOR_M6A != "[]":
                                            read_indexes_to_process.append(i)
                                        if FILTER_OUT_M6A != "[]":
                                            missing_cannonical_mods.append(r.query_name)
                                    else:
                                        # print("READ DID NOT HAVE ALL CANNONICAL MODS: id={}".format(r.query_name))
                                        if FILTER_FOR_M6A != "[]":
                                            missing_cannonical_mods.append(r.query_name)
                                        if FILTER_OUT_M6A != "[]":
                                            read_indexes_to_process.append(i)
                            else:
                                read_indexes_to_process.append(i)
                        else:
                            read_outside_3p_end.append(r.query_name)
                else:
                    read_on_different_strand += 1

            if GENERATE_FILTERED_BAM:
                if FILTER_FOR_M6A != "[]":
                    fs_str = "FILTER_FOR_M6A"
                else:
                    fs_str = "FILTER_OUT_M6A"

                filtered_bam_filename = "{}_rqc_{}_{}.bam".format(label, fs_str, str(cannonical_mods_genome_space_filter))
                print("GENERATING FILTERED BAM: {}".format(filtered_bam_filename))

                rqc_filtered = pysam.AlignmentFile(filtered_bam_filename, "wb", template=samfile)
                for i in read_indexes_to_process:
                    rqc_filtered.write(reads_in_region[i])

                rqc_filtered.close()