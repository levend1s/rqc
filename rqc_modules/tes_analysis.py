def tes_analysis(args):
    # load annotation file
    feature_id = args.inputs[0]

    # process input file. Each line contains a label, the type of file, and the filepath
    # label group type path
    input_files = {}
    # if len(args.inputs[1:]) % 4 != 0:
    #     print("ERROR: not enough information in specified files. Check each input follows the format [LABEL] [TYPE] [PATH]")
    #     sys.exit()
    # else:
    in_index = 1
    while in_index < len(args.inputs):
        if not args.inputs[in_index].startswith("#"):
            input_files[args.inputs[in_index]] = {
                'group': args.inputs[in_index+1], 
                'type': args.inputs[in_index+2], 
                'path': args.inputs[in_index+3]
            }
        in_index += 4

    # load annotation file and find indexes for all parent children
    ANNOTATION_FILE = gffpandas.read_gff3(ANNOTATION_FILE_PATH)
    GFF_DF = ANNOTATION_FILE.attributes_to_columns()

    for row_index, row in GFF_DF.iterrows():
        if row['Parent'] in GFF_PARENT_TREE:
            GFF_PARENT_TREE[row['Parent']].append(row_index)
        else:
            GFF_PARENT_TREE[row['Parent']] = [row_index]

    # Try to find matches of provided type, if not, assume that input is a list of IDs
    if os.path.isfile(feature_id):
        lines = []
        with open(feature_id) as f:
            lines = f.read().splitlines()

        matches = ANNOTATION_FILE.get_feature_by_attribute("ID", lines)

    else:
        matches = ANNOTATION_FILE.filter_feature_of_type([feature_id])
        if len(matches.df) == 0:
            evaluated_input = ast.literal_eval(feature_id)
            matches = ANNOTATION_FILE.get_feature_by_attribute("ID", evaluated_input)
            print("Looking for {} IDs, found {} matches: TES analysis for {}".format(len(evaluated_input) , len(matches.df), evaluated_input))
        else:
            print("Found {} matches for type {}. Calculating TES variance...".format(len(matches.df), feature_id))

    matches = matches.attributes_to_columns()

    SINGLE_GENE_ANALYSIS = False
    if len(matches) == 1:
        SINGLE_GENE_ANALYSIS = True


    num_bams = 0

    d_poly_a_lengths = {}
    d_tts = {}
    d_mod_info = {}
    
    d_not_beyond_3p = {}
    d_not_in_feature_counts = {}

    logfile_df_index = 0

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']
    bam_labels_control = [l for l in bam_labels if input_files[l]['group'] == 'control']
    bam_labels_treatment = [l for l in bam_labels if input_files[l]['group'] == 'knock-sideways']

    summary_df = pandas.DataFrame(columns=TES_SUMMARY_HEADER)

    gene_length = 0

    cannonical_mods_genome_space_filter = []

    if FILTER_FOR_M6A != "[]":
        cannonical_mods_genome_space_filter = ast.literal_eval(FILTER_FOR_M6A)
        print("FILTERING FOR READS CONTAINING: {}".format(cannonical_mods_genome_space_filter))
    if FILTER_OUT_M6A != "[]":
        cannonical_mods_genome_space_filter = ast.literal_eval(FILTER_OUT_M6A)
        print("FILTERING FOR READS CONTAINING: {}".format(cannonical_mods_genome_space_filter))

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
    print("label\tgene id\treads used\treads in region\tfiltered (strand)\tfiltered (fc)\tfiltered (3p)\tfiltered (mod)")

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


            gene_read_ids_fc = set(gene_reads['read_id'])

            found = 0
            not_found = 0
            missing_from_fc = 0
            no_poly_a = 0
            poly_a_lengths = []
            tts_sites = []

            for i in read_indexes_to_process:
                r = reads_in_region[i]
                if FILTER_FOR_FEATURE_COUNTS and (r.qname not in gene_read_ids_fc):
                    missing_from_fc += 1
                    continue
                else:
                    if row['strand'] == "-":
                        tts_sites.append(r.reference_start)
                    else:
                        tts_sites.append(r.reference_end)

                    # if SINGLE_GENE_ANALYSIS:
                    if r.has_tag('pt:i'):
                        poly_a_length = r.get_tag('pt:i')
                        poly_a_lengths.append(poly_a_length)
                    else:
                        # print("WARNING: {} does not have poly a tag".format(r.qname))
                        no_poly_a += 1
                        poly_a_lengths.append(0)


            d_poly_a_lengths[label][row['ID']] = poly_a_lengths
            d_tts[label][row['ID']] = tts_sites

            d_not_beyond_3p[label][row['ID']] = len(read_outside_3p_end)
            d_not_in_feature_counts[label][row['ID']] = missing_from_fc
            # STOP_CLOCK("row_start", "stop")

            # print("label\tgene id\treads used\treads in region\tfiltered (strand)\tfiltered (fc)\tfiltered (3p)\tfiltered (mod)")
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(label, row['ID'], len(tts_sites), len(reads_in_region), read_on_different_strand, missing_from_fc, len(read_outside_3p_end), len(missing_cannonical_mods)))


        samfile.close()

    # We are interested if the TES has multiple end sites 
    # OR
    # if the intratretment samples are the same (p > 0.05) and the intertreatment samples are different (p < 0.05) after correcting bonferri
    pairwise_combinations_same_treatment = list(itertools.combinations(bam_labels_control, 2)) + list(itertools.combinations(bam_labels_treatment, 2))
    
    pairwise_combinations_inter_treatment = list(itertools.combinations(bam_labels, 2))
    pairwise_combinations_inter_treatment = [c for c in pairwise_combinations_inter_treatment if c not in pairwise_combinations_same_treatment]
    first_label = bam_labels[0]

    ### CALCULTE METHYLATION CHANGE
    # d_mod_info
    # can calculate % difference in cannonical modified sites
    # can calculate the log ratio
    print("{} - Calculating methylation change...".format(row['ID']))
    weighted_mod_ratios_before = {}
    weighted_mod_ratios_after = {}
    d_wam_before = {}
    d_wam_after = {}
    groups = [(bam_labels_control, weighted_mod_ratios_before), (bam_labels_treatment, weighted_mod_ratios_after)]
    d_wam_change = {}
    d_x_ticks = {}
    d_poly_a_length_hists = {}
    d_tts_hist = {}
    d_cdfs = {}
    d_kdes = {}
    d_max_hist_count_poly_a = {}
    d_max_hist_count_tts = {}
    d_max_poly_a = {}
    d_min_tts = {}
    d_max_tts = {}
    d_max_density = {}

    for label in bam_labels:
        d_poly_a_length_hists[label] = {}
        d_tts_hist[label] = {}
        d_kdes[label] = {}
        d_cdfs[label] = {}

    for row_index, row in matches.iterrows():
        gene_id = row['ID']
        weighted_mod_ratios_before[gene_id] = []
        weighted_mod_ratios_after[gene_id] = []

        cannonical_mod_keys = d_mod_info[first_label][gene_id]['valid_cov'].keys()

        # calculate sum of total cannonical mods read depth across all samples for weighting
        for group_labels, group_weighted_outputs in groups:
            valid_cov_total = 0
            for label in group_labels:
                for cannonical_mod in cannonical_mod_keys:
                    if cannonical_mod in d_mod_info[label][gene_id]['valid_cov']:
                        valid_cov_total += d_mod_info[label][gene_id]['valid_cov'][cannonical_mod]
                    else:
                        print("WARNING: {} not in {}-{}".format(cannonical_mod, label, gene_id))
                        print(d_mod_info[label][gene_id])
                        valid_cov_total += 0

            for label in group_labels:            
                for cannonical_mod in cannonical_mod_keys:
                    num_valid = d_mod_info[label][gene_id]['valid_cov'][cannonical_mod]

                    if num_valid > 0:
                        num_mods = d_mod_info[label][gene_id]['num_mod'][cannonical_mod]
                        weight = num_valid / valid_cov_total

                        this_weighted_mod_proportion = (num_mods / num_valid) * weight
                        group_weighted_outputs[gene_id].append(this_weighted_mod_proportion)
                    else:
                        group_weighted_outputs[gene_id].append(0)

    for gene in weighted_mod_ratios_before:
        wam_before = sum(weighted_mod_ratios_before[gene])# / (len(weighted_mod_ratios_before[gene]))
        wam_after = sum(weighted_mod_ratios_after[gene])# / (len(weighted_mod_ratios_after[gene]))
        d_wam_before[gene] = wam_before
        d_wam_after[gene] = wam_after

        if wam_after == 0 or wam_before == 0:
            d_wam_change[gene] = 0
        else:
            wam_change = wam_after / wam_before
            d_wam_change[gene] = wam_change

        if DEBUG:
            print("gene {}: before {}, after: {}, change: {}".format(gene, wam_before, wam_after, d_wam_change[gene]))

    if DEBUG:
        pprint("weighted_mod_ratios_before: {}".format(weighted_mod_ratios_before))
        pprint("weighted_mod_ratios_after: {}".format(weighted_mod_ratios_after))
        pprint("d_wam_change: {}".format(d_wam_change))
    ### END METHYLATION CHANGE CALCULATION

    summary_df_index = 0
    for row_index, row in matches.iterrows():
        SAMPLE_HAS_LOW_EXP = False
        average_expression = 0
        average_not_beyond_3p = 0
        average_not_in_feature_counts = 0

        for label in bam_labels:
            if len(d_tts[label][row['ID']]) < READ_DEPTH_THRESHOLD:
                SAMPLE_HAS_LOW_EXP = True

            average_expression += len(d_tts[label][row['ID']])
            average_not_beyond_3p += d_not_beyond_3p[label][row['ID']]
            average_not_in_feature_counts += d_not_in_feature_counts[label][row['ID']]

        average_expression = math.floor(average_expression / len(bam_labels))
        average_not_beyond_3p = math.floor(average_not_beyond_3p / len(bam_labels))
        average_not_in_feature_counts = math.floor(average_not_in_feature_counts / len(bam_labels))

        if SAMPLE_HAS_LOW_EXP or average_expression < READ_DEPTH_THRESHOLD:
            # print pandas tsv row summary
            row_summary = [row['ID'], 0, 0, 0, 0, 0, 0, [], [], average_expression, cannonical_mods_start_pos[row['ID']], 0, 0, 0]
            summary_df.loc[summary_df_index] = row_summary
            summary_df_index += 1
            continue

        # calculate histograms for this gene
        # if SINGLE_GENE_ANALYSIS:
        gene_id = row['ID']

        max_hist_count_poly_a = 0
        max_hist_count_tts = 0
        max_poly_a = 0
        min_tts = 0
        max_tts = 0
        max_density = 0

        # find min, maxs, 
        for label in bam_labels:
            if min_tts > min(d_tts[label][gene_id]) or min_tts == 0:
                min_tts = min(d_tts[label][gene_id])

            if max_tts < max(d_tts[label][gene_id]):
                max_tts = max(d_tts[label][gene_id])

            if max_poly_a < max(d_poly_a_lengths[label][gene_id]):
                max_poly_a = max(d_poly_a_lengths[label][gene_id])
        
        x_ticks = range(min_tts, max_tts)

        # calculate hists
        print("{} - Generating transcript end site histograms...".format(row['ID']))
        for label in bam_labels:
            poly_a_length_range = list(range(1, max(d_poly_a_lengths[label][gene_id]) + 1))
            poly_a_hist = [d_poly_a_lengths[label][gene_id].count(i) for i in poly_a_length_range]
            
            unique_tes = set(d_tts[label][gene_id])
            tts_hist = [(d_tts[label][gene_id].count(i), i) for i in unique_tes]

            # split the tuple cause here we're interested in the biggest count in the hist
            e0 = [e[0] for e in tts_hist]
            if max_hist_count_tts < max(e0):
                max_hist_count_tts = max(e0)

            if max_hist_count_poly_a < max(poly_a_hist):
                max_hist_count_poly_a = max(poly_a_hist)

            d_poly_a_length_hists[label][row['ID']] = poly_a_hist
            d_tts_hist[label][row['ID']] = tts_hist



        # generate dennsity plots
        if SINGLE_GENE_ANALYSIS:
            print("{} - Generating transcript end site density information...".format(row['ID']))
            for label in bam_labels:
                kernel = scipy.stats.gaussian_kde(d_tts[label][gene_id])
                smoothed_tts_hist = kernel(x_ticks)
                cdf = numpy.cumsum(smoothed_tts_hist)

                d_kdes[label][row['ID']] = smoothed_tts_hist
                d_cdfs[label][row['ID']] = cdf

                if max_density < max(smoothed_tts_hist):
                    max_density = max(smoothed_tts_hist)

        d_max_hist_count_poly_a[gene_id] = max_hist_count_poly_a
        d_max_hist_count_tts[gene_id] = max_hist_count_tts
        d_max_poly_a[gene_id] = max_poly_a
        d_min_tts[gene_id] = min_tts
        d_max_tts[gene_id] = max_tts
        d_max_density[gene_id] = max_density
        d_x_ticks[row['ID']] = x_ticks

        # !!!! start of TES analysis, decide where the readthrough split point is
        # First, calculate the elbow for TES sites by count/frequency
        tes_variance_tests = ["z"]#, "x2", "mw-u", "ks"]

        d_read_through_counts = {}
        d_normal_read_counts = {}
        d_cannonical_tes = {}
        d_max_can_tes = {}
        d_sorted_tes = {}
        d_knee = {}
        d_tes_vs_prop = {}
        d_readthrough_split_points = {}
        d_fitted_curve_r_squared = {}
        d_fitted_curve_coeff = {}

        # calculate readthrough proportions for each sample
        print("{} - Finding max common transcript end site...".format(row['ID']))

        for label in bam_labels:
            d_tes_vs_prop[label] = []

            for c, p in d_tts_hist[label][row['ID']]:
                if row['strand'] == "-":
                    num_read_throughs = len([x for x in d_tts[label][gene_id] if x < p])
                    num_normal = len([x for x in d_tts[label][gene_id] if x >= p])
                else:
                    num_read_throughs = len([x for x in d_tts[label][gene_id] if x > p])
                    num_normal = len([x for x in d_tts[label][gene_id] if x <= p])
                
                rt_prop = num_read_throughs / num_normal

                # only keep this if there are fewer readthroughs than normals
                # NOTE this assumption doesn't work if there is one clean cut site for all reads
                # or the majority of reads
                if rt_prop < 1:
                    d_tes_vs_prop[label].append((p, rt_prop))

            # if we couldn't add any TES which had fewer readthroughs than normals
            # Then this must be a clean cut gene... The final TES must be the splitpoint
            if len(d_tes_vs_prop[label]) == 1:
                d_readthrough_split_points[label] = d_tes_vs_prop[label][0][0]
                d_fitted_curve_r_squared[label] = numpy.inf
                d_fitted_curve_coeff[label] = numpy.inf
                # TODO add type info to reporting dict
                continue

            # sort by genomic position
            sorted_tes_prop = sorted(d_tes_vs_prop[label], key=lambda a: a[0], reverse=False)
            pos = [x[0] for x in sorted_tes_prop]
            prop = [x[1] for x in sorted_tes_prop]

            # fit curve to TES data that will be used to to find the kneedle of the curve
            # data will have bumps and many false knees/elbows that we want to smooth out
            # so the kneedle function finds the correct knee

            # interpolate points between the few points we have so we can better approximate a
            # curve that fits our data
            pos_normalised = normalise_numpy_array(numpy.array(pos))
            x_interp = numpy.linspace(0, 1, 100)
            pos_y_interp = scipy.interpolate.interp1d(pos_normalised, prop)
            prop_normalised_interpolated = [pos_y_interp(x) for x in x_interp]

            # if strand is positive direction in kneedle function are reversed
            if row['strand'] == "-":
                elbow_direction = "increasing"
                initial_guess = [1, 0.1, 0]
            else:
                elbow_direction = "decreasing"
                initial_guess = [-1, 0.1, 1]

            # len(prop) < 100
            abc, pcov = scipy.optimize.curve_fit(
                power_func,
                x_interp,
                prop_normalised_interpolated,
                p0=initial_guess,
                maxfev=5000
            )

            y_fitted = power_func(x_interp, *abc)
            if DEBUG:
                print(abc)
                plt.scatter(x_interp, prop_normalised_interpolated)
                plt.scatter(x_interp, y_fitted)
                plt.show()
            # calculate R^2
            # https://stackoverflow.com/questions/19189362/getting-the-r-squared-value-using-curve-fit
            residuals = prop_normalised_interpolated - power_func(x_interp, *abc)
            residual_sum_squares = numpy.sum(residuals ** 2)
            total_sum_squares = numpy.sum((prop_normalised_interpolated - numpy.mean(prop_normalised_interpolated)) ** 2)
            r_squared = 1 - (residual_sum_squares / total_sum_squares)

            print("r_squared: {}".format(r_squared))

            fitted_curve_coeff = abc[1]
            d_fitted_curve_coeff[label] = fitted_curve_coeff

            # if we estimated a concave curve, take the final TES as end site
            if fitted_curve_coeff <= 1 and row['strand'] == "-":
                d_readthrough_split_points[label] = d_tes_vs_prop[label][0][0]
                d_fitted_curve_r_squared[label] = numpy.inf
                continue

            if fitted_curve_coeff >= 1 and row['strand'] == "+":
                d_readthrough_split_points[label] = d_tes_vs_prop[label][-1][0]
                d_fitted_curve_r_squared[label] = numpy.inf
                continue

            # as the genomic position increases, the splitpoint readthrough proportion gets bigger
            kneedle = KneeLocator(x_interp, y_fitted, S=1.0, curve='convex', direction=elbow_direction)

            # convert normalised elbow point back to genomic space
            if kneedle.knee:
                normalised_elbow_pos = kneedle.knee
            else:
                print("WARNING: {} - couldn't determine knee, setting as max TES".format(row['ID']))
                if row['strand'] == "-":
                    normalised_elbow_pos = 0
                else:
                    normalised_elbow_pos = 1

            genomic_elbow_pos = pos[0] + ((max(pos) - min(pos)) * normalised_elbow_pos)
            d_readthrough_split_points[label] = genomic_elbow_pos
            d_fitted_curve_r_squared[label] = r_squared

        print("{} - readthrough_split_points: {}".format(row['ID'], d_readthrough_split_points))

        # take the average of the control readthrough splitpoints and r^2
        readthrough_split_point = 0
        # average_r_squared = 0
        for label in bam_labels_control:
            readthrough_split_point += d_readthrough_split_points[label]
            # average_r_squared += d_fitted_curve_r_squared[label]

        readthrough_split_point = int(readthrough_split_point / len(bam_labels_control))
        # average_r_squared = int(average_r_squared / len(bam_labels_control))


        if DEBUG:
            fig, axes = plt.subplots()
            for label in bam_labels:
                # TODO add title of gene number, y axis labels etc, only do this if debug
                axes.scatter(*zip(*d_tes_vs_prop[label]), label=label, s=1)

                axes.axvline(x= readthrough_split_point, color='red', ls="--", linewidth=1.0)
                axes.legend()
            plt.show()

        # split into readthroughs and normals based on our TES
        print("{} - Finding readthrough proportions...".format(row['ID']))
        for label in bam_labels:
            if row['strand'] == "-":
                num_read_throughs = len([x for x in d_tts[label][gene_id] if x < readthrough_split_point])
                num_normal = len([x for x in d_tts[label][gene_id] if x >= readthrough_split_point])
            else:
                num_read_throughs = len([x for x in d_tts[label][gene_id] if x > readthrough_split_point])
                num_normal = len([x for x in d_tts[label][gene_id] if x <= readthrough_split_point])

            d_read_through_counts[label] = num_read_throughs
            d_normal_read_counts[label] = num_normal

        # CALCULATE WEIGHTED AVERAGE READTHROUGH RATIO
        weighted_rt_ratios_before = []
        weighted_rt_ratios_after = []
        d_wart_before = {}
        d_wart_after = {}
        groups = [(bam_labels_control, weighted_rt_ratios_before), (bam_labels_treatment, weighted_rt_ratios_after)]
        d_wart_change = {}

        # calculate sum of total cannonical mods read depth across all samples for weighting
        print("{} - Calculating weighted proportion change in readthroughs...".format(row['ID']))
        for group_labels, group_weighted_outputs in groups:
            valid_cov_total = 0
            for label in group_labels:
                valid_cov_total += d_normal_read_counts[label] + d_read_through_counts[label]

            for label in group_labels:
                weight = (d_normal_read_counts[label] + d_read_through_counts[label]) / valid_cov_total

                this_weighted_mod_proportion = (d_read_through_counts[label] / d_normal_read_counts[label]) * weight
                group_weighted_outputs.append(this_weighted_mod_proportion)

        rt_before = sum(weighted_rt_ratios_before)
        rt_after = sum(weighted_rt_ratios_after)
        d_wart_before[row['ID']] = rt_before
        d_wart_after[row['ID']] = rt_after

        if rt_before == 0 and rt_after > 0:
            d_wart_change[row['ID']] = numpy.inf
        elif rt_after == 0 or rt_before == 0:
            d_wart_change[row['ID']] = 0
        else:
            wart_change = rt_after / rt_before
            d_wart_change[row['ID']] = wart_change

        if DEBUG:
            print("gene {}: before {}, after: {}, change: {}".format(row['ID'], rt_before, rt_after, d_wart_change[row['ID']]))

        print("{} - Performing statistical test for change in readthrough proportions...".format(row['ID']))
        for test in tes_variance_tests:
            if average_expression < READ_DEPTH_THRESHOLD:
                p_inter_treatment = -1
                p_same_treatment = -1
                score = -1
                tests_passed = -1
            else:

                same_treatment_p_vals = []
                inter_treatment_p_vals = []

                # compare statistical tests
                # ks test
                # man whitney u test
                # chi squared test comparing KDEs

                for s1, s2 in pairwise_combinations_inter_treatment:
                    # two-sided: The null hypothesis is that the two distributions are identical, F(x)=G(x) for all x; the alternative is that they are not identical.
                    
                    # ks test
                    if test == "ks":
                        r = scipy.stats.ks_2samp(d_tts[s1][row['ID']], d_tts[s2][row['ID']])
                        print("{} vs {}: {}".format(s1, s2, r.pvalue))
                        print("{} vs {}".format(numpy.array(d_tts[s1][row['ID']]).mean(), numpy.array(d_tts[s2][row['ID']]).mean()))
                        pvalue = r.pvalue

                    if test == "x2":
                        print("{} and {}".format(len(d_tts_hist[s1][row['ID']]), len(d_tts_hist[s2][row['ID']])))
                        print("{} and {}".format(d_tts_hist[s1][row['ID']], d_tts_hist[s2][row['ID']]))

                        r = scipy.stats.chisquare(d_tts_hist[s1][row['ID']], f_exp=d_tts_hist[s2][row['ID']])
                        pvalue = r.pvalue
                    
                    if test == "z":
                        stat, pvalue = proportions_ztest(
                            count=[d_read_through_counts[s1], d_read_through_counts[s2]],
                            nobs=[d_normal_read_counts[s1], d_normal_read_counts[s2]],
                            alternative='two-sided'
                        )
                        if DEBUG:
                            print("{} vs {}: {}".format(s1, s2, pvalue))

                    # chi squared test
                    # print("{} and {}".format(len(d_tts[s1][row['ID']]), len(d_tts[s2][row['ID']])))
                    # r = scipy.stats.chisquare(d_tts_hist[s1], f_exp=d_tts_hist[s2])
                    inter_treatment_p_vals.append(pvalue)

                for s1, s2 in pairwise_combinations_same_treatment:
                    # two-sided: The null hypothesis is that the two distributions are identical, F(x)=G(x) for all x; the alternative is that they are not identical.
                    if test == "ks":
                        r = scipy.stats.ks_2samp(d_tts[s1][row['ID']], d_tts[s2][row['ID']])
                        print("{} vs {}: {}".format(s1, s2, r.pvalue))

                        print("{} vs {}".format(numpy.array(d_tts[s1][row['ID']]).mean(), numpy.array(d_tts[s2][row['ID']]).mean()))
                        pvalue = r.pvalue

                    if test == "x2":
                        r = scipy.stats.chisquare(d_tts_hist[s1][row['ID']], f_exp=d_tts_hist[s2][row['ID']])
                        pvalue = r.pvalue

                    if test == "z":
                        stat, pvalue = proportions_ztest(
                            count=[d_read_through_counts[s1], d_read_through_counts[s2]],
                            nobs=[d_normal_read_counts[s1], d_normal_read_counts[s2]],
                            alternative='two-sided'
                        )

                        if DEBUG:
                            print("{} vs {}: {}".format(s1, s2, pvalue))

                    same_treatment_p_vals.append(pvalue)

                # this is a debuious assumption of intratreatment ordering (idx = 0 and -1)
                alpha = 0.05
                num_pairs = 4
                bonferri_adjusted_alpha = alpha / num_pairs

                p_inter_treatment = scipy.stats.combine_pvalues(inter_treatment_p_vals).pvalue
                p_same_treatment = scipy.stats.combine_pvalues(same_treatment_p_vals).pvalue

                score = -1# math.log(p_same_treatment / p_inter_treatment, 10)

                tests_passed = 0
                if p_inter_treatment < (alpha / len(inter_treatment_p_vals)):
                    tests_passed += 1
                if p_same_treatment > (alpha / len(inter_treatment_p_vals)):
                    tests_passed += 1

            # print pandas tsv row summary
            # print("NUM READTHROUGHS: {}".format(d_read_through_counts))
            # print("NUM NORMAL: {}".format(d_normal_read_counts))

            row_summary = [row['ID'], d_wart_change[row['ID']], d_wart_before[row['ID']], d_wart_after[row['ID']], p_inter_treatment, p_same_treatment, readthrough_split_point, list(d_fitted_curve_r_squared.values()), list(d_fitted_curve_coeff.values()), average_expression, cannonical_mods_start_pos[row['ID']], d_wam_before[row['ID']], d_wam_after[row['ID']], d_wam_change[row['ID']]]
            summary_df.loc[summary_df_index] = row_summary
            summary_df_index += 1


    # TES_SUMMARY_PATH = "./tes_summary.tsv"
    print(summary_df)
    summary_df.to_csv(TES_SUMMARY_PATH, sep='\t', index=False)


    # --------- PLOT ---------- #
    if SINGLE_GENE_ANALYSIS and row['ID'] in d_tts_hist[first_label]:
        gene_id = matches.iloc[0]['ID']

        NUM_VERT_PLOTS = 3
        fig, axes = plt.subplots(NUM_VERT_PLOTS, num_bams)
        axes_index = 0
        for label in bam_labels:
            # scatter plot tts vs poly-a length
            axes[0, axes_index].scatter(d_tts[label][gene_id], d_poly_a_lengths[label][gene_id], s=1)
            axes[0, axes_index].set_ylim(ymin=0, ymax=d_max_poly_a[gene_id]*1.1)
            axes[0, axes_index].set_xlim(xmin=d_x_ticks[row['ID']][0], xmax=d_x_ticks[row['ID']][-1])
            axes[0, axes_index].get_xaxis().set_visible(False)
            
            sorted_tes_counts_by_pos = sorted(d_tts_hist[label][row['ID']], key=lambda a: a[1])
            d_tts_hist_y = [e[0] for e in sorted_tes_counts_by_pos]
            d_tts_hist_x = [e[1] for e in sorted_tes_counts_by_pos]
            axes[1, axes_index].plot(d_tts_hist_x, d_tts_hist_y)
            axes[1, axes_index].set_ylim(ymin=0, ymax=d_max_hist_count_tts[gene_id]*1.1)
            axes[1, axes_index].set_xlim(xmin=d_x_ticks[row['ID']][0], xmax=d_x_ticks[row['ID']][-1])
            axes[1, axes_index].get_xaxis().set_visible(False)

            # axes[2, axes_index].plot(d_poly_a_length_hists[label])
            # axes[2, axes_index].set_ylim(ymin=0, ymax=max_hist_count*1.1)
            # axes[2, axes_index].set_xlim(xmin=0, xmax=max_poly_a)
            # axes[2, axes_index].get_xaxis().set_visible(False)

            axes[2, axes_index].plot(d_x_ticks[row['ID']], d_kdes[label][gene_id])
            axes[2, axes_index].set_ylim(ymin=0, ymax=d_max_density[gene_id]*1.1)
            axes[2, axes_index].set_xlim(xmin=d_x_ticks[row['ID']][0], xmax=d_x_ticks[row['ID']][-1])
            axes[2, axes_index].set(xlabel='transcription end site (nt)')


            # PLOT GENE END AS VERT LINE
            axes[0, axes_index].axvline(x= gene_length, color='darkgray', ls="--", linewidth=1.0)
            axes[1, axes_index].axvline(x= gene_length, color='darkgray', ls="--", linewidth=1.0)
            axes[2, axes_index].axvline(x= gene_length, color='darkgray', ls="--", linewidth=1.0)

            # if SHOW_CANNONICAL_M6A:
            #     # # PLOT MOD LOCATION AS VERT LINES
            #     for mod_location in d_cannonical_mod_locations[row_name]:
            #         axes[0, axes_index].axvline(x= mod_location, color='red', ls="--", linewidth=1.0)
            #         axes[1, axes_index].axvline(x= mod_location, color='red', ls="--", linewidth=1.0)
            #         axes[2, axes_index].axvline(x= mod_location, color='red', ls="--", linewidth=1.0)

            # add axis labels
            if axes_index == 0:
                axes[0, axes_index].set(ylabel='poly-A length (nt)')
                axes[1, axes_index].set(ylabel='count')
                # axes[2, axes_index].set(xlabel='poly a length (nt)', ylabel='count')
                axes[2, axes_index].set(xlabel='transcription end site (nt)', ylabel='density (au)')
            else:
                axes[0, axes_index].get_yaxis().set_visible(False)
                axes[1, axes_index].get_yaxis().set_visible(False)
                # axes[2, axes_index].get_yaxis().set_visible(False)
                axes[2, axes_index].get_yaxis().set_visible(False)


            axes_index += 1

        fig.subplots_adjust(hspace=0, wspace=0.1)

        # from kneed import KneeLocator

        # NUM_VERT_PLOTS = 2
        # fig, axes = plt.subplots(NUM_VERT_PLOTS, num_bams)
        # axes_index = 0

        # for label in bam_labels:
        #     # scatter plot tts vs poly-a length
        #     tes = list(range(1, max(d_tts[label][gene_id]) + 1))
        #     paired_tes_hist = list(zip(d_tts_hist[label], tes))
        #     elbow = sorted([x for x in paired_tes_hist if x[0] > 0], key=lambda a: a[0], reverse=True)
        #     print(elbow)

        #     e1 = [x[0] for x in elbow]

        #     kneedle = KneeLocator(e1, list(range(len(e1))), S=1.0, curve='convex', direction='decreasing')
        #     cannonical_tes = elbow[0:kneedle.knee]
            
        #     max_cannonical_tes = 0
        #     for t in cannonical_tes:
        #         if t[1] > max_cannonical_tes:
        #             max_cannonical_tes = t[1]

        #     print(max_cannonical_tes)
        #     num_read_throughs = len([x for x in d_tts[label][gene_id] if x > max_cannonical_tes])
        #     num_normal = len([x for x in d_tts[label][gene_id] if x <= max_cannonical_tes])

        #     print("read throughs: {}, normal: {}".format(num_read_throughs, num_normal))

        #     axes[0, axes_index].plot(e1)

        #     axes_index += 1

            # axes[0, axes_index].set_ylim(ymin=0, ymax=max_poly_a*1.1)
            # axes[0, axes_index].set_xlim(xmin=min_tts, xmax=max_tts)
            # axes[0, axes_index].get_xaxis().set_visible(False)

        # TODO: also plot density violin plot of poly-A lengths. Violin plot of transcript termination sites
        # poly_a_labels = ["{}\nn={}".format(key, len(d_poly_a_lengths[key])) for key in d_poly_a_lengths.keys()]
        # fig, axes = plt.subplots(1, 1)
        # d_violins = axes.violinplot(d_poly_a_lengths.values())
        # axes.set_xticks(range(1,len(d_poly_a_lengths.keys())+1))
        # axes.set_xticklabels(poly_a_labels)
        # axes.set_ylabel('poly-A length (nt)')

        # tts_labels = ["{}\nn={}".format(key, len(d_tts[key])) for key in d_tts.keys()]
        # fig, axes = plt.subplots(1, 1)
        # axes.violinplot(d_tts.values())
        # axes.set_xticks(range(1,len(d_tts.keys())+1))
        # axes.set_xticklabels(tts_labels)
        # axes.set_ylabel('TTS (nt)')

        plt.show()

    else:
        print("plotting gene batch methylation changes and transcript end site changes...")

