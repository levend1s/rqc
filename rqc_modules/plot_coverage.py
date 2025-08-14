def plot_coverage(args):
    # load annotation file
    feature_id = INPUT[0]

    # process input file. Each line contains a label, the type of file, and the filepath
    input_files = {}
    # if len(INPUT[1:]) % 4 != 0:
    #     print("ERROR: not enough information in specified files. Check each input follows the format [LABEL] [TYPE] [PATH]")
    #     sys.exit()
    # else:
    in_index = 1
    while in_index < len(INPUT):
        if not INPUT[in_index].startswith("#"):
            input_files[INPUT[in_index]] = {
                'group': INPUT[in_index+1], 
                'type': INPUT[in_index+2], 
                'path': INPUT[in_index+3]
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
    matches = ANNOTATION_FILE.filter_feature_of_type([feature_id])
    if len(matches.df) == 0:
        evaluated_input = ast.literal_eval(feature_id)
        matches = ANNOTATION_FILE.get_feature_by_attribute("ID", evaluated_input)
        print("Looking for {} IDs, found {} matches. Plotting gene coverage for {}".format(len(evaluated_input) , len(matches.df), evaluated_input))
    else:
        print("Found {} matches for type {}. Plotting gene coverage...".format(len(matches.df), feature_id))

    num_matches = len(matches.df)
    PYSAM_PILEUP_MAX_DEPTH = 8000 # default
    subfeature_names = []
    subfeature_info = {}
    matches = matches.attributes_to_columns()
    matches['strand'] = matches['strand'].astype('category')
    matches['type'] = matches['type'].astype('category')
    matches['seq_id'] = matches['seq_id'].astype('category')


    index = 0
    sites_of_interest = None

    coverages = {}
    normalised_coverages = {}
    feature_coverages = {}
    normalised_feature_coverages = {}
    density_coverages = {}
    tx_lengths = {}
    mod_peaks = {}


    LOGFILE_PATH = "rqc_plot_coverage_stats.log"
    log_header = ["label", "type", "id", "max_depth", "total_coverage", "length", "AUC", "num peaks", "peak locations"]
    print("\t".join(log_header))
    logfile_df = pandas.DataFrame(columns=log_header)
    logfile_df_index = 0

    max_num_subfeatures = 0

    for label in input_files.keys():
        type = input_files[label]['type']
        path = input_files[label]['path']

        if type in ['bam', 'bedmethyl', 'bed']:

            feature_coverages[label] = {}
            normalised_feature_coverages[label] = {}

            mod_peaks[label] = {}

            if type == "bam":
                samfile = pysam.AlignmentFile(path, 'rb')
            elif type == "bedmethyl":
                mods_file_df = pandas.read_csv(path, sep='\t', names=MODKIT_BEDMETHYL_HEADER)
            elif type == "bed":
                site_file_df = pandas.read_csv(path, sep='\t', skiprows=1, names=GENERIC_BED_HEADER)
                site_file_df['strand'] = site_file_df['strand'].astype('category')
                site_file_df['contig'] = site_file_df['contig'].astype('category')

                # site_file_df['start'] = site_file_df['start'].astype('int32')
                # site_file_df['end'] = site_file_df['end'].astype('int32')

            # else:
            #     print("ERROR UNKNOWN FILE TYPE {}".format(type))

            # generate coverage for all matches in this bam file
            mal_annotation = []
            for row_index, row in matches.iterrows():
                # find subfeatures
                START_CLOCK("row_start")

                row_subfeatures = getSubfeatures(row['ID'], COVERAGE_TYPE, COVERAGE_PADDING)

                subfeature_names = row_subfeatures['type'].to_list()

                SKIP_MALANNOTATIONS = False
                # if SKIP_MALANNOTATIONS and ("five_prime_UTR" not in subfeature_names) or ("three_prime_UTR" not in subfeature_names) or \
                # (len(row_subfeatures[row_subfeatures.type == "five_prime_UTR"]) > 1) or \
                # (len(row_subfeatures[row_subfeatures.type == "three_prime_UTR"]) > 1):
                #     print("ERROR: {} has no or misannotated UTR's, skipping...".format(row.ID))
                #     mal_annotation.append(row.ID)
                #     continue

                # gen coverage for each subfeature in a gene
                subfeature_index = 0
                num_subfeatures = len(row_subfeatures.index)
                subfeature_base_coverages = [None] * num_subfeatures

                # if max_num_subfeatures == 0:
                #     max_num_subfeatures = num_subfeatures
                # else:
                #     if num_subfeatures != max_num_subfeatures:
                #         print("ERROR: trying to calculate subfeature coverage for genes with different number of subfeatures")
                #         print("{} has {} subfeatures, but previous genes had {} subfeatures. Exiting...".format(row['ID'], num_subfeatures, max_num_subfeatures))

                if num_subfeatures > 1:
                    if 'UTR' in subfeature_names[0]:
                        subfeature_names[0] = "5'UTR"
                    if 'UTR' in subfeature_names[1]:
                        subfeature_names[1] = "5'UTR"
                    if 'UTR' in subfeature_names[-1]:
                        subfeature_names[-1] = "3'UTR"
                    if 'UTR' in subfeature_names[-2]:
                        subfeature_names[-2] = "3'UTR"

                mod_peaks[label][row['ID']] = []

                exon_idx = 1
                for i in range(num_subfeatures):
                    if subfeature_names[i] == 'CDS':
                        subfeature_names[i] = "E{}".format(exon_idx)
                        exon_idx += 1

                tx_lengths[row['ID']] = (row['end'] - row['start'])

                if type == "bedmethyl":
                    row_mods_file_df = mods_file_df[
                        (mods_file_df.contig == row['seq_id']) & 
                        (mods_file_df.start > (row['start']) - COVERAGE_PADDING) & 
                        (mods_file_df.start < (row['end']) + COVERAGE_PADDING) & 
                        (mods_file_df.strand == row['strand'])
                    ]
                elif type == "bed":
                    if IGNORE_STRAND:
                        row_match_condition = (site_file_df.contig == row['seq_id']) & \
                            (site_file_df.start > (row['start']) - COVERAGE_PADDING) & \
                            (site_file_df.start < (row['end']) + COVERAGE_PADDING)
                    else:
                        row_match_condition = (site_file_df.contig == row['seq_id']) & \
                            (site_file_df.start > (row['start']) - COVERAGE_PADDING) & \
                            (site_file_df.start < (row['end']) + COVERAGE_PADDING) & \
                            (site_file_df.strand == row['strand'])
                    row_site_matches_df = site_file_df[
                        row_match_condition
                    ]

                row_flag_filters = BAM_PILEUP_DEFAULT_FLAGS
                row_flags_requires = 0

                if row['strand'] == '+':
                    row_flag_filters = BAM_PILEUP_DEFAULT_FLAGS | BAM_REVERSE_STRAND
                else:
                    row_flags_requires = BAM_REVERSE_STRAND

                for _, subfeature in row_subfeatures.iterrows():
                    subfeature_length = subfeature['end'] - subfeature['start'] + 1
                    subfeature_base_coverages[subfeature_index] = numpy.zeros(subfeature_length)

                    if type == "bam":
                        # pysam indexes are zero indexed but gff are 1-indexed, so pysam index = gffindex-1
                        for column in samfile.pileup(
                            contig=subfeature['seq_id'], 
                            start=subfeature['start'] - 1, 
                            stop=subfeature['end'],
                            # min_mapping_quality=MIN_MAPQ,
                            max_depth=PYSAM_PILEUP_MAX_DEPTH,
                            flag_require=row_flags_requires,
                            flag_filter=row_flag_filters,
                            truncate = True
                        ):
                            # take the number of aligned reads at this column position (column.n) minus the number of aligned reads which have either a skip or delete base at this column position (r.query_position)
                            read_depth = len(list(filter(None, column.get_query_sequences())))
                            # print(column.get_query_sequences(mark_matches=True, mark_ends=True, add_indels=True))
                            # reference pos is 0 indexed, gff (subfeature) is 1 indexed, add one to bring it back to zero
                            # TODO: this method of read depth shows only aligned bases. For reads which have mismatches/indels those bases do not contribute to read depth.
                            subfeature_base_coverages[subfeature_index][column.reference_pos - subfeature['start'] + 1] = read_depth

                    elif type == "bedmethyl":
                        mod_matches = row_mods_file_df[
                            (row_mods_file_df.end >= subfeature['start']) & 
                            (row_mods_file_df.end <= subfeature['end'])
                        ]
                        # convert from genome to transcript space
                        # in a bedmethyl file, the end position is the gff exact position
                        subfeature_mod_positions = mod_matches['end'].to_numpy() - subfeature['start']
                        num_mods_at_pos = mod_matches['num_mod'].to_list()

                        for mod_pos_index in range(len(num_mods_at_pos)):
                            subfeature_base_coverages[subfeature_index][subfeature_mod_positions[mod_pos_index]] = num_mods_at_pos[mod_pos_index]

                        # valid_cov and percent_mod determine 'cannonical mods'
                        if CANNONICAL_MOD_PROP_THRESHOLD > 0 and CANNONICAL_MOD_READ_DEPTH_THRESHOLD > 0:
                            mod_peak_matches = mod_matches[
                                (mod_matches.percent_mod >= (CANNONICAL_MOD_PROP_THRESHOLD * 100)) & 
                                (mod_matches.valid_cov >= CANNONICAL_MOD_READ_DEPTH_THRESHOLD)
                            ]

                            peak_positions = mod_peak_matches['end'].to_numpy() - row['start']

                            if (row["strand"] == "-"):
                                peak_positions = tx_lengths[row['ID']] - peak_positions

                            mod_peaks[label][row['ID']] += peak_positions.tolist()

                    elif type == "bed":
                        site_matches = row_site_matches_df[
                            ((row_site_matches_df.start - 3) >= subfeature['start']) & 
                            ((row_site_matches_df.start + 3) <= subfeature['end'])
                        ]
                        # convert from genome to transcript space
                        # start position is the 0-indexed start position of the 5mer, so add 3 to get the gff exact position of the central A to DRACH motif. 
                        subfeature_site_positions = site_matches['start'].to_numpy() - subfeature['start'] + 3

                        for site_pos_index in subfeature_site_positions:
                            subfeature_base_coverages[subfeature_index][site_pos_index] = 1
                    # else:
                    #     print("WARNING: unknown type: {}".format(type))


                    subfeature_index += 1

                total_count_depth_gene = 0
                # squish the cds subfeatures into a single one
                if COVERAGE_TYPE == "subfeature_cds":
                    num_exons = 0
                    first_exon_idx = -1
                    last_exon_idx = -1

                    temp_subfeature_names = []
                    for i in range(len(subfeature_names)):
                        total_count_depth_gene += sum(subfeature_base_coverages[i])
                        if subfeature_names[i].startswith("E"):
                            if first_exon_idx == -1:
                                temp_subfeature_names.append("CDS")
                                first_exon_idx = i
                            
                            last_exon_idx = i
                            num_exons += 1
                        else:
                            temp_subfeature_names.append(subfeature_names[i])
                    
                    num_subfeatures = num_subfeatures - num_exons + 1
                    temp_subfeature_base_coverages = [None] * num_subfeatures

                    offset = 0
                    for i in range(num_subfeatures):
                        if i == first_exon_idx:
                            temp_subfeature_base_coverages[i] = numpy.concatenate(subfeature_base_coverages[first_exon_idx:last_exon_idx+1]).ravel()
                            offset = num_exons-1
                        else:
                            temp_subfeature_base_coverages[i] = subfeature_base_coverages[i+offset]

                    subfeature_base_coverages = temp_subfeature_base_coverages
                    corrected_subfeature_names = temp_subfeature_names
                else:
                    corrected_subfeature_names = subfeature_names


                sf_base_coverage_list = [None] * num_subfeatures
                STOP_CLOCK("row_start", "coverage_stop")

                if COVERAGE_PADDING:
                    num_bins_cds = int(COVERAGE_BINS * (1 - (2 * PADDING_RATIO)))
                    num_bins_padding = int(COVERAGE_BINS * PADDING_RATIO)
                else:
                    num_bins_cds = COVERAGE_BINS
                    num_bins_padding = 0

                running_sf_bin_count = 0

                # implicitly reset subfeature_index
                for subfeature_index in range(num_subfeatures):
                    # resample coverage
                    if subfeature_index == 0 and COVERAGE_PADDING:
                        sf_bin_size = num_bins_padding
                    elif subfeature_index == (num_subfeatures - 1):
                        # sf_bin_size = num_bins_cds - (math.floor(num_bins_cds / num_subfeatures) * (num_subfeatures-1))
                        sf_bin_size = COVERAGE_BINS - running_sf_bin_count
                    elif COVERAGE_PADDING:
                        sf_bin_size = math.floor(num_bins_cds / (num_subfeatures - 2))
                    else:
                        sf_bin_size = math.floor(num_bins_cds / num_subfeatures)
                    
                    subfeature_info[corrected_subfeature_names[subfeature_index]] = sf_bin_size

                    # if type == "bed":
                    #     sf_resampled_coverage = resample_coverage(subfeature_base_coverages[subfeature_index], sf_bin_size, "sum")
                    # else:
                    sf_resampled_coverage = resample_coverage(subfeature_base_coverages[subfeature_index], sf_bin_size, COVERAGE_METHOD)

                    sf_base_coverage_list[subfeature_index] = sf_resampled_coverage

                    running_sf_bin_count += sf_bin_size

                # flatten resampled subfeature coverages into a single array
                resampled_base_coverage = numpy.concatenate(sf_base_coverage_list).ravel()
                # reverse coverages if necessary
                if (row["strand"] == "-"):
                    resampled_base_coverage = numpy.flip(resampled_base_coverage)

                # # find out how many mod peaks there are based off thresholds
                # if MOD_PROP_THRESHOLD > 0 and READ_DEPTH_THRESHOLD > 0:
                #     num_prop_threshold_peaks = 0
                #     for i in range(len(feature_coverages[read_cov_label][feature_index])):
                #         if feature_coverages[read_cov_label][feature_index][i] >= READ_DEPTH_THRESHOLD and mod_base_proportion[i] >= MOD_PROP_THRESHOLD:
                #             num_prop_threshold_peaks += 1

                #     additional_info += "\tmod peaks: {}".format(num_prop_threshold_peaks)

                STOP_CLOCK("row_start", "resample_subfeature_stop")
                additional_info = [len(mod_peaks[label][row['ID']]), ",".join(str(x) for x in mod_peaks[label][row['ID']])]


                # if this is modification coverage, we'll 'normalise' it against the gene read depth coverage
                if type == "bedmethyl":
                    read_cov_label = label.split("_")[0] + "_read_depth"

                    # weighted probability of m6A function
                    # for each given site, we have P(m6A) = num_m6A / read_depth
                    # P(m6A) prior = 0.05, which is the abundance of m6A / A in entire RNA-seq
                    # formula for weighted probability is P_weighted = (N * P_observed) + (Weight_prior * P_prior) / (N + Weight_prior)
                    # Weight_prior = feature_coverages[read_cov_label][feature_index].max()
                    # P_prior = 0.01

                    # resampled_base_coverage = feature_coverages[read_cov_label][feature_index] / 2
                    # denom = (feature_coverages[read_cov_label][feature_index] + Weight_prior)
                    # normalised_feature_coverages[feature_index] = numpy.nan_to_num( (resampled_base_coverage + (Weight_prior * P_prior)) / denom)

                    # normalise against itself
                    #normalised_feature_coverages[feature_index] = normalise_coverage(resampled_base_coverage)

                    # normalise against read depth (fraction of bases methylated * normalised coverage)
                    mod_base_proportion = numpy.nan_to_num(resampled_base_coverage / feature_coverages[read_cov_label][row['ID']])

                    # NOTE: this is due to how we calculate read depth, mentioned in the bam section above
                    # There are cases where a read may have indels/mismatches (which do not contribute to read depth) but within those sections a mod is detected
                    # this leads to numbers greater than 1 when calculating mod proportion
                    # So for now we'll just clamp those numbers down to 1
                    mod_base_proportion[mod_base_proportion > 1.0] = 1.0

                    if MOD_NORMALISATION == "raw":
                        normalised_feature_coverages[label][row['ID']] = mod_base_proportion
                    else:
                        normalised_feature_coverages[label][row['ID']] = mod_base_proportion * normalised_feature_coverages[read_cov_label][row['ID']]
                else:
                    normalised_feature_coverages[label][row['ID']] = normalise_coverage(resampled_base_coverage)

                feature_coverages[label][row['ID']] = resampled_base_coverage
                AUC = round(numpy.sum(normalised_feature_coverages[label][row['ID']]) / COVERAGE_BINS, 2) # gives score between 0 and 1

                STOP_CLOCK("row_start", "row_end")

                # label, gene id, max coverage, gene length, auc, num mod peaks, mod peaks
                row_coverage_summary = [label, type, row['ID'], int(max(resampled_base_coverage)), total_count_depth_gene, row['end'] - row['start'], AUC]

                if additional_info:
                    row_coverage_summary += additional_info

                print("\t".join([str(x) for x in row_coverage_summary]))
                logfile_df.loc[logfile_df_index] = row_coverage_summary
                logfile_df_index += 1

            if type == "bam":
                samfile.close()

    # drop all coverages which don't meet a coverage threshold across ALL samples
    # this could be AUC or read depth

    # for each gene, get how many rows 
    if READ_DEPTH_THRESHOLD > 0:
        bam_labels = []
        for label in input_files.keys():
            type = input_files[label]['type']

            if type == "bam":
                bam_labels.append(label)

        num_low_coverage = 0

        for row_index, row in matches.iterrows():
            gene_matches_below_read_threshold = logfile_df[
                (logfile_df.id == row['ID']) & 
                (logfile_df.max_depth < READ_DEPTH_THRESHOLD) &
                (logfile_df.type == "bam")]
            
            if len(gene_matches_below_read_threshold) > 0:
                # print("ignoring {} since it has low read depth (<{})...".format(row['ID'], READ_DEPTH_THRESHOLD))
                # print(gene_matches_below_read_threshold)

                for label, file in input_files.items():
                    # print(feature_coverages[label][row['ID']])
                    feature_coverages[label].pop(row['ID'], None)
                    normalised_feature_coverages[label].pop(row['ID'], None)

                num_low_coverage += 1 

        print("REMOVED {} DUE TO LOW COVERAGE (<{})".format(num_low_coverage, READ_DEPTH_THRESHOLD))
        print("REMAINING ID's: {}".format(feature_coverages[label].keys()))

    print("REMOVED {} DUE TO MISSING UTR annotation".format(len(mal_annotation)))

    for label in input_files.keys():
        label_total_depth = logfile_df[logfile_df['label'] == label].total_coverage.sum()
        print("{} TOTAL DEPTH: {}".format(label, label_total_depth))

    x_ticks = list(range(COVERAGE_BINS))

    for label in input_files.keys():
        type = input_files[label]['type']
        if type in ['bam', 'bedmethyl', 'bed']:

            # flatten down all resampled coverages for this label and store dict under label
            total_coverage = numpy.array([sum(i) for i in zip(*list(feature_coverages[label].values()))], dtype=numpy.uint32)
            all_normalised_total_coverage = numpy.array([sum(i) for i in zip(*list(normalised_feature_coverages[label].values()))])

            # normalised mod coverage is the average weighted proportion of a modification against read depth across all genes
            if type == "bedmethyl":
                normalised_total_coverage = all_normalised_total_coverage / len(normalised_feature_coverages[label].keys())
            else:
                normalised_total_coverage = normalise_coverage(all_normalised_total_coverage)

            if PLOT_DENSITY:
                base_idx_counts = []
                for i in range(COVERAGE_BINS):
                    if total_coverage[i] > 0:
                        for j in range(total_coverage[i]):
                            base_idx_counts.append(i)

                kernel = scipy.stats.gaussian_kde(base_idx_counts, bw_method=0.1)
                # HACK PAM
                # this is not a true KDE anymore now that we're multiplying it by total coverage
                # it does however allow us to plot and compare multiple KDEs on the same plot
                smoothed_tts_hist = kernel(x_ticks) * sum(total_coverage)
                density_coverages[label] = smoothed_tts_hist

            # cdf = numpy.cumsum(smoothed_tts_hist)

            # if type == "bed":
            #     sites_of_interest = total_coverage # normalised_total_coverage
            # else:
            coverages[label] = total_coverage
            normalised_coverages[label] = normalised_total_coverage


    # print("\nsummary:\nnum matches: {}\nnum bins: {}".format(num_matches, COVERAGE_BINS))
    additional_text = "num transcripts: {}\naverage transcript length: {}".format(len(tx_lengths.keys()), int(sum(tx_lengths.values()) / len(tx_lengths.values())))
    logfile_df.to_csv(LOGFILE_PATH, sep='\t', index=False)

    # plot coverages
    # this looks at coverage for each gene, resamples and normalises the coverage and adds it to a list
    # then takes the average of all those resampled and normalised coverages
    # this smooths out the cases where some genes might have read depth in the 1000's, and others in the 10's
    # so our data isn't skewed toward genes that are higher expressed
    # coverage: dict of read depths of transcript: eg {'sample1': [coverage...], 'sample2': [coverage...]}
    # mod_coverage: dict of mod coverages to plot: eg {'sample1m6As': [coverage...], 'sample2m6as': [coverage...], 'sample1pseU': [coverage...]}
    # sites_of_interest: dict of sites of interest to plot as vertical lines??? But how to do this for aggregate transcript searches?
    #       Maybe plot vlines for individual transcript plots and areas shaded with intensity according to how often motifs appear
    coverage_dict = {
        "coverages": coverages,
        "method": COVERAGE_METHOD,
        "num_matches": num_matches,
        "sites_of_interest": sites_of_interest,
        "num_bins": COVERAGE_BINS,
        "subfeature_names": corrected_subfeature_names,
        "subfeature_info": subfeature_info,
        "y_label": "count (nt)",
        "additional_text": additional_text
    }

    plot_subfeature_coverage(coverage_dict)

    coverage_dict['coverages'] = normalised_coverages
    coverage_dict['y_label'] = "normalised coverage (au)"

    plot_subfeature_coverage(coverage_dict)

    if PLOT_DENSITY:
        coverage_dict['coverages'] = density_coverages
        coverage_dict['y_label'] = "density (au)"

        plot_subfeature_coverage(coverage_dict)

    plt.tight_layout()
    plt.show()

    if (OUTFILE):
        plt.savefig("coverage_{}".format(OUTFILE))
