import pandas
import pysam
import numpy
import itertools
import math
import scipy.stats
import scipy.interpolate
import scipy.optimize
import matplotlib.pyplot as plt
from kneed import KneeLocator
from statsmodels.stats.proportion import proportions_ztest

from rqc_modules.utils import process_annotation_file, process_genome_file, normalise_numpy_array, power_func, reverse_complement, process_input_files
from rqc_modules.utils import START_CLOCK, STOP_CLOCK, getSubfeatures

from rqc_modules.constants import TES_SUMMARY_HEADER, FEATURECOUNTS_HEADER

def approximate_tes(args):
    # load annotation file
    ANNOTATION_FILE = args.annotation
    INPUT = args.input
    OUTFILE = args.output
    IDS = args.ids
    FEATURE_TYPE = args.type
    FILTER_FOR_M6A = args.filter_for_m6A
    FILTER_OUT_M6A = args.filter_out_m6A
    VERBOSE = args.verbose
    FILTER_FOR_FEATURE_COUNTS = args.feature_counts
    PADDING = args.padding
    READ_DEPTH_THRESHOLD = args.read_depth

    # process input file. Each line contains a label, the type of file, and the filepath
    input_files = process_input_files(INPUT)

    # load annotation file and find indexes for all parent children
    gff_df = process_annotation_file(ANNOTATION_FILE)

    print(input_files)

    if FEATURE_TYPE:
        matches = gff_df[gff_df['type'] == FEATURE_TYPE]
    else:
        matches = gff_df[gff_df['ID'].isin(IDS)]
    if matches.empty:
        print("ERROR: no matches found for type {} and ids {}".format(FEATURE_TYPE, IDS))

    num_bams = 0

    d_poly_a_lengths = {}
    d_tts = {}
    d_mod_info = {}
    
    d_not_beyond_3p = {}
    d_not_in_feature_counts = {}

    logfile_df_index = 0

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']
    summary_df = pandas.DataFrame(columns=TES_SUMMARY_HEADER)

    gene_length = 0

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
        if FILTER_FOR_FEATURE_COUNTS:
            feature_counts_sample_label = label.split("_")[0] + "_featureCounts"
            feature_counts_df = pandas.read_csv(input_files[feature_counts_sample_label]['path'], sep='\t', names=FEATURECOUNTS_HEADER)
            feature_counts_df['targets'] = feature_counts_df['targets'].astype('category')

        # generate coverage for all matches in this bam file
        for row_index, row in matches.iterrows():
            # find subfeatures
            # 0.8235001564025879s - 1.4s
            # START_CLOCK("row_start")

            summary_df_index = 0
            read_on_different_strand = 0

            # 0.30355286598205566s
            # START_CLOCK("fetch")

            # 0.0003178119659423828s
            reads_in_region = samfile.fetch(
                contig=row['seq_id'], 
                start=row['start'] - PADDING, 
                stop=row['end'] + PADDING
            )
            reads_in_region = list(reads_in_region)

            gene_length = row['end'] - row['start']
            row_name = row['ID']

            # !!!!! START NANOPORE SPECIFIC !!!!!
            # filter out reads where the 3' end is not in or beyond the last feature (3'UTR or last exon) of the target gene
            row_subfeatures = getSubfeatures(gff_df, row['ID'], "subfeature", 0)

            read_indexes_to_process = []

            missing_cannonical_mods = []
            read_outside_3p_end = []

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
                            read_indexes_to_process.append(i)
                        else:
                            read_outside_3p_end.append(r.query_name)

                    else:
                        read_3p_end = r.reference_end

                        if read_3p_end >= most_3p_subfeature.start:
                            # read_indexes_to_process.append(this_index)
                            
                                read_indexes_to_process.append(i)
                        else:
                            read_outside_3p_end.append(r.query_name)
                else:
                    read_on_different_strand += 1

            if FILTER_FOR_FEATURE_COUNTS:
                gene_reads = feature_counts_df[feature_counts_df.targets == row['ID'].split(".")[0]]
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



    # can calculate the log ratio
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

    

    summary_df_index = 0
    for row_index, row in matches.iterrows():
        SAMPLE_HAS_LOW_EXP = False
        average_expression = 0
        average_not_beyond_3p = 0
        average_not_in_feature_counts = 0

        if row['strand'] == "+":
            annotation_row_3p_end = row['end']
        else:   
            annotation_row_3p_end = row['start']

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
            row_summary = [row['ID'], annotation_row_3p_end, 0, average_expression, 0]
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
        if len(matches) == 1:
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
            if VERBOSE:
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
        for label in bam_labels:
            readthrough_split_point += d_readthrough_split_points[label]

        readthrough_split_point = int(readthrough_split_point / len(bam_labels))

        if VERBOSE:
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
        weighted_rt_ratios = []
        d_wart = {}

        # calculate sum of total cannonical mods read depth across all samples for weighting
        print("{} - Calculating weighted proportion change in readthroughs...".format(row['ID']))
        valid_cov_total = 0
        weighted_outputs = []
        for label in bam_labels:
            valid_cov_total += d_normal_read_counts[label] + d_read_through_counts[label]

        for label in bam_labels:
            weight = (d_normal_read_counts[label] + d_read_through_counts[label]) / valid_cov_total

            this_weighted_mod_proportion = (d_read_through_counts[label] / d_normal_read_counts[label]) * weight
            weighted_outputs.append(this_weighted_mod_proportion)

        rt = sum(weighted_outputs)
        d_wart[row['ID']] = rt

        # FIXME: row end should be start for negative strand?
        row_summary = [row['ID'], annotation_row_3p_end, readthrough_split_point, average_expression, d_wart[row['ID']]]
        summary_df.loc[summary_df_index] = row_summary
        summary_df_index += 1

    print(summary_df)

    if OUTFILE:
        summary_df.to_csv(OUTFILE, sep='\t', index=False)

    # --------- PLOT ---------- #
    if len(matches) == 1:
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

