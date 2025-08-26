import pandas
import pysam
import numpy
import math
import scipy.stats
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel, wilcoxon
import statsmodels.formula.api as smf
import statsmodels.api as sm

NUM_DPS = 4

from rqc_modules.utils import process_annotation_file, process_genome_file, normalise_numpy_array, power_func, reverse_complement, process_input_files
from rqc_modules.utils import START_CLOCK, STOP_CLOCK, getSubfeatures

from rqc_modules.constants import TES_SUMMARY_HEADER, FEATURECOUNTS_HEADER, UNIQUE_APA_DISTANCE

def process_row():
    print("process row")

def find_num_read_throughs(strand, tx_end_sites, split_point, tolerance):
    if strand == "-":
        num_read_throughs = len([x for x in tx_end_sites if x < (split_point - tolerance)])
    else:
        num_read_throughs = len([x for x in tx_end_sites if x > (split_point + tolerance)])

    return num_read_throughs

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
    COMPARE_APA_BETWEEN_TREATMENTS = args.compare_apa_between_treatments

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

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']
    summary_df = pandas.DataFrame(columns=TES_SUMMARY_HEADER)

    gene_length = 0

    row_header = [
        "label",
        "ID",
        "strand",
        "reads(mapped)",
        "reads_in_region",
        "filtered(strand)",
        "filtered(3p)"
    ]
    print("\t".join(map(str, row_header)))

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
        for _, row in matches.iterrows():
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
            # row_subfeatures = getSubfeatures(gff_df, row['ID'], "subfeature", 0)

            read_indexes_to_process = []

            read_outside_3p_end = []
            read_on_different_strand = []

            # example: PF3D7_0709050.1
            # if len(row_subfeatures) == 0:
            #     most_3p_subfeature = row
            # else:
            #     if row['strand'] == "-":
            #         most_3p_subfeature = row_subfeatures.iloc[0]
            #     else:
            #         most_3p_subfeature = row_subfeatures.iloc[-1]

            # NOTE: now the longest function in the TES analysis
            for i in range(len(reads_in_region)):
                r = reads_in_region[i]

                # keep only reads in the same direction as this strand
                if (row['strand'] == "+" and r.is_forward) or (row['strand'] == "-" and r.is_reverse):
                    if row['strand'] == "-":
                        read_3p_end = r.reference_start

                        if read_3p_end <= row.end:
                            read_indexes_to_process.append(i)
                        else:
                            read_outside_3p_end.append(r.query_name)

                    else:
                        read_3p_end = r.reference_end

                        if read_3p_end >= row.start:
                            read_indexes_to_process.append(i)
                        else:
                            read_outside_3p_end.append(r.query_name)
                else:
                    read_on_different_strand.append(r.query_name)

            no_poly_a = 0
            poly_a_lengths = []
            tts_sites = []

            for i in read_indexes_to_process:
                r = reads_in_region[i]
            
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


            row_summary = [label, row.ID, row.strand, len(reads_in_region), len(tts_sites), len(read_on_different_strand), len(read_outside_3p_end)]
            print("\t".join(map(str, row_summary)))

        samfile.close()



    # can calculate the log ratio
    d_x_ticks = {}
    d_poly_a_length_hists = {}
    d_tts_hist = {}
    d_kdes = {}
    d_max_hist_count_poly_a = {}
    d_max_hist_count_tts = {}
    d_max_poly_a = {}
    d_min_tts = {}
    d_max_tts = {}
    d_max_density = {}
    d_genomic_apa_sites = {}

    raw_summary_header = [
        "label",
        "ID",
        "type",
        "strand",
        "3p_annotation end",
        "tts_count",
        "canonical_pa_site",
        "apa_score",
        "genomic_apa_sites",
        "pa_site_counts",
        "pa_site_proportions"
    ]
    print("\t".join(map(str, raw_summary_header)))
    raw_results = []

    for label in bam_labels:
        d_poly_a_length_hists[label] = {}
        d_tts_hist[label] = {}
        d_kdes[label] = {}
        d_genomic_apa_sites[label] = {}

    for row_index, row in matches.iterrows():
        if row['strand'] == "+":
            annotation_row_3p_end = row['end']
        else:   
            annotation_row_3p_end = row['start']

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
        # print("{} - Generating transcript end site histograms...".format(row['ID']))
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
            # print("{} - Generating transcript end site density information...".format(row['ID']))
            for label in bam_labels:
                kernel = scipy.stats.gaussian_kde(d_tts[label][gene_id], bw_method=0.1)
                smoothed_tts_hist = kernel(x_ticks)

                d_kdes[label][row['ID']] = smoothed_tts_hist

                if max_density < max(smoothed_tts_hist):
                    max_density = max(smoothed_tts_hist)

        d_max_hist_count_poly_a[gene_id] = max_hist_count_poly_a
        d_max_hist_count_tts[gene_id] = max_hist_count_tts
        d_max_poly_a[gene_id] = max_poly_a
        d_min_tts[gene_id] = min_tts
        d_max_tts[gene_id] = max_tts
        d_max_density[gene_id] = max_density
        d_x_ticks[row['ID']] = x_ticks

        for label in bam_labels:
            tts_hist = [d_tts[label][gene_id].count(i) for i in range(min_tts, max_tts)]
            peaks, peak_dict = find_peaks(tts_hist, distance=UNIQUE_APA_DISTANCE, height=READ_DEPTH_THRESHOLD)
            
            genomic_apa_sites = peaks + min_tts
            max_height_apa = 0
            canonical_pa_site = 0
            second_max_height_apa = 0

            if len(peaks) > 1:
                sorted_genomic_peaks_by_height = sorted(zip(genomic_apa_sites, peak_dict['peak_heights']), key=lambda x: x[1])
                max_height_apa = sorted_genomic_peaks_by_height[-1][1]
                canonical_pa_site = sorted_genomic_peaks_by_height[-1][0]

                second_max_height_apa = sorted_genomic_peaks_by_height[-2][1]
                apa_score = second_max_height_apa / max_height_apa
            elif len(peaks) == 1:
                canonical_pa_site = genomic_apa_sites[0]
                apa_score = 0
            else:
                canonical_pa_site = 0
                apa_score = 0

            pa_site_proportions = []
            POLY_ADENYLATION_TOLERANCE = 10

            for p in genomic_apa_sites:
                rt_prop = find_num_read_throughs(row.strand, d_tts[label][gene_id], p, POLY_ADENYLATION_TOLERANCE) / len(d_tts[label][gene_id])
                
                # print("p: {}".format(p))

                # print("num_read_throughs: {}".format(num_read_throughs))
                # print("num_normal: {}".format(num_normal))
                #rt_prop = num_read_throughs / num_normals

                pa_site_proportions.append(rt_prop)

            num_tts = len(d_tts[label][gene_id])
            pa_site_counts = list(peak_dict['peak_heights'])
            d_genomic_apa_sites[label][row.ID] = genomic_apa_sites

            # TODO add max_read_depth
            row_summary = [label, row.ID, row.type, row.strand, annotation_row_3p_end, num_tts, canonical_pa_site, round(apa_score, NUM_DPS), genomic_apa_sites, [int(x) for x in pa_site_counts], [round(float(x), NUM_DPS) for x in pa_site_proportions]]
            print("\t".join(map(str, row_summary)))
            raw_results.append(row_summary)

    raw_summary_df = pandas.DataFrame(raw_results, columns=raw_summary_header)

    if OUTFILE:
        raw_summary_df.to_csv("raw_{}".format(OUTFILE), sep='\t', index=False)

    # go through d_kdes, find all local max's with count > read_depth threshold and call these poly_adenylation sites
    # The max PA is the canonical poly_adenylation site, and belongs in it's own column
    if COMPARE_APA_BETWEEN_TREATMENTS:
        summary_header = [
            "ID",
            "type",
            "strand",
            "3p_annotation_end",
            "canonical_pa_site",
            "average_rt_prop_g1",
            "average_rt_prop_g2",
            "t_stat",
            "p_val"
        ]
        print("\t".join(map(str, summary_header)))
        results = []
        group_1_prefix = COMPARE_APA_BETWEEN_TREATMENTS[0]
        group_2_prefix = COMPARE_APA_BETWEEN_TREATMENTS[1]

        group1_rows = raw_summary_df[raw_summary_df['label'].str.startswith(group_1_prefix)]
        group2_rows = raw_summary_df[raw_summary_df['label'].str.startswith(group_2_prefix)]

        group1_bam_labels = group1_rows.label.to_list()
        group2_bam_labels = group2_rows.label.to_list()

        # go through each row, determine quantify change in 
        for row_index, row in matches.iterrows():
            gene_id = row.ID
            if row['strand'] == "+":
                annotation_row_3p_end = row['end']
            else:   
                annotation_row_3p_end = row['start']
            these_rows_g1 = group1_rows[group1_rows['ID'] == row.ID]
            # these_rows_g2 = group2_rows[group2_rows['ID'] == row.ID]

            # look through both treatments and keep only apa sites common to all of them
            canonical_pas = these_rows_g1.canonical_pa_site.to_list()
            canonical_pas = [x for x in canonical_pas if x > 0]
            average_canonical_pa = int(sum(canonical_pas) / len(canonical_pas))

            group1_num_read_throughs = []
            group2_num_read_throughs = []

            group1_read_through_props = []
            group2_read_through_props = []

            for label in group1_bam_labels:
                num_read_throughs = find_num_read_throughs(row.strand, d_tts[label][gene_id], average_canonical_pa, POLY_ADENYLATION_TOLERANCE)
                group1_num_read_throughs.append(num_read_throughs)
                group1_read_through_props.append(num_read_throughs / len(d_tts[label][gene_id]))
            for label in group2_bam_labels:
                num_read_throughs = find_num_read_throughs(row.strand, d_tts[label][gene_id], average_canonical_pa, POLY_ADENYLATION_TOLERANCE)
                group2_num_read_throughs.append(find_num_read_throughs(row.strand, d_tts[label][gene_id], average_canonical_pa, POLY_ADENYLATION_TOLERANCE))
                group2_read_through_props.append(num_read_throughs / len(d_tts[label][gene_id]))

            group1_average_rt_prop = sum(group1_read_through_props) / len(group1_read_through_props)
            group2_average_rt_prop = sum(group2_read_through_props) / len(group2_read_through_props)

            num_tes = [len(d_tts[l][gene_id]) for l in bam_labels]

            if len(canonical_pas) > 0:
                # from chatGPT
                data = pandas.DataFrame({
                    'group': [group_1_prefix] * len(group1_bam_labels) + [group_2_prefix] * len(group2_bam_labels),
                    'num_read_throughs': group1_num_read_throughs + group2_num_read_throughs,
                    'num_tes': num_tes
                })

                data['group_binary'] = (data['group'] == group_2_prefix).astype(int)

                # Failures column
                data['num_normal'] = data['num_tes'] - data['num_read_throughs']
                response = data[['num_read_throughs', 'num_normal']]

                # Fit Binomial GLM
                model = sm.GLM(response, sm.add_constant(data['group_binary']), family=sm.families.Binomial())
                result = model.fit()
                # print(result.summary())
                test_stat = result.params.group_binary
                pval = result.pvalues.group_binary
            else:
                test_stat = 0
                pval = 0

            row_summary = [row.ID, row.type, row.strand, annotation_row_3p_end, average_canonical_pa, round(group1_average_rt_prop, NUM_DPS), round(group2_average_rt_prop, NUM_DPS), test_stat, pval]
            results.append(row_summary)
            print("\t".join(map(str, row_summary)))

        summary_df = pandas.DataFrame(results, columns=summary_header)

        if OUTFILE:
            summary_df.to_csv(OUTFILE, sep='\t', index=False)


    # --------- PLOT ---------- #
    if len(matches) == 1:
        row = matches.iloc[0]
        gene_id = row.ID
        if row.strand == "-":
            gene_end = row.start
        else:
            gene_end = row.end

        NUM_VERT_PLOTS = 3
        fig, axes = plt.subplots(NUM_VERT_PLOTS, num_bams)
        axes_index = 0
        for label in bam_labels:
            # scatter plot tts vs poly-a length
            if len(bam_labels) == 1:
                poly_a_plot = axes[0]
                tts_count_plot = axes[1]
                tts_density_plot = axes[2]
            else:
                poly_a_plot = axes[0, axes_index]
                tts_count_plot = axes[1, axes_index]
                tts_density_plot = axes[2, axes_index]

            poly_a_plot.scatter(d_tts[label][gene_id], d_poly_a_lengths[label][gene_id], s=1)
            poly_a_plot.set_ylim(ymin=0, ymax=d_max_poly_a[gene_id]*1.1)
            poly_a_plot.set_xlim(xmin=d_x_ticks[row['ID']][0], xmax=d_x_ticks[row['ID']][-1])
            poly_a_plot.get_xaxis().set_visible(False)
            
            sorted_tes_counts_by_pos = sorted(d_tts_hist[label][row['ID']], key=lambda a: a[1])
            d_tts_hist_y = [e[0] for e in sorted_tes_counts_by_pos]
            d_tts_hist_x = [e[1] for e in sorted_tes_counts_by_pos]
            tts_count_plot.plot(d_tts_hist_x, d_tts_hist_y)
            tts_count_plot.set_ylim(ymin=0, ymax=d_max_hist_count_tts[gene_id]*1.1)
            tts_count_plot.set_xlim(xmin=d_x_ticks[row['ID']][0], xmax=d_x_ticks[row['ID']][-1])
            tts_count_plot.get_xaxis().set_visible(False)

            tts_density_plot.plot(d_x_ticks[row['ID']], d_kdes[label][gene_id])
            tts_density_plot.set_ylim(ymin=0, ymax=d_max_density[gene_id]*1.1)
            tts_density_plot.set_xlim(xmin=d_x_ticks[row['ID']][0], xmax=d_x_ticks[row['ID']][-1])
            tts_density_plot.set(xlabel='transcription end site (nt)')

            # PLOT GENE END AS VERT LINE
            poly_a_plot.axvline(x= gene_length, color='darkgray', ls="--", linewidth=1.0)
            tts_count_plot.axvline(x= gene_length, color='darkgray', ls="--", linewidth=1.0)
            tts_density_plot.axvline(x= gene_length, color='darkgray', ls="--", linewidth=1.0)

            for apa in d_genomic_apa_sites[label][gene_id]:
                tts_density_plot.axvline(x= apa, color='red', ls="--", linewidth=1.0)

            tts_density_plot.axvline(x= gene_end, color='grey', ls="--", linewidth=1.0)


            # add axis labels
            if axes_index == 0:
                poly_a_plot.set(ylabel='poly-A length (nt)')
                tts_count_plot.set(ylabel='count')
                tts_density_plot.set(xlabel='tes (nt)', ylabel='density (au)')
            else:
                poly_a_plot.get_yaxis().set_visible(False)
                tts_count_plot.get_yaxis().set_visible(False)
                tts_density_plot.get_yaxis().set_visible(False)

            axes_index += 1

        fig.subplots_adjust(hspace=0, wspace=0.1)
        plt.show()

