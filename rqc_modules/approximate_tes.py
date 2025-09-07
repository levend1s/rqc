import pandas
import pysam
import numpy
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel, wilcoxon, kruskal, gaussian_kde
import statsmodels.api as sm
import os

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
    EXCLUDE_CONTIGS = args.exclude_contigs
    POLY_A_FILTER = args.poly_a_filter

    if OUTFILE:
        output_dir = os.path.dirname(OUTFILE)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        OUTPUT_COMPARE = OUTFILE.replace(".tsv", "_compare.tsv")
        OUTPUT_RAW = OUTFILE.replace(".tsv", "_raw.tsv")

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

    num_bams = 0

    d_poly_a_lengths = {}
    d_tts = {}
    d_mod_info = {}
    d_utr_lengths = {}
    
    d_not_beyond_3p = {}

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
        "filtered(3p)",
        "filtered(poly_a)"
    ]
    print("\t".join(map(str, row_header)))

    for label in bam_labels:
        samfile = pysam.AlignmentFile(input_files[label]['path'], 'rb')
        num_bams += 1
        d_poly_a_lengths[label] = {}
        d_tts[label] = {}
        d_mod_info[label] = {}
        d_utr_lengths[label] = {}
        
        d_not_beyond_3p[label] = {}

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

            read_indexes_to_process = []

            read_outside_3p_end = []
            read_on_different_strand = []

            row_subfeatures = getSubfeatures(gff_df, row['ID'], "subfeature", 0)
            row_subfeatures = row_subfeatures[row_subfeatures.type == "CDS"]

            if len(row_subfeatures) == 0:
                if row['strand'] == "-":
                    cds_three_p_end = row.start
                else:
                    cds_three_p_end = row.end
            else:
                if row['strand'] == "-":
                    cds_three_p_end = row_subfeatures.iloc[0].start
                else:
                    cds_three_p_end = row_subfeatures.iloc[-1].end

            # NOTE: now the longest function in the TES analysis
            for i in range(len(reads_in_region)):
                r = reads_in_region[i]

                # keep only reads in the same direction as this strand
                if (row['strand'] == "+" and r.is_forward) or (row['strand'] == "-" and r.is_reverse):
                    if row['strand'] == "-":
                        read_3p_end = r.reference_start
                        read_3p_start = r.reference_end

                        if read_3p_end <= cds_three_p_end and read_3p_start >= cds_three_p_end:
                            read_indexes_to_process.append(i)
                        else:
                            read_outside_3p_end.append(r.query_name)

                    else:
                        read_3p_end = r.reference_end
                        read_3p_start = r.reference_start

                        if read_3p_end >= cds_three_p_end and read_3p_start <= cds_three_p_end:
                            read_indexes_to_process.append(i)
                        else:
                            read_outside_3p_end.append(r.query_name)
                else:
                    read_on_different_strand.append(r.query_name)

            filtered_poly_a = 0
            if POLY_A_FILTER > 0:
                reads_with_poly_a = []
                for i in read_indexes_to_process:
                    r = reads_in_region[i]
                    if r.has_tag('pt:i'):
                        poly_a_length = r.get_tag('pt:i')
                        if poly_a_length >= POLY_A_FILTER:
                            reads_with_poly_a.append(i)

                # print("INFO: filtered {} reads with poly a filter <= {}".format(len(read_indexes_to_process) - len(reads_with_poly_a), POLY_A_FILTER))
                filtered_poly_a = len(read_indexes_to_process) - len(reads_with_poly_a)
                read_indexes_to_process = reads_with_poly_a

            no_poly_a = 0
            poly_a_lengths = []
            tts_sites = []
            utr_lengths = []

            for i in read_indexes_to_process:
                r = reads_in_region[i]

                # TODO: does r.reference start / end give the true transcipt end site, or just where it finishes mapping to reference?
                if row['strand'] == "-":
                    tts_sites.append(r.reference_start)
                    # TODO: account for soft clipping here
                    utr_lengths.append(cds_three_p_end - r.reference_start)
                else:
                    tts_sites.append(r.reference_end)
                    utr_lengths.append(r.reference_end - cds_three_p_end)

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
            d_utr_lengths[label][row['ID']] = utr_lengths

            d_not_beyond_3p[label][row['ID']] = len(read_outside_3p_end)


            row_summary = [label, row.ID, row.strand, len(reads_in_region), len(tts_sites), len(read_on_different_strand), len(read_outside_3p_end), filtered_poly_a]
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
        "contig",
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
        "pa_site_proportions",
        "mean_utr_length",
        "mean_poly_a_length"
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
            if len(d_tts[label][row.ID]) < READ_DEPTH_THRESHOLD:
                min_tts = 0
                max_tts = 0
                max_poly_a = 0
            else:
                if min_tts > min(d_tts[label][gene_id]) or min_tts == 0:
                    min_tts = min(d_tts[label][gene_id])

                if max_tts < max(d_tts[label][gene_id]):
                    max_tts = max(d_tts[label][gene_id])

                if max_poly_a < max(d_poly_a_lengths[label][gene_id]):
                    max_poly_a = max(d_poly_a_lengths[label][gene_id])
        
        x_ticks = range(min_tts, max_tts)

        # calculate hists
        # print("{} - Generating transcript end site histograms...".format(row['ID']))
        max_hist_count_poly_a = 0
        max_hist_count_tts = 0

        for label in bam_labels:
            if len(d_tts[label][row.ID]) < READ_DEPTH_THRESHOLD or len(d_poly_a_lengths[label][gene_id]) < READ_DEPTH_THRESHOLD:
                d_poly_a_length_hists[label][row['ID']] = []
                d_tts_hist[label][row['ID']] = []
            else:
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
        max_density = 0
        if len(matches) == 1:
            # print("{} - Generating transcript end site density information...".format(row['ID']))
            for label in bam_labels:
                if len(d_tts[label][gene_id]) < READ_DEPTH_THRESHOLD:
                    d_kdes[label][gene_id] = [0] * len(x_ticks)
                else:
                    kernel = gaussian_kde(d_tts[label][gene_id], bw_method=0.1)
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
        d_x_ticks[gene_id] = x_ticks

        for label in bam_labels:
            if len(d_tts[label][row.ID]) < READ_DEPTH_THRESHOLD:
                d_kdes[label][row['ID']] = [0] * len(x_ticks)
                num_tts = 0
                canonical_pa_site = 0
                apa_score = 0
                genomic_apa_sites = []
                pa_site_counts = []
                pa_site_proportions = []
                d_genomic_apa_sites[label][row.ID] = []
                mean_utr_length = 0
                mean_poly_a_length = 0
            else:
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
                    pa_site_proportions.append(rt_prop)

                num_tts = len(d_tts[label][gene_id])
                pa_site_counts = list(peak_dict['peak_heights'])
                d_genomic_apa_sites[label][row.ID] = genomic_apa_sites
                mean_utr_length = round(numpy.mean(d_utr_lengths[label][gene_id]), NUM_DPS)
                mean_poly_a_length = round(numpy.mean(d_poly_a_lengths[label][gene_id]), NUM_DPS)

            # TODO add max_read_depth
            row_summary = [row.seq_id, label, row.ID, row.type, row.strand, annotation_row_3p_end, num_tts, canonical_pa_site, round(apa_score, NUM_DPS), genomic_apa_sites, [int(x) for x in pa_site_counts], [round(float(x), NUM_DPS) for x in pa_site_proportions], mean_utr_length, mean_poly_a_length]
            print("\t".join(map(str, row_summary)))
            raw_results.append(row_summary)

    raw_summary_df = pandas.DataFrame(raw_results, columns=raw_summary_header)

    if OUTFILE:
        raw_summary_df.to_csv(OUTPUT_RAW, sep='\t', index=False)

    # go through d_kdes, find all local max's with count > read_depth threshold and call these poly_adenylation sites
    # The max PA is the canonical poly_adenylation site, and belongs in it's own column
    if COMPARE_APA_BETWEEN_TREATMENTS:
        summary_header = [
            "contig",
            "ID",
            "type",
            "strand",
            "3p_annotation_end",
            "canonical_pa_site",
            "num_tx_used",
            "mean_rt_prop_g1",
            "mean_rt_prop_g2",
            "test_stat_rt_prop",
            "p_val_rt_prop",
            "mean_utr_length_g1",
            "mean_utr_length_g2",
            "test_stat_utr_length",
            "p_val_utr_length",
            "mean_poly_a_length_g1",
            "mean_poly_a_length_g2",
            "test_stat_poly_a_length",
            "p_val_poly_a_length"
        ]
        print("\t".join(map(str, summary_header)))
        results = []
        group_1_prefix = COMPARE_APA_BETWEEN_TREATMENTS[0]
        group_2_prefix = COMPARE_APA_BETWEEN_TREATMENTS[1]

        group1_rows = raw_summary_df[raw_summary_df['label'].str.startswith(group_1_prefix)]
        group2_rows = raw_summary_df[raw_summary_df['label'].str.startswith(group_2_prefix)]

        group1_bam_labels = [l for l in bam_labels if l.startswith(group_1_prefix)]
        group2_bam_labels = [l for l in bam_labels if l.startswith(group_2_prefix)]

        # go through each row, determine quantify change in 
        for row_index, row in matches.iterrows():
            gene_id = row.ID
            if row['strand'] == "+":
                annotation_row_3p_end = row['end']
            else:   
                annotation_row_3p_end = row['start']
            these_rows_g1 = group1_rows[group1_rows['ID'] == row.ID]
            these_rows_g2 = group2_rows[group2_rows['ID'] == row.ID]

            # look through both treatments and keep only apa sites common to all of them
            canonical_pas = these_rows_g1.canonical_pa_site.to_list()
            canonical_pas = [x for x in canonical_pas if x > 0]

            num_tx_used = []
            for label in bam_labels:
                num_tx_used.append(len(d_tts[label][gene_id]))

            if len(canonical_pas) > 0:
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
                    group2_num_read_throughs.append(num_read_throughs)
                    group2_read_through_props.append(num_read_throughs / len(d_tts[label][gene_id]))

                group1_average_rt_prop = sum(group1_read_through_props) / len(group1_read_through_props)
                group2_average_rt_prop = sum(group2_read_through_props) / len(group2_read_through_props)

                num_tes = [len(d_tts[l][gene_id]) for l in bam_labels]

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
                rt_prop_test_stat = result.params.group_binary
                rt_prop_pval = result.pvalues.group_binary

                # test for change in UTR length
                group1_utr_lengths = []
                group2_utr_lengths = [] 

                for label in group1_bam_labels:
                    group1_utr_lengths += d_utr_lengths[label][gene_id]
                for label in group2_bam_labels:
                    group2_utr_lengths += d_utr_lengths[label][gene_id]

                utr_length_test_stat, utr_length_pval = kruskal(group1_utr_lengths, group2_utr_lengths)
                group1_mean_utr = numpy.mean(group1_utr_lengths)
                group2_mean_utr = numpy.mean(group2_utr_lengths)

                # test for change in poly_a length
                group1_poly_a_lengths = []
                group2_poly_a_lengths = [] 

                for label in group1_bam_labels:
                    group1_poly_a_lengths += d_poly_a_lengths[label][gene_id]
                for label in group2_bam_labels:
                    group2_poly_a_lengths += d_poly_a_lengths[label][gene_id]

                test_stat_poly_a_length, p_val_poly_a_length = kruskal(group1_poly_a_lengths, group2_poly_a_lengths)
                mean_poly_a_length_g1 = numpy.mean(group1_poly_a_lengths)
                mean_poly_a_length_g2 = numpy.mean(group2_poly_a_lengths)

            else:
                average_canonical_pa = 0
                group1_average_rt_prop = 0
                group2_average_rt_prop = 0
                rt_prop_test_stat = 0
                rt_prop_pval = 0

                group1_mean_utr = 0
                group2_mean_utr = 0
                utr_length_test_stat = 0
                utr_length_pval = 0

                mean_poly_a_length_g1 = 0
                mean_poly_a_length_g2 = 0
                test_stat_poly_a_length = 0
                p_val_poly_a_length = 0

            row_summary = [
                row.seq_id,
                row.ID, 
                row.type, 
                row.strand, 
                annotation_row_3p_end, 
                average_canonical_pa, 
                num_tx_used,
                round(group1_average_rt_prop, NUM_DPS), 
                round(group2_average_rt_prop, NUM_DPS), 
                rt_prop_test_stat, 
                rt_prop_pval,
                int(group1_mean_utr),
                int(group2_mean_utr),
                utr_length_test_stat,
                utr_length_pval,
                round(mean_poly_a_length_g1, NUM_DPS),
                round(mean_poly_a_length_g2, NUM_DPS),
                test_stat_poly_a_length,
                p_val_poly_a_length
            ]
            results.append(row_summary)
            print("\t".join(map(str, row_summary)))

        summary_df = pandas.DataFrame(results, columns=summary_header)

        if OUTFILE:
            summary_df.to_csv(OUTPUT_COMPARE, sep='\t', index=False)


    # --------- PLOT ---------- #
    if len(matches) == 1 and len(d_tts[label][row.ID]) >= READ_DEPTH_THRESHOLD:
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

