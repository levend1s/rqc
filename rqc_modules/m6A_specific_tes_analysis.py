import matplotlib.pyplot as plt
import scipy.stats
import os
import numpy

from rqc_modules.utils import process_input_files, process_annotation_file, filter_gff_for_target_features, find_canonical_mods, get_filtered_reads_ids


def m6A_specific_tes_analysis(args):
    ANNOTATION_FILE = args.annotation
    MOD_PROP_THRESHOLD = args.mod_prob_threshold
    MOD_RATIO = args.mod_ratio
    COVERAGE_PADDING = args.coverage_padding
    INPUT = args.input
    IDS = args.ids
    OUTPUT = args.output
    READ_DEPTH_THRESHOLD = args.read_depth
    OFFSET = args.offset
    SEPARATE_MOD_TRACKS = args.separate_mod_tracks
    OFFSET_PADDING = args.offset_padding

    print("m6A_specific_tes_analysis...")

    
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

 
    matches = gff_df[gff_df['ID'].isin(IDS)]

    if matches.empty:
        print("ERROR: no matches found for ids {}".format(IDS))

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']

    canonical_mods = find_canonical_mods(
        matches, 
        input_files, 
        bam_labels, 
        MOD_RATIO,
        READ_DEPTH_THRESHOLD, 
        COVERAGE_PADDING
    )
    print(canonical_mods)

    # approximated_tes = approximate_tes(target_features, input_files, bam_labels)
    approximated_tes = None

    filtered_read_ids = get_filtered_reads_ids(
        matches, 
        input_files,
        bam_labels, 
        MOD_PROP_THRESHOLD,
        READ_DEPTH_THRESHOLD,
        COVERAGE_PADDING,
        canonical_mods
        # approximated_tes = approximated_tes
    )

    # HACK manually remove total coverage for all mods
    for _, row in matches.iterrows():
        for label in bam_labels:
            filtered_read_ids[label][row.ID]['tx_end_sites'].pop('None', None)

    # TODO: make more readable, refactor
    for _, row in matches.iterrows():
        # create hist / kde for all tes, and separate ones for reads containing each mod
        min_tes = 0
        max_tes = 0
        max_hist_count_tes = 0
        max_density = 0
        d_tes_hist = {}
        d_kdes = {}

        for label in bam_labels:
            d_tes_hist[label] = {}
            for key, tx_end_sites in filtered_read_ids[label][row.ID]['tx_end_sites'].items():

                if min_tes > min(tx_end_sites) or min_tes == 0:
                    min_tes = min(tx_end_sites)

                if max_tes < max(tx_end_sites):
                    max_tes = max(tx_end_sites)

                # calculate hists
                unique_tes = set(tx_end_sites)
                tes_hist = [(i, tx_end_sites.count(i)) for i in unique_tes]

                # split the tuple cause here we're interested in the biggest count in the hist
                e1 = [e[1] for e in tes_hist]
                if max_hist_count_tes < max(e1):
                    max_hist_count_tes = max(e1)


                d_tes_hist[label][key] = tes_hist

        min_x = min_tes
        max_x = max_tes
        if min(canonical_mods[row.ID]) < min_x:
            min_x = min(canonical_mods[row.ID])
        if max(canonical_mods[row.ID]) > max_x:
            max_x = min(canonical_mods[row.ID])

        x_width = max_x - min_x
        x_ticks = range(min_x - int(x_width * 0.1), max_x + int(x_width * 0.1))

        for label in bam_labels:
            # generate dennsity plots
            print("{} - Generating transcript end site density information...".format(row.ID))
            d_kdes[label] = {}
            
            for key, tx_end_sites in filtered_read_ids[label][row.ID]['tx_end_sites'].items():
                
                kernel = scipy.stats.gaussian_kde(tx_end_sites, bw_method=0.3)
                smoothed_tes_hist = kernel(x_ticks)
                d_kdes[label][key] = smoothed_tes_hist

                if max_density < max(smoothed_tes_hist):
                    max_density = max(smoothed_tes_hist)

        if row.strand == "-":
            row_annotation_3p_end = row.start
            # x_ticks = x_ticks[::-1]
        else:
            row_annotation_3p_end = row.end

        offset = row_annotation_3p_end

        if OFFSET:
            offset = OFFSET

        if row.strand == "-":
            x_tick_offset = numpy.array([offset - x for x in x_ticks])
        else:
            x_tick_offset = numpy.array([x - offset for x in x_ticks])

        # assign mod colours to mod locations
        d_mod_colours = {
            'all': '#cbcbcb',
            'None': 'black'
        }
        mod_colours = ['green', 'orange', 'olive', 'cyan', 'red', 'purple', 'magenta', 'brown', 'blue', 'yellow']
        mod_colours = ['#2e0064', '#2e0064', '#2e0064', '#2e0064', '#2e0064', '#2e0064', '#2e0064']
        mod_colours = ['#cbcbcb', '#cbcbcb', '#cbcbcb', '#cbcbcb', '#cbcbcb', '#cbcbcb', '#cbcbcb']

        for mod_location in canonical_mods[row.ID]:
            d_mod_colours[mod_location] = mod_colours.pop(0)

        # HACK manually remove total coverage for all mods
        d_tes_hist[label].pop('all', None)
        d_kdes[label].pop('all', None)

        plot_data = d_tes_hist

        if len(bam_labels) == 1 and SEPARATE_MOD_TRACKS:
                NUM_VERT_PLOTS = 1 + len(filtered_read_ids[label][row.ID]['tx_end_sites'].keys())
                # HISTS
                fig, axes = plt.subplots(len(plot_data[label].items()), sharex=True, sharey=False)
                axes_index = 0
                this_axes = axes

                for label in bam_labels:
                    if len(bam_labels) == 1:
                        this_axes = axes
                    else:
                        this_axes = axes[axes_index]


                hists_to_plot = list(plot_data[label].items())
                if row.strand == "-":
                    hists_to_plot = hists_to_plot[::-1]

                for key, hist in hists_to_plot:
                    this_axes = axes[axes_index]
                    if key != "all":
                        label_offset = "offset: {}".format(int(key) - offset)
                    else:
                        label_offset = key
                    l = "{}\nnum reads={}".format(label_offset, len(filtered_read_ids[label][row.ID]['tx_end_sites'][key]))
                    sorted_tes_counts_by_pos = sorted(hist, key=lambda a: a[0])
                    x_coords, y_coords = zip(*sorted_tes_counts_by_pos)
                    # add y=0
                    y_coords_w_0 = []
                    for x in x_ticks:
                        try:
                            idx = x_coords.index(x)
                            y_coords_w_0.append(y_coords[idx])
                        except ValueError:
                            y_coords_w_0.append(0)


                    # mask = (x_tick_offset >= -OFFSET_PADDING) & (x_tick_offset <= OFFSET_PADDING)
                    # x_filtered = x_tick_offset[mask]
                    # y_filtered = numpy.array(y_coords_w_0)[mask]

                    this_axes.step(x_tick_offset, y_coords_w_0, label=l, where="mid", color=d_mod_colours[key])
                    this_axes.fill_between(x_tick_offset, y_coords_w_0, alpha=1, step="mid", color=d_mod_colours[key])

                    if key != "all":
                        all_pos = int(key) - offset

                        if row.strand == "-":
                            all_pos = offset - int(key)
                        this_axes.axvline(x=all_pos, color='#2e0064', ls="--", linewidth=1.0)
                    else:
                        # # PLOT MOD LOCATION AS VERT LINES
                        for mod_location in canonical_mods[row.ID]:
                            mod_pos = mod_location - offset
                            if row.strand == "-":
                                mod_pos = offset - mod_location
                            this_axes.axvline(x=mod_pos, color='#2e0064', ls="--", linewidth=1.0)

                    if OFFSET:
                        # also plot the real gene end
                        three_p_pos = row_annotation_3p_end - offset
                        if row.strand == "-":
                            three_p_pos = offset - row_annotation_3p_end
                        this_axes.axvline(x=three_p_pos, color='green', ls="--", linewidth=1.0)

                    # this_axes.axvline(x=0, color='darkgray', ls="--", linewidth=1.0)

                    # this_axes.legend()
                    this_axes.set_ylim(ymin=0)

                    axes_index += 1

                this_axes.set_xlim(xmin=-200, xmax=1000)

                # PLOT APPROXIMATED TES AS VERT LINE
                # if approximated_tes:
                #     this_axes.axvline(x=approximated_tes[row.ID], color='darkgray', ls=":", linewidth=1.0)

                # PLOT GENE END AS VERT LINE
                # this_axes.axvline(x=0, color='darkgray', ls="--", linewidth=1.0)

                # if OFFSET:
                #     # also plot the real gene end
                #     this_axes.axvline(x=row_annotation_3p_end - offset, color='darkgray', ls="--", linewidth=1.0)

                # this_axes.axvline(x=offset, color='darkgray', ls="--", linewidth=1.0)

                # # PLOT MOD LOCATION AS VERT LINES
                # for mod_location in canonical_mods[row.ID]:
                #     this_axes.axvline(x=mod_location - offset, color=d_mod_colours[mod_location], ls="--", linewidth=1.0)

                # add axis labels
                sample_name = label.split('_')[0]
                this_axes.set(xlabel='offset from canonical PAS (nt)', ylabel='{} - count'.format(sample_name))

        else:
            NUM_VERT_PLOTS = 1 + len(filtered_read_ids[label][row.ID]['tx_end_sites'].keys())
            # HISTS
            fig, axes = plt.subplots(len(bam_labels), sharex=True, sharey=True)
            axes_index = 0
            for label in bam_labels:
                if len(bam_labels) == 1:
                    this_axes = axes
                else:
                    this_axes = axes[axes_index]

                for key, hist in d_tes_hist[label].items():
                    l = "{}, n={}".format(key, len(filtered_read_ids[label][row.ID]['tx_end_sites'][key]))
                    sorted_tes_counts_by_pos = sorted(hist, key=lambda a: a[0])
                    x_coords, y_coords = zip(*sorted_tes_counts_by_pos)
                    # add y=0
                    y_coords_w_0 = []
                    for x in x_ticks:
                        try:
                            idx = x_coords.index(x)
                            y_coords_w_0.append(y_coords[idx])
                        except ValueError:
                            y_coords_w_0.append(0)

                    this_axes.plot(x_tick_offset, y_coords_w_0, label=l, color=d_mod_colours[key])
                    this_axes.fill_between(x_tick_offset, y_coords_w_0, alpha=0.2, color=d_mod_colours[key])

                this_axes.set_ylim(ymin=0, ymax=max_hist_count_tes*1.1)
                this_axes.set_xlim(xmin=-OFFSET_PADDING, xmax=OFFSET_PADDING)
                # this_axes.legend()

                # PLOT APPROXIMATED TES AS VERT LINE
                if approximated_tes:
                    this_axes.axvline(x=approximated_tes[row.ID], color='darkgray', ls=":", linewidth=1.0)

                # PLOT GENE END AS VERT LINE
                this_axes.axvline(x=0, color='darkgray', ls="--", linewidth=1.0)

                if OFFSET:
                    # also plot the real gene end
                    this_axes.axvline(x=row_annotation_3p_end - offset, color='darkgray', ls="--", linewidth=1.0)

                this_axes.axvline(x=offset, color='darkgray', ls="--", linewidth=1.0)

                # # PLOT MOD LOCATION AS VERT LINES
                for mod_location in canonical_mods[row.ID]:
                    this_axes.axvline(x=mod_location - offset, color=d_mod_colours[mod_location], ls="--", linewidth=1.0)

                # add axis labels
                sample_name = label.split('_')[0]
                this_axes.set(xlabel='offset from canonical PAS (nt)', ylabel='{} - count'.format(sample_name))
                    
                axes_index += 1

        fig.tight_layout()
        # fig.suptitle("transcript end sites filtered for reads containing a cm6A")
        # fig.subplots_adjust(hspace=0, wspace=0.1)

        # KDEs
        # fig, axes = plt.subplots(len(bam_labels), sharex=True, sharey=True)
        # axes_index = 0
        # for label in bam_labels:
        #     for key, kde in d_kdes[label].items():
        #         l = "{}, n={}".format(key, len(d_tes_hist[label][key]))
        #         axes[axes_index].plot(x_ticks, kde, label=l, color=d_mod_colours[key])
        #         axes[axes_index].fill_between(x_ticks, kde, alpha=0.2, color=d_mod_colours[key])

        #     axes[axes_index].set_ylim(ymin=0, ymax=max_density*1.1)
        #     axes[axes_index].set_xlim(xmin=x_ticks[0], xmax=x_ticks[-1])
        #     axes[axes_index].legend()

        #     # PLOT GENE END AS VERT LINE
        #     axes[axes_index].axvline(x=row_annotation_3p_end, color='darkgray', ls="--", linewidth=1.0)

        #     # # PLOT MOD LOCATION AS VERT LINES
        #     for mod_location in canonical_mods[row.ID]:
        #         axes[axes_index].axvline(x= mod_location, color=d_mod_colours[mod_location], ls="--", linewidth=1.0)

        #     # add axis labels
        #     sample_name = label.split('_')[0]
        #     axes[axes_index].set(xlabel='transcription end site (nt)', ylabel='density (au) ({})'.format(sample_name))
                
        #     axes_index += 1

        # fig.tight_layout()
        # # fig.suptitle("transcript end sites filtered for reads containing a cm6A")
        # fig.subplots_adjust(hspace=0, wspace=0.1)

    if OUTPUT:
        OUTPUT_FORMAT = OUTPUT.split(".")[-1] if OUTPUT else "png"

        if OUTPUT_FORMAT not in ["png", "eps", "pdf"]:
            raise ValueError("Output format must be one of: png, eps, pdf")
        
        plt.savefig("m6A_specific_tes_analysis_{}".format(OUTPUT), transparent=True, dpi=300, format=OUTPUT_FORMAT)

    plt.show()
