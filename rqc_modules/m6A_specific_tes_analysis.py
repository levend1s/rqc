def m6A_specific_tes_analysis(args):
    print("m6A_specific_tes_analysis...")

    input_files = process_input_files()

    gff = process_annotation_file()
    target_features = filter_gff_for_target_features(gff)

    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']
    bam_labels_control = [l for l in bam_labels if input_files[l]['group'] == 'control']
    bam_labels_treatment = [l for l in bam_labels if input_files[l]['group'] == 'knock-sideways']

    canonical_mods = find_canonical_mods(
        target_features, 
        input_files, 
        bam_labels, 
        CANNONICAL_MOD_PROP_THRESHOLD, 
        READ_DEPTH_THRESHOLD, 
        COVERAGE_PADDING
    )
    print(canonical_mods)

    # approximated_tes = approximate_tes(target_features, input_files, bam_labels)
    approximated_tes = None

    filtered_read_ids = get_filtered_reads_ids(
        target_features, 
        input_files,
        bam_labels, 
        0.95,
        READ_DEPTH_THRESHOLD,
        COVERAGE_PADDING,
        canonical_mods
        # approximated_tes = approximated_tes
    )

    # HACK manually remove total coverage for all mods
    for _, row in target_features.iterrows():
        for label in bam_labels:
            filtered_read_ids[label][row.ID]['tx_end_sites'].pop('None', None)

    # TODO: make more readable, refactor
    for _, row in target_features.iterrows():
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
        else:
            row_annotation_3p_end = row.end

        # assign mod colours to mod locations
        d_mod_colours = {
            'all': 'lightsteelblue',
            'None': 'black'
        }
        mod_colours = ['green', 'orange', 'olive', 'cyan']
        for mod_location in canonical_mods[row.ID]:
            d_mod_colours[mod_location] = mod_colours.pop(0)


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

                this_axes.plot(x_ticks, y_coords_w_0, label=l, color=d_mod_colours[key])
                this_axes.fill_between(x_ticks, y_coords_w_0, alpha=0.2, color=d_mod_colours[key])

            this_axes.set_ylim(ymin=0, ymax=max_hist_count_tes*1.1)
            this_axes.set_xlim(xmin=x_ticks[0], xmax=x_ticks[-1])
            this_axes.legend()

            # PLOT APPROXIMATED TES AS VERT LINE
            if approximated_tes:
                this_axes.axvline(x=approximated_tes[row.ID], color='darkgray', ls=":", linewidth=1.0)

            # PLOT GENE END AS VERT LINE
            this_axes.axvline(x=row_annotation_3p_end, color='darkgray', ls="--", linewidth=1.0)

            # # PLOT MOD LOCATION AS VERT LINES
            for mod_location in canonical_mods[row.ID]:
                this_axes.axvline(x=mod_location, color=d_mod_colours[mod_location], ls="--", linewidth=1.0)

            # add axis labels
            sample_name = label.split('_')[0]
            this_axes.set(xlabel='transcription end site (nt)', ylabel='{} - count'.format(sample_name))
                
            axes_index += 1

        fig.tight_layout()
        # fig.suptitle("transcript end sites filtered for reads containing a cm6A")
        fig.subplots_adjust(hspace=0, wspace=0.1)

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
        plt.show()


