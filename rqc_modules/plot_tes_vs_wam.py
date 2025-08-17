def plot_tes_vs_wam(args):
    tes_tsv_file_path = args.inputs[0]
    print("LOADING: {}".format(tes_tsv_file_path))
    tes_file_df = pandas.read_csv(tes_tsv_file_path, sep='\t')

    neighbour_file_df = {}

    if NEIGHBOUR_FILE:
        with open(NEIGHBOUR_FILE) as json_data:
            neighbour_file_df = json.load(json_data)

    # print(neighbour_file_df.keys())

    # FILTER_BY_NEIGHBOUR_TYPE = "lonely"
    tes_file_df_len_raw = len(tes_file_df)
    print("args.inputs: {}".format(tes_file_df_len_raw))

    if FILTER_BY_NEIGHBOUR_TYPE != "all":
        gene_list_filter = neighbour_file_df[FILTER_BY_NEIGHBOUR_TYPE]
        tes_file_df["parent_id"] = tes_file_df.gene_id.apply(lambda s: s.split('.')[0])

        # flatten list if required
        if len(gene_list_filter) > 0 and isinstance(gene_list_filter[0], list):
            gene_list_filter = [x for xs in gene_list_filter for x in xs]

        gene_list_filter = set(gene_list_filter)

        # remove all entries from tes_file if the gene isn't in the gene list
        tes_file_df = tes_file_df[
                (tes_file_df['parent_id'].isin(gene_list_filter))
        ]

        print("REMOVED {} DUE TO FILTER (GENE NEIGHOUR={})".format(tes_file_df_len_raw - len(tes_file_df), FILTER_BY_NEIGHBOUR_TYPE))

    tes_file_df_len_raw = len(tes_file_df)
    tes_file_df = tes_file_df[
        (tes_file_df.wart_after < 1.0) & 
        (tes_file_df.wart_before < 1.0)
    ]
    print("REMOVED {} DUE TO 'IMPOSSIBLE' WEIRD TES CHANGE".format(tes_file_df_len_raw - len(tes_file_df)))



    tes_file_df['minus_log10_p_inter_treatment'] = (numpy.log10(tes_file_df['p_inter_treatment']) * -1)
    tes_file_df['log2_average_expression'] = (numpy.log2(tes_file_df['average_expression']))
    tes_file_df['-log2_wam_change'] = (numpy.log2(tes_file_df['wam_change']) * -1)
    tes_file_df['log2_wart_change'] = numpy.log2(tes_file_df['wart_change'])

    tes_file_df['wam_diff'] = (tes_file_df['wam_before'] - tes_file_df['wam_after']) * 100
    tes_file_df['tes_diff'] = (tes_file_df['wart_after'] - tes_file_df['wart_before']) * 100


    # drop all genes where p_same_treatment < 0.05 (ie the same conditions don't have same TES)
    # drop all genes where wam_change == 0
    p_same_treatment_cutoff = 0.05
    MIN_GAP_BETWEEN_M6A = 1
    num_smeared = 0

    num_canonical_mods = []
    for _, row in tes_file_df.iterrows():
        canonical_mods = sorted([int(s) for s in ast.literal_eval(row['cannonical_mods'])])
        this_num_canonical_mods = len(canonical_mods)

        if this_num_canonical_mods <= 1:
            num_canonical_mods.append(this_num_canonical_mods)
        else:
            mod_distances = []
            prev = 0
            for x in canonical_mods:
                if prev == 0:
                    prev = x
                else:
                    mod_distances.append(x - prev)
                    prev = x

            mod_distances = [x for x in mod_distances if x > MIN_GAP_BETWEEN_M6A]

            num_canonical_mods.append(len(mod_distances) + 1)

            if this_num_canonical_mods > 1 and len(mod_distances) != this_num_canonical_mods - 1:
                #print(canonical_mods)
                #print(mod_distances)
                num_smeared += 1
                #print("NOTE: {} had {} m6As that were too close (<={}nt), ...".format(row['gene_id'], this_num_canonical_mods - len(mod_distances), MIN_GAP_BETWEEN_M6A))

    print("NOTE: {} GENES HAD m6A SMEARING".format(num_smeared))

    tes_file_df["num_cannonical_mods"] = num_canonical_mods

    filtered_genes_tes_wam = tes_file_df[
        # (tes_file_df.p_same_treatment >= p_same_treatment_cutoff) &
        (tes_file_df.num_cannonical_mods > 0) & 
        (tes_file_df.average_expression >= READ_DEPTH_THRESHOLD)
    ]
    print("REMOVING {} DUE TO FILTER (num canonical mods > 0, avg expression > {})".format(len(tes_file_df) - len(filtered_genes_tes_wam), READ_DEPTH_THRESHOLD))

    import mplcursors


    axes = None
    if SPLIT_BY_CANONICAL_MODS:
        for i in range(1, max(filtered_genes_tes_wam["num_cannonical_mods"].to_list()) + 1):
            filtered_genes_tes_wam_mods = filtered_genes_tes_wam[
                (filtered_genes_tes_wam.num_cannonical_mods == i)
            ]

            if len(filtered_genes_tes_wam_mods) == 0:
                continue

            x_col = 'wam_diff'
            y_col = 'tes_diff'

            print("REMOVING {} DUE TO FILTER (MODS={})".format(len(filtered_genes_tes_wam) - len(filtered_genes_tes_wam_mods), i))

            axes = filtered_genes_tes_wam_mods.plot.scatter(
                x='wam_diff',
                y='tes_diff',
                c='log2_average_expression'
            )

            m, c, r_value, p_value, std_err = scipy.stats.linregress(filtered_genes_tes_wam_mods[x_col], filtered_genes_tes_wam_mods[y_col])
            # axes.plot(filtered_genes_tes_wam_mods[x_col], m * filtered_genes_tes_wam_mods[x_col] + c)
            # axes.text(1, 1, "R^2: {}".format(round(r_value ** 2, 2)), transform=axes.transAxes, horizontalalignment='right', verticalalignment='top')

            axes.set_title("{} genes with {} cannonical m6A (n={})".format(FILTER_BY_NEIGHBOUR_TYPE, i, len(filtered_genes_tes_wam_mods)))
    else:
        if NUM_CANONICAL_MODS_FILTER > 0:
            filtered_genes_tes_wam = filtered_genes_tes_wam[
                (filtered_genes_tes_wam.num_cannonical_mods == NUM_CANONICAL_MODS_FILTER)
            ]
        
        x_col = 'wam_diff'
        y_col = 'tes_diff'

        axes = filtered_genes_tes_wam.plot.scatter(
            x='wam_diff',
            y='tes_diff',
            # c='log2_average_expression',
            s=10
        )
        m, c, r_value, p_value, std_err = scipy.stats.linregress(filtered_genes_tes_wam[x_col], filtered_genes_tes_wam[y_col])
        # axes.plot(filtered_genes_tes_wam[x_col], m * filtered_genes_tes_wam[x_col] + c)
        axes.set_ylim(ymin=0, ymax=100)
        axes.set_xlim(xmin=0, xmax=100)
        axes.set_xlabel("")
        axes.set_ylabel("")


        # axes.text(1, 1, "R: {}".format(round(r_value, 2)), transform=axes.transAxes, horizontalalignment='right', verticalalignment='top')
        # axes.set_xticks([round(x*0.1, 1) for x in range(0, 11)])
        print("R: {}".format(round(r_value, 2)))
        print("{} genes with cannonical m6A (n={})".format(FILTER_BY_NEIGHBOUR_TYPE, len(filtered_genes_tes_wam)))

        # if NUM_CANONICAL_MODS_FILTER > 0:
        #     axes.set_title("{} genes with {} cannonical m6A (n={})".format(FILTER_BY_NEIGHBOUR_TYPE, NUM_CANONICAL_MODS_FILTER, len(filtered_genes_tes_wam)))
        # else:
        #     axes.set_title("{} genes with cannonical m6A (n={})".format(FILTER_BY_NEIGHBOUR_TYPE, len(filtered_genes_tes_wam)))

        def show_label(sel):
            index = sel.index
            sel.annotation.set_text(filtered_genes_tes_wam['gene_id'].to_list()[index])
            print(filtered_genes_tes_wam['gene_id'].to_list()[index])
            
        mplcursors.cursor(axes, hover=True).connect("add", show_label)

    plt.savefig("plot_tes_vs_wam.png", dpi=300)
    plt.show()

