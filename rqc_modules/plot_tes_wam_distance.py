def plot_tes_wam_distance(args):
    tes_tsv_file_path = args.inputs[0]
    print("LOADING: {}".format(tes_tsv_file_path))
    tes_file_df = pandas.read_csv(tes_tsv_file_path, sep='\t')

    # drop all genes where p_same_treatment < 0.05 (ie the same conditions don't have same TES)
    # drop all genes where wam_change == 0
    p_same_treatment_cutoff = 0.05

    tes_file_df["num_cannonical_mods"] = tes_file_df.cannonical_mods.apply(lambda s: len(list(ast.literal_eval(s))))
    tes_file_df["cannonical_mods"] = tes_file_df.cannonical_mods.apply(lambda s: list(ast.literal_eval(s)))

    # flatten 2d list of cannonical mods
    all_cannonical_mods = [x for xs in tes_file_df["cannonical_mods"].to_list() for x in xs]

    tes_split_sites = tes_file_df["tes"].to_list()

    # load annotation file and find indexes for all parent children
    ANNOTATION_FILE = gffpandas.read_gff3(ANNOTATION_FILE_PATH)
    GFF_DF = ANNOTATION_FILE.attributes_to_columns()
    GFF_DF['ID'] = GFF_DF['ID'].astype('category')

    cannonical_mod_offsets = []
    annotation_start_offsets = []
    annotation_end_offsets = []
    tes_end_offsets = []


    # TODO also plot DRACH sites
    # TODO also plot DRACH sites
    # TODO also plot DRACH sites

    for row_index, row in tes_file_df.iterrows():
        this_row_gff = GFF_DF[GFF_DF['ID'] == row['gene_id']]
        num_matches = len(this_row_gff)
        if num_matches != 1:
            print("ERROR: found {} matches for {}".format(num_matches, row['ID']))
            continue

        row_strand = this_row_gff.iloc[0]['strand']
        row_start = this_row_gff.iloc[0]['start']
        row_end = this_row_gff.iloc[0]['end']

        if REFERENCE_POINT == "3_PRIME":
            if row_strand == "-":
                reference_point = row_start
            else:
                reference_point = row_end

            reference_label = "annotated 3' end"
        if REFERENCE_POINT == "TES":
            reference_point = row["tes"]
            reference_label = "approximated TES"

        # assume it's -ve
        row_mod_offsets = reference_point - numpy.array(row["cannonical_mods"])
        row_start_offset = reference_point - row_end
        row_end_offset = reference_point - row_start
        row_tes_offset = reference_point - row['tes']

        if row_strand == "+":
            row_mod_offsets *= -1
            row_start_offset = row_start - reference_point
            row_end_offset = row_end - reference_point
            row_tes_offset = row['tes'] - reference_point

        for x in row_mod_offsets:
            cannonical_mod_offsets.append(x)

        annotation_start_offsets.append(row_start_offset)
        annotation_end_offsets.append(row_end_offset)
        tes_end_offsets.append(row_tes_offset)

    min_x = int(min([
        min(cannonical_mod_offsets),
        # min(annotation_start_offsets),
        min(tes_end_offsets),
        min(annotation_end_offsets)
    ]) * 1.1)
    max_x = int(max([
        max(cannonical_mod_offsets),
        # max(annotation_start_offsets),
        max(tes_end_offsets),
        max(annotation_end_offsets)
    ]) * 1.1)

    x_ticks = numpy.linspace(
        min_x, 
        max_x, 
        max_x - min_x
    )
    if len(tes_file_df) > 1:
        kernel = scipy.stats.gaussian_kde(cannonical_mod_offsets)
        cannonical_mod_offset_kde = kernel(x_ticks)
        kernel = scipy.stats.gaussian_kde(annotation_start_offsets)
        annotation_start_offset_kde = kernel(x_ticks)
        if REFERENCE_POINT == "TES":
            kernel = scipy.stats.gaussian_kde(annotation_end_offsets)
            annotation_end_offset_kde = kernel(x_ticks)
        if REFERENCE_POINT == "3_PRIME":
            kernel = scipy.stats.gaussian_kde(tes_end_offsets)
            tes_offset_kde = kernel(x_ticks)

    cannonical_mod_offsets_hist = [cannonical_mod_offsets.count(i) for i in range(min_x, max_x)]
    annotation_start_offsets_hist = [annotation_start_offsets.count(i) for i in range(min_x, max_x)]
    annotation_end_offsets_hist = [annotation_end_offsets.count(i) for i in range(min_x, max_x)]
    tes_offsets_hist = [tes_end_offsets.count(i) for i in range(min_x, max_x)]

    d_colors = {
        'mods': 'green',
        'start': 'red',
        'end': 'blue',
    }
    if len(tes_file_df) > 1:
        fig, axes = plt.subplots()
        axes.plot(x_ticks, cannonical_mod_offset_kde, label='cannonical m6A', color=d_colors['mods'])
        axes.fill_between(x_ticks, cannonical_mod_offset_kde, alpha=0.2, color=d_colors['mods'])

        if REFERENCE_POINT == "TES":
            axes.plot(x_ticks, annotation_end_offset_kde, label='annotation end 3\'', color=d_colors['end'])
            axes.fill_between(x_ticks, annotation_end_offset_kde, alpha=0.2, color=d_colors['end'])
            axes.set_xlabel('distance from TES (nt)')

        if REFERENCE_POINT == "3_PRIME":
            axes.plot(x_ticks, tes_offset_kde, label='TES', color=d_colors['end'])
            axes.fill_between(x_ticks, tes_offset_kde, alpha=0.2, color=d_colors['end'])
            axes.set_xlabel('distance from 3\' (nt)')

        axes.axvline(x=0, color='grey', label=reference_label, ls="--", linewidth=1.0)
        axes.set_ylabel('density (au)')
        axes.legend()
        plt.legend(loc="upper right")


    fig, axes = plt.subplots()
    axes.plot(x_ticks, cannonical_mod_offsets_hist, label='cannonical m6A', color=d_colors['mods'])
    axes.fill_between(x_ticks, cannonical_mod_offsets_hist, alpha=0.2, color=d_colors['mods'])

    if REFERENCE_POINT == "TES":
        axes.plot(x_ticks, annotation_end_offsets_hist, label='annotation end 3\'', color=d_colors['end'])
        axes.fill_between(x_ticks, annotation_end_offsets_hist, alpha=0.2, color=d_colors['end'])
        axes.set_xlabel('distance from TES (nt)')

    if REFERENCE_POINT == "3_PRIME":
        axes.plot(x_ticks, tes_offsets_hist, label='TES', color=d_colors['end'])
        axes.fill_between(x_ticks, tes_offsets_hist, alpha=0.2, color=d_colors['end'])
        axes.set_xlabel('distance from 3\' (nt)')

    axes.axvline(x=0, color='grey', label=reference_label, ls="--", linewidth=1.0)
    axes.set_ylabel('count')
    axes.legend()
    plt.legend(loc="upper right")

    plt.show()

