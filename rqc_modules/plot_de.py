import json
import ast
import numpy
import pandas
import matplotlib.pyplot as plt
import mplcursors

def plot_de(args):

    de = pandas.read_csv(args.inputs[0], sep='\t')

    pval_cutoff = 0.05
    log10_pval_cutoff = 10 ** pval_cutoff
    fc_cutoff = 1

    # de_filtered = de[de["adj.P.Val"] < 0.05]
    de_filtered = de
    # print(de_filtered)

    de_analysis_type = None

    if 'adj.P.Val' in de_filtered.columns:
        de_analysis_type = "limma"
        pval_col_name = 'adj.P.Val'
    else:
        de_analysis_type = "edgeR"
        pval_col_name = 'FDR'
        de_filtered['gene_id'] = de_filtered.index
        print(de)

    de_filtered['-log10_adj_pval'] = (numpy.log10(de_filtered[pval_col_name]) * -1)
    de_filtered['gene_id'] = de_filtered['gene_id'].astype('category')
    de_filtered['parent_id'] = de_filtered['gene_id'].astype('category')
    # de_filtered = de_filtered.set_index('gene_id')
    
    de_filtered_raw_size = len(de_filtered)

    tes_file_df = {}

    if TES_ANALYSIS_FILE:
        print("LOADING: {}".format(TES_ANALYSIS_FILE))
        tes_file_df = pandas.read_csv(TES_ANALYSIS_FILE, sep='\t')
        tes_file_df['gene_id'] = tes_file_df['gene_id'].astype('category')
        tes_file_df["parent_id"] = tes_file_df.gene_id.apply(lambda s: s.split('.')[0])
        tes_file_df_raw_size = len(tes_file_df)

    tes_file_df = tes_file_df.sort_values(by="gene_id").drop_duplicates(subset=['parent_id'])
    print("REMOVED ISOFORMS FROM TES FILE: {}".format(tes_file_df_raw_size - len(tes_file_df)))

    neighbour_file_df = {}

    if NEIGHBOUR_FILE:
        with open(NEIGHBOUR_FILE) as json_data:
            neighbour_file_df = json.load(json_data)

    # remove mitochondrial and api genes
    de_filtered = de_filtered[
        (de_filtered['gene_id'].str.contains("MIT|API") == False)
    ]

    # if READ_DEPTH_THRESHOLD:
    #     log2rdf = numpy.log2(READ_DEPTH_THRESHOLD)
    #     de_filtered_prior_size = len(de_filtered)

    #     if de_analysis_type == "edgeR":
    #         de_filtered = de_filtered[de_filtered.logCPM >= log2rdf]
    #     else:
    #         de_filtered = de_filtered[de_filtered.AveExpr >= log2rdf]

    #     print("REMOVED {} DUE TO FILTER (READ DEPTH THRESHOLD >= {})".format(de_filtered_prior_size - len(de_filtered), READ_DEPTH_THRESHOLD))

    MIN_GAP_BETWEEN_M6A = 1
    num_smeared = 0

    num_canonical_mods = []
    if TES_ANALYSIS_FILE:
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
                    num_smeared += 1
                    #print("NOTE: {} had {} m6As that were too close (<={}nt), ...".format(row['gene_id'], this_num_canonical_mods - len(mod_distances), MIN_GAP_BETWEEN_M6A))

        print("NOTE: {} GENES HAD m6A SMEARING".format(num_smeared))

        tes_file_df["num_cannonical_mods"] = num_canonical_mods
        tes_file_df['tes_change'] = tes_file_df['wart_after'] = tes_file_df['wart_before']
        tes_file_df['wam_change'] = tes_file_df['wam_after'] = tes_file_df['wam_before']

        de_filtered = de_filtered.merge(
            tes_file_df[['parent_id', 'tes_change', 'wam_change', 'num_cannonical_mods']], 
            left_on='gene_id',
            right_on='parent_id', 
            how="outer"
        )

    de_filtered['is_downstream'] = 'lightgray'
    de_filtered = de_filtered.set_index('gene_id')

    # connecting_lines = []
    # de_filtered = de_filtered.set_index('gene_id')
    # de_dict = de_filtered.to_dict('index')
    # x coords are logFC, y coords are -log10_adj_pval
    # for each pair, get the xy coords of each element and plot a connecting line
    # if FILTER_BY_NEIGHBOUR_TYPE != "all" and SHOW_NEIGHBOURS:
    #     for a, b in neighbour_file_df[FILTER_BY_NEIGHBOUR_TYPE]:
    #         if a in de_dict.keys() and b in de_dict.keys():
    #             x = [de_dict[a]['logFC'], de_dict[b]['logFC']]
    #             y = [de_dict[a]['-log10_adj_pval'], de_dict[b]['-log10_adj_pval']]

    #             xy1 = [de_dict[a]['logFC'], de_dict[a]['-log10_adj_pval']]
    #             xy2 = [de_dict[b]['logFC'], de_dict[b]['-log10_adj_pval']]

    #             if y[0] >= log10_pval_cutoff or y[1] >= log10_pval_cutoff:
    #                 axes.plot(x, y, color='red', ls='-', linewidth=1.0)

    de_filtered["colorby"] = "lightgray"
    COLOR_BY = "updown"

    # neighbor_type_counts = {}
    # for t in neighbour_file_df.keys():
    #     neighbor_type_counts[t] = {}
    #     # flatten list if required
    #     if len(neighbour_file_df[t]) > 0 and isinstance(neighbour_file_df[t][0], list):
    #         neighbours_of_type = [x for xs in neighbour_file_df[t] for x in xs]
    #     else:
    #         neighbours_of_type = neighbour_file_df[t]

    #     for g in neighbours_of_type:
    #         neighbor_type_counts[t][g] = neighbours_of_type.count(g)
    genes_to_color_file_path = "edgeR_28hpi_overlaps_high_conf"

    if COLOR_BY == "methylation_discrete":
        de_filtered.loc[de_filtered['num_cannonical_mods'] == 0, "colorby"] = "lightgray"
        de_filtered.loc[de_filtered['num_cannonical_mods'] == 1, "colorby"] = "blue"
        de_filtered.loc[de_filtered['num_cannonical_mods'] > 1, "colorby"] = "red"
    if COLOR_BY == "wam_change":
        de_filtered.loc[
            (de_filtered['wam_change'] > 0)
            , "colorby"] = "red"
    if COLOR_BY == "tes_change":
        de_filtered.loc[
            (de_filtered['tes_change'] > 0)
            , "colorby"] = "red"
    if COLOR_BY == "updown":
        # pval_cutoff = 0.05
        fc_cutoff = 1
        de_filtered.loc[
            (de_filtered['FDR'] < pval_cutoff) & 
            (de_filtered['logFC'] > fc_cutoff)
            , "colorby"] = "red"
        de_filtered.loc[
            (de_filtered['FDR'] < pval_cutoff) & 
            (de_filtered['logFC'] < -fc_cutoff)
            , "colorby"] = "blue"
    if COLOR_BY == "neighbor_type":
        print("hehe")
    if COLOR_BY == "neighbor_methylation_change":
        for a, b in neighbour_file_df['co_directional']:
            if ("MIT" not in a) and ("API" not in a):
                # if a not in de_filtered.index.to_list() or b not in tes_file_df.parent_id.to_list():
                #     print("WARNING: neighbour {} or {} not found in analysis!".format(a, b))
                # else:
                # print("setting {} to {}'s TES change...".format(a, b))
                # neighbor_tes_change = tes_file_df[tes_file_df.parent_id == b].iloc[0].tes_change
                # print("setting {} to {}'s TES change ({})...".format(a, b, neighbor_tes_change))
                
                # a is always upstream
                # b is always downstream

                # neighbor_counts = 0
                # for key in ['co_directional', 'convergent', 'divergent', 'embedded_co_directional', 'embedded_antiparallel']:
                #     if b in neighbor_type_counts[key]:
                #         neighbor_counts += neighbor_type_counts[key][b]
                
                if a in de_filtered.index and b in de_filtered.index:
                    # de_filtered.at[b, 'colorby'] = de_filtered.at[a, 'tes_change']
                    # if de_filtered.at[a, 'tes_change'] > 0:
                    if de_filtered.at[b, 'tes_change'] == 0:
                        de_filtered.at[b, 'colorby'] = 'salmon'
                # de_filtered.at[a, 'is_downstream'] = 'salmon'
    if COLOR_BY == "neighbor_major_minor":
        # color the 
        for a, b in neighbour_file_df['convergent']:
            if ("MIT" not in a) and ("API" not in a):
                if a in de_filtered.index and b in de_filtered.index and (de_filtered.at[a, 'FDR'] < pval_cutoff or de_filtered.at[b, 'FDR'] < pval_cutoff):
                    if de_filtered.at[a, 'logCPM'] > de_filtered.at[b, 'logCPM']:
                        de_filtered.at[a, 'colorby'] = 'green'
                        de_filtered.at[b, 'colorby'] = 'salmon'
                    else:
                        de_filtered.at[a, 'colorby'] = 'salmon'
                        de_filtered.at[b, 'colorby'] = 'green'
                # if a in de_filtered.index and b not in de_filtered.index:
                #     de_filtered.at[b, 'colorby'] = 'green'
                # if a not in de_filtered.index and b in de_filtered.index:
                #     de_filtered.at[b, 'colorby'] = 'green'
    if COLOR_BY == "gene_list":

        goi = []
        with open(genes_to_color_file_path) as file:
            goi = [line.rstrip() for line in file]

        de_filtered.loc[de_filtered.index.isin(goi), 'colorby'] = 'salmon'

        print(goi)

    if COLOR_BY == "gene_list_downstream_codirectional":
        # genes_to_color_file_path = "28hpi_no_overlaps_ctrl_overlaps_ks.txt"
        goi = []
        with open(genes_to_color_file_path) as file:
            goi = [line.rstrip() for line in file]

        for a, b in neighbour_file_df['co_directional']:
            if ("MIT" not in a) and ("API" not in a):
                if a in de_filtered.index and b in de_filtered.index and b in goi:
                    de_filtered.at[b, 'colorby'] = 'salmon'

    background_points = de_filtered[de_filtered["colorby"] == "lightgray"]
    others = de_filtered[de_filtered["colorby"] != "lightgray"]

    print(de_filtered)

    fig, axes = plt.subplots()
    # y_ax_label = '-log10_adj_pval'
    # x_ax_label = 'logFC'
    # axes.set_ylabel('adjusted p values (-log10)')
    # axes.set_xlabel('fold change (log2 KS/C)')
    # log10_pval_cutoff = numpy.log10(0.05)
    # axes.axhline(y=-log10_pval_cutoff, color='grey', ls="--", linewidth=1.0)

    y_ax_label = 'logFC'
    x_ax_label = 'logCPM'
    axes.set_ylabel('fold change (log2 KS/C)')
    axes.set_xlabel('transcript abundance (logCPM)')
    axes.axhline(y=0, color='grey', ls="--", linewidth=1.0)



    axes.scatter(
        background_points[x_ax_label].to_list(), 
        background_points[y_ax_label].to_list(),
        c=background_points['colorby'].to_list(),
        s=10
    )

    axes.scatter(
        others[x_ax_label].to_list(), 
        others[y_ax_label].to_list(),
        c=others['colorby'].to_list(),
        s=10
    )

    axes.set_ylim(ymin=-6, ymax=8)
    axes.set_xlim(xmin=2, xmax=16)

    # # traditional volcano
    axes = de_filtered.plot.scatter(
        x=x_ax_label,
        y=y_ax_label,
        c='colorby',
        s=10
    )

    # mettl3 ks 28,32,36 hpi limits
    axes.set_ylim(ymin=-6, ymax=8)
    axes.set_xlim(xmin=4, xmax=16)

    # axes.set_ylim(ymin=-4, ymax=8)
    # axes.set_xlim(xmin=2, xmax=16)


    def show_label(sel):
        index = sel.index
        print(index)
        if index < len(de_filtered):
            sel.annotation.set_text(de_filtered.index[index])
            print(de_filtered.index[index])
            
    mplcursors.cursor(axes, hover=True).connect("add", show_label)

    plt.show()

