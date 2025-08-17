def gene_neighbour_analysis(args):
    gene_neighbour_tsv_file_path = args.inputs[0]
    print("LOADING: {}".format(gene_neighbour_tsv_file_path))
    gene_neighbour_df = pandas.read_csv(gene_neighbour_tsv_file_path, sep='\t')
    gene_neighbour_df['ID'] = gene_neighbour_df['ID'].astype('category')
    gene_neighbour_df['type'] = gene_neighbour_df['type'].astype('category')
    gene_neighbour_df['seq_id'] = gene_neighbour_df['seq_id'].astype('category')
    gene_neighbour_df['strand'] = gene_neighbour_df['strand'].astype('category')

    gene_neighbour_df["neighbours"] = gene_neighbour_df.neighbours.apply(lambda s: ast.literal_eval(s))

    print("gene_neighbour_df size: {}".format(len(gene_neighbour_df)))

    all_neighbour_pairs = []
    lonely_genes = []

    # returns true if a is greater than b and less than c
    def is_between(a, b, c):
        if (a >= b) and (a <= c):
            return True
        else:
            return False 

    for a_idx, a in gene_neighbour_df.iterrows():
        if len(a['neighbours']) == 0:
            lonely_genes.append(a.ID)

        else:
            for b_id in a['neighbours']:

                smaller = min(a.ID, b_id)
                bigger = max(a.ID, b_id)

                all_neighbour_pairs.append((smaller, bigger))

    unique_neighbour_pairs = set(all_neighbour_pairs)

    d_neighbour_types = {}
    d_neighbour_types['co_directional'] = []
    d_neighbour_types['embedded_antiparallel'] = []
    d_neighbour_types['divergent'] = []
    d_neighbour_types['convergent'] = []
    d_neighbour_types['warning'] = []
    d_neighbour_types['embedded_co_directional'] = []
    d_neighbour_types['lonely'] = lonely_genes

    for pair in unique_neighbour_pairs:
        print(pair)
        a = gene_neighbour_df[gene_neighbour_df['ID'] == pair[0]].iloc[0]
        b = gene_neighbour_df[gene_neighbour_df['ID'] == pair[1]].iloc[0]

        if a.strand == "+":
            a_5p_pos = a.start
            a_3p_pos = a.end
        else:
            a_5p_pos = a.end
            a_3p_pos = a.start

        if b.strand == "+":
            b_5p_pos = b.start
            b_3p_pos = b.end
        else:
            b_5p_pos = b.end
            b_3p_pos = b.start

        # parallel / co-directional
        # if a and b are on same strand
        if (b.strand == a.strand) and is_between(a.start, b.start, b.end) and is_between(a.end, b.start, b.end):
            # embedded co directional (complete)
            print("{} and {}: EMBEDDED CO DIRECTIONAL".format(a.ID, b.ID))
            d_neighbour_types['embedded_co_directional'].append(pair)
        elif (b.strand == a.strand) and is_between(b.start, a.start, a.end) and is_between(b.end, a.start, a.end):
            # embedded co directional (complete)
            print("{} and {}: EMBEDDED CO DIRECTIONAL".format(a.ID, b.ID))
            d_neighbour_types['embedded_co_directional'].append(pair)
        elif b.strand == a.strand:
            print("{} and {}: CO DIRECTIONAL".format(a.ID, b.ID))

            if a.strand == "+":
                d_neighbour_types['co_directional'].append((pair[0], pair[1]))
            else:
                d_neighbour_types['co_directional'].append((pair[1], pair[0]))

        elif is_between(a.start, b.start, b.end) and is_between(a.end, b.start, b.end):
            # embedded anti parallel (complete)
            print("{} and {}: EMBEDDED ANTIPARALLEL".format(a.ID, b.ID))
            d_neighbour_types['embedded_antiparallel'].append(pair)

        elif is_between(b.start, a.start, a.end) and is_between(b.end, a.start, a.end):
            # embedded anti parallel (complete)
            print("{} and {}: EMBEDDED ANTIPARALLEL".format(a.ID, b.ID))
            d_neighbour_types['embedded_antiparallel'].append(pair)

        elif is_between(a_5p_pos, b.start - NEIGHBOUR_DISTANCE, b.end + NEIGHBOUR_DISTANCE):
            # divergent (partial, 5' ends overlap)
            print("{} and {}: DIVERGENT".format(a.ID, b.ID))
            d_neighbour_types['divergent'].append(pair)

        elif is_between(a_3p_pos, b.start - NEIGHBOUR_DISTANCE, b.end + NEIGHBOUR_DISTANCE):
            # convergent (partial, 3' ends overlap)
            print("{} and {}: CONVERGENT".format(a.ID, b.ID))
            d_neighbour_types['convergent'].append(pair)

        else:
            print("WARNING: COULDN'T DETERMINE NEIGHBOUR NATURE OF {} and {}".format(a.ID, b.ID))
            d_neighbour_types['warning'].append(pair)


    outfile = "./gene_neighbour_analysis.json"
    with open(outfile, 'w') as f:
        json.dump(d_neighbour_types, f)

    d_neighbour_types_counts = d_neighbour_types.copy()
    for k, v in d_neighbour_types_counts.items():
        d_neighbour_types_counts[k] = len(v)

    pprint(d_neighbour_types)
    print(d_neighbour_types_counts)

    d_neighbour_types_counts.pop("warning")
    plt.bar(*zip(*d_neighbour_types_counts.items()))
    plt.show()

