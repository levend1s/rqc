def find_gene_neighbours(args):
    type = args.inputs[0]

    ANNOTATION_FILE = gffpandas.read_gff3(ANNOTATION_FILE_PATH)
    GFF_DF = ANNOTATION_FILE.attributes_to_columns()
    GFF_DF['ID'] = GFF_DF['ID'].astype('category')
    GFF_DF['type'] = GFF_DF['type'].astype('category')
    GFF_DF['seq_id'] = GFF_DF['seq_id'].astype('category')

    # keep only entries that don't have a parent (removes exons, utrs etc)
    # print(GFF_DF['Parent'])
    GFF_DF = GFF_DF[GFF_DF['Parent'].isna()]

    print(GFF_DF)
    
    if type == 'all':
        all_types = set(GFF_DF['type'])
        print(all_types)
        gff_matching_type = GFF_DF[GFF_DF['type'].isin(all_types)]
    else:
        gff_matching_type = GFF_DF[GFF_DF['type'] == type]

    neighbours_series = [[]] * len(gff_matching_type)
    i = 0

    print("FOUND {} MATCHES FOR TYPE {}".format(len(gff_matching_type), TYPE))

    for a_idx, a in gff_matching_type.iterrows():
        print("processing: {}".format(a['ID']))
        neighbours = []

        same_contig = gff_matching_type[gff_matching_type['seq_id'] == a['seq_id']]
        for b_idx, b in same_contig.iterrows():
            
            #skip checking the same gene
            if (a_idx != b_idx) and (a.seq_id == b.seq_id):

                # do the to genes overlap? add it to the list of neighbours
                if (a.end <= (b.end + NEIGHBOUR_DISTANCE) and a.end >= (b.start - NEIGHBOUR_DISTANCE)) \
                    or (a.start <= (b.end + NEIGHBOUR_DISTANCE) and a.start >= (b.start - NEIGHBOUR_DISTANCE)):
                    
                    neighbours.append(b['ID'])

        print("neighbours: {}".format(neighbours))
        neighbours_series[i] = neighbours
        i += 1

    neighbours_df = gff_matching_type[['ID', 'strand', 'type', 'start', 'end', 'seq_id']].copy()
    neighbours_df['neighbours'] = neighbours_series

    TES_SUMMARY_PATH = "./gene_neighbours.tsv"
    print(neighbours_df)
    neighbours_df.to_csv(TES_SUMMARY_PATH, sep='\t', index=False)

