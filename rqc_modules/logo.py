def logo(args):
    import logomaker

    MINUS_OFFSET = int(args.inputs[0]) + 1
    PLUS_OFFSET = int(args.inputs[1]) - 1
    filename = args.inputs[2]

    gff = process_annotation_file()
    gff_df = gff.attributes_to_columns()
    gff_df['strand'] = gff_df['strand'].astype('category')
    gff_df['seq_id'] = gff_df['seq_id'].astype('category')
    gff_df['ID'] = gff_df['ID'].astype('category')
    gff_df['type'] = gff_df['type'].astype('category')

    HAS_LOCUS_TAG = 'locus_tag' in gff_df.columns.to_list()
    if HAS_LOCUS_TAG:
        gff_df['locus_tag'] = gff_df['locus_tag'].astype('category')


    print("GFF types:")
    types = set(gff_df['type'].to_list())
    type_counts = {}
    for t in types:
        of_this_type = gff_df[gff_df.type == t]
        type_counts[t] = len(of_this_type)

    contig_lengths = {}
    if 'region' in types:
        type_rows = gff_df[gff_df.type == 'region']

        for _, row in type_rows.iterrows():
            contig_lengths[row.seq_id] = row.end

    fasta = process_genome_file(contig_lengths)

    site_file_df = pandas.read_csv(filename, index_col=False, sep='\t', header=0, names=GENERIC_BED_HEADER)
    print(site_file_df)
    site_file_df['strand'] = site_file_df['strand'].astype('category')
    site_file_df['contig'] = site_file_df['contig'].astype('category')
    site_file_df['start'] = site_file_df['start'].astype('int32')
    site_file_df['end'] = site_file_df['end'].astype('int32')

    bf_count = len(site_file_df)
    # drop any entries which have contig not in our fasta
    site_file_df = site_file_df[site_file_df['contig'].isin(fasta.keys())]
    print("REMOVED {} DUE TO MISSING CONTIG".format(bf_count - len(site_file_df)))

    seqs = []
    mal_seqs = []

    for row_idx, row in site_file_df.iterrows():
        if (row.start-MINUS_OFFSET) < 0 or (row.end+PLUS_OFFSET) > fasta[row.contig]['length']:
            mal_seqs.append(row)
        else:
            found_seq = fasta[row.contig]['sequence'][row.start-MINUS_OFFSET:row.end+PLUS_OFFSET]

            if row.strand == "-":
                found_seq = reverse_complement(found_seq)
            
            seqs.append(found_seq.upper())
            print("{}: {}:{}-{} ({})".format(found_seq, row.contig, row.start, row.end, row.strand))

    print("NUM MAL SEQS: {}".format(len(mal_seqs)))
    print(mal_seqs)
    counts_matrix = logomaker.alignment_to_matrix(seqs)
    print(counts_matrix)
    l = logomaker.Logo(counts_matrix)
    l.ax.set_xticklabels(["", "-2", "-1", "0", "+1", "+2"])
    l.ax.set_ylabel('count')
    plt.show()
