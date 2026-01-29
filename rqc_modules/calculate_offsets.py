import pandas
from rqc_modules.constants import GENERIC_BED_HEADER_BASE, GENERIC_BED_HEADER

def calculate_offsets(args):
    REFERENCE_BED = args.reference
    DISTANCE = args.distance
    OUTFILE = args.output
    SAME_STRAND = args.same_strand
    INPUTS = args.inputs
    BED_REFERENCE = args.bed_reference

    if len(INPUTS) % 2 != 0:
        raise ValueError("Inputs must be provided in pairs: name and file path.")

    print("LOG - calculate_offsets started")
    print("LOG - reference bed file: {}".format(REFERENCE_BED))
    print("LOG - distance: {}".format(DISTANCE))
    print("LOG - output file: {}".format(OUTFILE))
    print("LOG - inputs: {}".format(INPUTS))
    print("LOG - same strand: {}".format(SAME_STRAND))

    d_offset_files_by_contig = {}
    d_coverages = {}
    d_num_pam_sites = {}
    print(INPUTS)
    
    i = 0
    while i < len(INPUTS):
        key = INPUTS[i]
        file_path = INPUTS[i+1]
        print("LOADING: {}".format(file_path))
        file_df = pandas.read_csv(file_path, sep='\t', header=None)

        if file_df.iloc[0][0] != "contig":
            file_df = pandas.read_csv(file_path, sep='\t', header=None)
            file_df.columns = GENERIC_BED_HEADER_BASE
        else:
            file_df = pandas.read_csv(file_path, sep='\t')


        file_df['contig'] = file_df['contig'].astype('category')
        file_df['strand'] = file_df['strand'].astype('category')

        # d_offset_files[key] = file_df
        if key not in d_num_pam_sites:
            d_num_pam_sites[key] = []

        for contig in set(file_df.contig.to_list()):
            # if "KE" in contig:
                # continue

            if contig not in d_offset_files_by_contig:
                d_offset_files_by_contig[contig] = {}
                d_offset_files_by_contig[contig][key] = file_df[file_df.contig == contig]
            else:
                d_offset_files_by_contig[contig][key] = file_df[file_df.contig == contig]

        i += 2

    # HACK for malformed bed file that my pipeline produces
    LAZY = False
    if LAZY:
        special_bed_header = GENERIC_BED_HEADER + ["base_ID"]
        reference_df = pandas.read_csv(REFERENCE_BED, sep='\t', names=special_bed_header, header=None)
    else:
        reference_df = pandas.read_csv(REFERENCE_BED, sep='\t')
    
    if BED_REFERENCE:
        reference_df = pandas.read_csv(REFERENCE_BED, sep='\t', names=GENERIC_BED_HEADER_BASE, header=None)
        reference_df["ID"] = reference_df[["contig", "strand", "start"]].astype(str).agg("_".join, axis=1)



    reference_df['contig'] = reference_df['contig'].astype('category')
    reference_df['strand'] = reference_df['strand'].astype('category')

    # TEST_LIMIT = 1000
    for row_index, row in reference_df.iterrows():
        # if row_index > TEST_LIMIT:
            # break

        if row.strand == "+":
            offset_point = int(row.start)
        else:
            offset_point = int(row.end)

        print("processing row {} (of {})...".format(row_index, len(reference_df)))

        # TODO; this should not be necessary? why are we bound by IDs...
        d_coverages[row.ID] = {}

        d_coverages[row.ID]["position"] = offset_point

        # for k, v in d_offset_files.items():
        for k, v in d_offset_files_by_contig[row.contig].items():
            if SAME_STRAND:
                if row.strand == "+":
                    in_range = v[(v.start >= (offset_point - DISTANCE))
                                & (v.start <= (offset_point + DISTANCE))
                                & (v.contig == row.contig)
                                & (v.strand == row.strand)]
                    
                    offsets = [x - offset_point for x in in_range.start.to_list()]
                else:
                    in_range = v[(v.end >= (offset_point - DISTANCE))
                                & (v.end <= (offset_point + DISTANCE))
                                & (v.contig == row.contig)
                                & (v.strand == row.strand)]
                    
                    offsets = [x - offset_point for x in in_range.end.to_list()]
                    offsets = [-x for x in offsets]
            else:
                in_range = v[(v.end >= (offset_point - DISTANCE))
                                & (v.start <= (offset_point + DISTANCE))]

                for_in_range = in_range[in_range.strand == "+"]
                rev_in_range = in_range[in_range.strand == "-"]

                for_offsets = [x - offset_point for x in for_in_range.start.to_list()]
                rev_offsets = [x - offset_point for x in rev_in_range.end.to_list()]

                offsets = for_offsets + rev_offsets

                if row.strand == "-":
                    offsets = [-x for x in offsets]
                
            # print("offset_point: {}".format(offset_point))
            # print(in_range)
            # print(offsets)
            d_coverages[row.ID][k] = offsets
            d_coverages[row.ID]["{}_count".format(k)] = len(offsets)
            d_num_pam_sites[k].append(len(offsets))
            # if (len(offsets)) > 70:
                # print(row)

        # if row_index > 100:
            # break

    # TODO: write coverages as tsv file for plotting without recalculating offsets
    df = pandas.DataFrame.from_dict(d_coverages, orient='index')
    df = df.rename_axis('gene_id', axis='index')
    df.to_csv(OUTFILE, sep='\t')
    print(df)
    # gene_id genomic_position file_1_count file_1_list_of_positions file_2_count file_2_list_of_positions

