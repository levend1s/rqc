import re
import pandas
from pprint import pprint

from rqc_modules.utils import process_annotation_file, process_genome_file, reverse_complement
from rqc_modules.constants import GENERIC_BED_HEADER

def motif_finder(args):
    MOTIF = args.motif
    OUTFILE = args.output
    ANNOTATION_FILE = args.annotation
    GENOME_FILE = args.genome
    FEATURE_FILTER = None

    if args.feature_filter:
        FEATURE_FILTER = args.feature_filter

    print("LOG - motif finder started")
    print("LOG - motif: {}".format(MOTIF))
    print("LOG - annotation file: {}".format(ANNOTATION_FILE))
    print("LOG - genome file: {}".format(GENOME_FILE))
    print("LOG - output file: {}".format(OUTFILE))

    # if filter by genomic regions
    # filter=exons
    # filter=first_exon
    # filter=last_exon
    gff = process_annotation_file(ANNOTATION_FILE)
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

    fasta = process_genome_file(GENOME_FILE, contig_lengths)

    MOTIF_RC = reverse_complement(args.motif)
    print("LOG - looking for motif: {} (rc={})".format(MOTIF, MOTIF_RC))
    forward_lookahead_regex = re.compile("(?=({}))".format(MOTIF), re.IGNORECASE)
    reverse_lookahead_regex = re.compile("(?=({}))".format(MOTIF_RC), re.IGNORECASE)

    # create output file
    rows = []


    # build a dict containing the number of CDS for each gene, needed so we can number our CDS
    lt_cds_counts = {}
    if HAS_LOCUS_TAG:
        uniq_locus_tags = set(gff_df.locus_tag.to_list())
        print("LOG - locus tag in GFF file, finding number of CDS for each gene ({})...".format(len(uniq_locus_tags)))
        for lt in uniq_locus_tags:
            cds_this_gene = gff_df[
                (gff_df.locus_tag == lt) & 
                (gff_df.type == "CDS")]

            if len(cds_this_gene) > 0 and cds_this_gene.iloc[0].strand == "+":
                lt_cds_counts[lt] = 1
            else:
                lt_cds_counts[lt] = len(cds_this_gene)

    for contig in fasta.keys():
        print("LOG - searching: {}".format(contig))

        forward_count = 0
        reverse_count = 0

        # if feature filter the output is different
        # the contig is the parent gene ID
        if FEATURE_FILTER:
            this_contig_genes = gff_df[
                (gff_df.seq_id == contig) & 
                (gff_df.type == FEATURE_FILTER)
            ]

            for _, row in this_contig_genes.iterrows():
                if row.strand == "+":
                    this_regex = forward_lookahead_regex
                    if HAS_LOCUS_TAG:
                        this_cds_idx = lt_cds_counts[row.locus_tag]
                        lt_cds_counts[row.locus_tag] += 1
                else:
                    this_regex = reverse_lookahead_regex
                    if HAS_LOCUS_TAG:
                        this_cds_idx = lt_cds_counts[row.locus_tag]
                        lt_cds_counts[row.locus_tag] -= 1
                    
                if row.phase == ".":
                    phase = 0
                else:
                    phase = int(row.phase)

                # HACK for CryptoBGF, Toxo ME49
                # ID is a different naming scheme and CDS aren't numbered
                if HAS_LOCUS_TAG:
                    row_id = "{}-CDS{}".format(row.locus_tag, this_cds_idx)
                else:
                    row_id = row.ID


                # gff files are 1-indexed, so subtract 1 from the start co-ord for accurate python list slicing
                matches_in_gene = re.finditer(this_regex, fasta[contig]['sequence'][row.start-1:row.end])

                for m in matches_in_gene:
                    # if row.locus_tag == "TGME49_222160":
                    #     print(fasta[contig]['sequence'][row.start-1:row.end])
                    #     print(this_regex)
                    #     print(m.group(1))

                    match = m.group(1)
                    if row.strand == "-":
                        match = reverse_complement(m.group(1))
                
                    # add 1 since an index of 0 in the fasta substring is actually a 1
                    row_summary = [contig, row.start + m.start(), row.start + m.start()+len(m.group(1)), match, 0, row.strand, row_id]

                    rows.append(row_summary)

                    if row.strand == "+":
                        forward_count += 1
                    else:
                        reverse_count += 1
        else:
            forward_matches = re.finditer(forward_lookahead_regex, fasta[contig]['sequence'])
            reverse_matches = re.finditer(reverse_lookahead_regex, fasta[contig]['sequence'])
            
            strand = "+"
            for m in forward_matches:
                row_summary = [contig, m.start()+1, m.start()+len(m.group(1))+1, m.group(1), 0, strand, ""]
                rows.append(row_summary)
                forward_count += 1

            strand = "-"
            for m in reverse_matches:
                row_summary = [contig, m.start()+1, m.start()+len(m.group(1))+1, reverse_complement(m.group(1)), 0, strand, ""]
                rows.append(row_summary)
                reverse_count += 1

        print("LOG - {} matches: {} forward, {} reverse".format(contig, forward_count, reverse_count))

    motif_matches_df = pandas.DataFrame(columns=GENERIC_BED_HEADER, data=rows)

    motif_matches_df.to_csv(OUTFILE, sep='\t', index=False)

    bases_counts = {
        'A': 0,
        'T': 0,
        'C': 0,
        'G': 0
    }

    genome_size = 0

    for contig in fasta.keys():
        genome_size += len(fasta[contig]['sequence'])
        for b in bases_counts.keys():
            bases_counts[b] += fasta[contig]['sequence'].upper().count(b)

    print("genome size: {}".format(genome_size))
    pprint(bases_counts)
    