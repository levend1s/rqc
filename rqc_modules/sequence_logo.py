import logomaker
import pandas
import matplotlib.pyplot as plt
import numpy
import os


from rqc_modules.utils import process_annotation_file, process_genome_file, reverse_complement
from rqc_modules.constants import GENERIC_BED_HEADER_BASE

def sequence_logo(args):
    PADDING_DISTANCE = args.padding_distance
    LENGTH = args.length
    GENOME_FILE = args.genome
    ANNOTATION_FILE = args.annotation
    INPUT = args.input
    ADJUST = args.adjust
    OUTPUT = args.output

    gff_df = process_annotation_file(ANNOTATION_FILE)

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

    site_file_df = pandas.read_csv(INPUT, index_col=False, sep='\t', header=0, names=GENERIC_BED_HEADER_BASE)
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

    for _, row in site_file_df.iterrows():
        if (row.start-PADDING_DISTANCE) < 0 or (row.end+PADDING_DISTANCE) > fasta[row.contig]['length']:
            mal_seqs.append(row)
        else:
            found_seq = fasta[row.contig]['sequence'][row.start-PADDING_DISTANCE+ADJUST:row.end+PADDING_DISTANCE+ADJUST]

            if row.strand == "-":
                found_seq = reverse_complement(found_seq)

            seqs.append(found_seq.upper())
            print("{}: {}:{}-{} ({})".format(found_seq, row.contig, row.start, row.end, row.strand))
                
    print("NUM MAL SEQS: {}".format(len(mal_seqs)))
    print(mal_seqs)

    plasmo_bg = {"A":0.4, "C":0.1, "G":0.1, "T":0.4}
    plasmo_bg_values = numpy.array([plasmo_bg[nt] for nt in 'ACGT'])
    print(plasmo_bg_values)
    counts_matrix = logomaker.alignment_to_matrix(seqs, to_type='probability')

    # CHATGPT START
    H_max = -(plasmo_bg_values * numpy.log2(plasmo_bg_values)).sum()
    # Observed entropy per position
    entropy = -(counts_matrix * numpy.log2(counts_matrix + 1e-9)).sum(axis=1)
    # Information content per position (bits)
    R = H_max - entropy


    bits_mat = counts_matrix.multiply(R, axis=0)

    # Make the logo
    logo = logomaker.Logo(bits_mat)
    logo.ax.set_ylabel("Bits")

    # CHATGPT END

    OUTPUT_FORMAT = OUTPUT.split(".")[-1] if OUTPUT else "png"

    if OUTPUT:
        output_path = os.path.abspath(OUTPUT)

        # Split into directory + basename
        outdir, base = os.path.split(output_path)
        # Prepend prefix to the filename
        new_base_bits = f"logo_bits_{base}"
        new_base_counts = f"logo_counts_{base}"

        # Build final path
        plt.savefig(os.path.join(outdir, new_base_bits), transparent=True, dpi=300, format=OUTPUT_FORMAT)

    counts_matrix = logomaker.alignment_to_matrix(seqs, to_type='counts')

    l = logomaker.Logo(counts_matrix)
    labels = ["-{}".format(x) for x in range(-PADDING_DISTANCE, 0)] + (["0"] * LENGTH) + ["+{}".format(x) for x in range(1, PADDING_DISTANCE+1)]

    print(labels)
    l.ax.set_xticklabels([""] + labels)
    l.ax.set_ylabel('count')

    if OUTPUT:
        plt.savefig(os.path.join(outdir, new_base_counts), transparent=True, dpi=300, format=OUTPUT_FORMAT)


    plt.show()
