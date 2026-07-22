import pandas
import pysam
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import MultiLabelBinarizer, StandardScaler

from rqc_modules.constants import PYSAM_MOD_TUPLES
from rqc_modules.utils import process_input_files, process_annotation_file

def no_values_within(s, x, tol=100):
    for v in s:
        if abs(v - x) <= tol:
            return False
    return True

# scan region and pileup mods (list of genomic positions and their mod / unmod ratios)
# create a table where read_ids are row, columns are mods (m6A, m5C, pseU, m6A_inosine) with a list of mod positions (genomic space)
# other columns include read_start, read_end, read_length, read_strand, poly_A length, average_read_quality (could we cluster by individual base quality?)
# Perform dimensionality reduction (PCA, tSNE, UMAP) on this table and cluster reads based on their mod positions and other features
# create cartoon representation of read type that each cluster represents (e.g. m6A at position 100, m5C at position 150, etc.)
# TODO consider how this will work for bigger regions?

def cluster_transcripts(args):
    INPUT = args.input
    ANNOTATION_FILE = args.annotation
    IDS = args.ids
    COVERAGE_PADDING = args.padding
    MOD_PROB_THRESHOLD = args.mod_prob_threshold
    OUTFILE = args.outfile
    MIN_DELETION_LENGTH = args.min_deletion_length

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)

    # load annotation file and find indexes for all parent children
    gff_df = process_annotation_file(ANNOTATION_FILE)
    if COVERAGE_PADDING:
        gff_df["type"] = gff_df["type"].cat.add_categories(["{}bp".format(COVERAGE_PADDING)])

    matches = gff_df[gff_df['ID'].isin(IDS)]

    if matches.empty:
        print("ERROR: no matches found for ids {}".format(IDS))

    MODS = ['m6A', 'm5C', 'pseU', 'm6A_inosine']

    read_table_header = [
        "read_id",
        "label",
        "bamfile_path",
        "contig",
        "read_start",
        "read_end",
        "read_strand",
        "read_length",
        "poly_a_length",
        "average_quality",
        "introns"
    ]

    for mod in MODS:
        read_table_header.append("{}_positions".format(mod))
        read_table_header.append("{}_num_mods".format(mod))

    read_table = pandas.DataFrame(columns=read_table_header)
    read_table_index = 0

    input_files = process_input_files(INPUT)

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)
    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']

    for label in bam_labels:
        samfile_path = input_files[label]['path']
        print("PROCESSING BAM: {}".format(samfile_path))

        samfile = pysam.AlignmentFile(samfile_path, 'rb')

        for _, row in matches.iterrows():
            READS_IN_REGION = list(samfile.fetch(
                contig=row['seq_id'], 
                start=row['start']-COVERAGE_PADDING, 
                stop=row['end']+COVERAGE_PADDING
            ))

            # filter reads
            for i in range(len(READS_IN_REGION)):
                r = READS_IN_REGION[i]
                if r.is_secondary or r.is_supplementary:
                    continue


                mod_positions = {}

                for mod in MODS:
                    mod_positions["{}_positions".format(mod)] = []
                    if r.is_forward:
                        pysam_mod_tuple_code = '{}_for'.format(mod)
                    else:
                        pysam_mod_tuple_code = '{}_rev'.format(mod)

                    if pysam_mod_tuple_code in PYSAM_MOD_TUPLES:

                        ref_pos = r.get_reference_positions(full_length=True)
                        mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES[pysam_mod_tuple_code])

                        if mods_probs:
                            # keep only mod positions which are above mod prob threshold
                            read_mod_positions = [x[0] for x in mods_probs if x[1] >= PYSAM_MOD_THRESHOLD]
                            genomic_mod_positions = [ref_pos[mod] for mod in read_mod_positions if ref_pos[mod] is not None]

                            mod_positions["{}_positions".format(mod)] = genomic_mod_positions


                # introns
                introns = []
                ref_pos = r.reference_start

                for op, length in r.cigartuples:
                    if op == 2:  # D = 
                        if length >= MIN_DELETION_LENGTH:
                            introns.append((ref_pos, ref_pos + length))

                        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
                            ref_pos += length

                # phred quality
                avg_quality = (
                    sum(r.query_qualities) / len(r.query_qualities)
                    if r.query_qualities
                    else 0
                )                

                poly_a_length = 0

                if r.has_tag('pt:i'):
                    poly_a_length = r.get_tag('pt:i')

                read_strand = '+' if r.is_forward else '-'

                read_entry = [
                    r.query_name,
                    label,
                    samfile_path,
                    row['seq_id'],
                    r.reference_start,
                    r.reference_end,
                    read_strand,
                    r.query_length,
                    poly_a_length,
                    avg_quality,
                    introns,
                ]

                for mod in MODS:
                    read_entry.append(mod_positions["{}_positions".format(mod)])
                    read_entry.append(len(mod_positions["{}_positions".format(mod)]))

                read_table.loc[read_table_index] = read_entry
                read_table_index += 1
            
        # print(read_table)
        samfile.close()

    read_table.to_csv(OUTFILE, sep='\t', index=False)