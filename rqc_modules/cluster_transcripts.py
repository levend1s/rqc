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
# TODO differential analysis of clusters between conditions (wt vs kd)
# TODO include large indels (ie introns) in the clustering analysis

def cluster_transcripts(args):
    INPUT = args.input
    ANNOTATION_FILE = args.annotation
    IDS = args.ids
    COVERAGE_PADDING = args.padding
    MOD_PROB_THRESHOLD = args.mod_prob_threshold
    BAMFILE = args.bamfile
    OUTFILE = args.outfile

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)

    # load annotation file and find indexes for all parent children
    gff_df = process_annotation_file(ANNOTATION_FILE)
    if COVERAGE_PADDING:
        gff_df["type"] = gff_df["type"].cat.add_categories(["{}bp".format(COVERAGE_PADDING)])

    matches = gff_df[gff_df['ID'].isin(IDS)]

    if matches.empty:
        print("ERROR: no matches found for ids {}".format(IDS))

    read_table_header = [
        "read_id",
        "source_file",
        "read_start",
        "read_end",
        "read_strand",
        "read_length",
        "poly_a_length",
        "mod_positions",
        "num_mods"
    ]
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

            # TODO fix
            MODS = ['m6A', 'm5C', 'pseU', 'm6A_inosine']
            MOD = 'm6A'

            # filter reads
            for i in range(len(READS_IN_REGION)):
                r = READS_IN_REGION[i]
                if r.is_forward:
                    pysam_mod_tuple_code = '{}_for'.format(MOD)
                else:
                    pysam_mod_tuple_code = '{}_rev'.format(MOD)
                
                ref_pos = r.get_reference_positions(full_length=True)
                mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES[pysam_mod_tuple_code])

                genomic_mod_positions = []
                num_mods = 0
                poly_a_length = 0

                if mods_probs:
                    # keep only mod positions which are above mod prob threshold
                    read_mod_positions = [x[0] for x in mods_probs if x[1] >= PYSAM_MOD_THRESHOLD]
                    genomic_mod_positions = [ref_pos[mod] for mod in read_mod_positions if ref_pos[mod] is not None]
                    num_mods = len(genomic_mod_positions)

                if r.has_tag('pt:i'):
                    poly_a_length = r.get_tag('pt:i')

                read_strand = '+' if r.is_forward else '-'

                read_entry = [
                    r.query_name,
                    label,
                    r.reference_start,
                    r.reference_end,
                    read_strand,
                    r.query_length,
                    poly_a_length,
                    genomic_mod_positions,
                    num_mods
                ]

                read_table.loc[read_table_index] = read_entry
                read_table_index += 1
            
        # print(read_table)
        samfile.close()

    read_table.to_csv(OUTFILE, sep='\t', index=False)