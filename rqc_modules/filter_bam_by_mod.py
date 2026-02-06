import pysam

from rqc_modules.constants import PYSAM_MOD_TUPLES
from rqc_modules.utils import process_input_files

def no_values_within(s, x, tol=100):
    for v in s:
        if abs(v - x) <= tol:
            return False
    return True

def filter_bam_by_mod(args):
    INPUT = args.input
    MOD_PROB_THRESHOLD = args.mod_prob_threshold
    INCLUDE_SITES = args.include
    CONTIG = args.contig
    STRAND = args.strand
    EXCLUDE_TOL = args.exclude_tol
    MOD = args.mod # m6A, m5C, pseU, m6A_inosine

    print("INCLUDE SITES: {}".format(INCLUDE_SITES))

    input_files = process_input_files(INPUT)

    PYSAM_MOD_THRESHOLD = int(256 * MOD_PROB_THRESHOLD)
    bam_labels = [l for l in input_files.keys() if input_files[l]['type'] == 'bam']

    for label in bam_labels:
        samfile_path = input_files[label]['path']
        print("PROCESSING BAM: {}".format(samfile_path))

        samfile = pysam.AlignmentFile(samfile_path, 'rb')
        READS_IN_REGION = list(samfile.fetch(
            contig=CONTIG, 
            start=INCLUDE_SITES-100, 
            stop=INCLUDE_SITES+100
        ))

        read_indexes_to_process = []

        if STRAND == '+':
            pysam_mod_tuple_code = '{}_for'.format(MOD)
        else:
            pysam_mod_tuple_code = '{}_rev'.format(MOD)

        read_indexes_to_process = []
        read_indexes_to_process_exclude = []
        # filter reads
        for i in range(len(READS_IN_REGION)):
            r = READS_IN_REGION[i]

            # keep only reads in the same direction as this strand
            if (STRAND == "+" and r.is_forward) or (STRAND == "-" and r.is_reverse):

                ref_pos = r.get_reference_positions(full_length=True)
                mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES[pysam_mod_tuple_code])

                if mods_probs:
                    # keep only mod positions which are above mod prob threshold
                    read_mod_positions = [x[0] for x in mods_probs if x[1] >= PYSAM_MOD_THRESHOLD]
                    genomic_mod_positions = [ref_pos[mod] for mod in read_mod_positions if ref_pos[mod] is not None]

                    # if r.query_name == "588c2c93-f1c8-423e-8c1f-a850f8a5f495":
                    #     print(genomic_mod_positions)

                    if INCLUDE_SITES in genomic_mod_positions:
                        read_indexes_to_process.append(i)

                    if (INCLUDE_SITES <= r.reference_end and INCLUDE_SITES >= r.reference_start) and no_values_within(genomic_mod_positions, INCLUDE_SITES, tol=EXCLUDE_TOL):
                        read_indexes_to_process_exclude.append(i)
                else:
                    if (INCLUDE_SITES <= r.reference_end and INCLUDE_SITES >= r.reference_start):
                        read_indexes_to_process_exclude.append(i)

        filtered_bam_filename = "{}_{}.bam".format(label, str(INCLUDE_SITES))
        filtered_bam_filename_exclude = "{}_{}_exclude.bam".format(label, str(INCLUDE_SITES))

        print("GENERATING FILTERED BAM: {}".format(filtered_bam_filename))

        rqc_filtered = pysam.AlignmentFile(filtered_bam_filename, "wb", template=samfile)
        for i in read_indexes_to_process:
            rqc_filtered.write(READS_IN_REGION[i])

        rqc_filtered.close()

        rqc_filtered_exclude = pysam.AlignmentFile(filtered_bam_filename_exclude, "wb", template=samfile)
        for i in read_indexes_to_process_exclude:
            rqc_filtered_exclude.write(READS_IN_REGION[i])

        rqc_filtered_exclude.close()
        samfile.close()