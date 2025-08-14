# read unmapped (0x4)
# read reverse strand (0x10)
# not primary alignment (0x100)
# read fails platform/vendor quality checks (0x200)
# read is PCR or optical duplicate (0x400)
# https://broadinstitute.github.io/picard/explain-flags.html
BAM_UNMAPPED = 0x4
BAM_REVERSE_STRAND = 0x10
BAM_SECONDARY_ALIGNMENT = 0x100
BAM_FAIL_QC = 0x200
BAM_DUPLICATE = 0x400

BAM_PILEUP_DEFAULT_FLAGS = BAM_UNMAPPED | BAM_SECONDARY_ALIGNMENT | BAM_FAIL_QC | BAM_DUPLICATE

d_phred = {}
d_mapq = {}
d_tlen = {}
d_read_ids = {}
dataframes = {}

GFF_DF = None
GFF_PARENT_TREE = {}
ANNOTATION_FILE = None
GFF_DF = None
CLOCKS = {}

TES_SUMMARY_HEADER = ["gene_id", "wart_change", "wart_before", "wart_after", "p_inter_treatment", "p_same_treatment", "tes", "tes_curve_r2", "tes_curve_coeff", "average_expression", "cannonical_mods", "wam_before", "wam_after", "wam_change"]

MODKIT_BEDMETHYL_HEADER = [
    "contig", "start", "end", "code", "score", "strand", 
    "start_2", "end_2", "color", "valid_cov", "percent_mod", "num_mod", 
    "num_canonical", "num_other_mod", "num_delete", "num_fail", "num_diff", "num_nocall"
]

GENERIC_BED_HEADER = [
    "contig",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "ID"
]

FEATURECOUNTS_HEADER = [
    "read_id", "status", "number of targets", "targets"
]

PYSAM_MOD_TUPLES = {
    'm6A_rev': ('A', 1, 'a'),
    'm6A_inosine_rev': ('A', 1, 17596),
    'pseU_rev': ('T', 1, 17802),
    'm5C_rev': ('C', 1, 'm'),
    'm6A_for': ('A', 0, 'a'),
    'm6A_inosine_for': ('A', 0, 17596),
    'pseU_for': ('T', 0, 17802),
    'm5C_for': ('C', 0, 'm')
}