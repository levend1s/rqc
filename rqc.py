
import argparse
import numpy
import sys

numpy.set_printoptions(threshold=sys.maxsize)
numpy.seterr(divide='ignore', invalid='ignore')

from rqc_modules import motif_finder
from rqc_modules import calculate_offsets
from rqc_modules import plot_relative_distance
from rqc_modules import plot_coverage
from rqc_modules import sequence_logo
from rqc_modules import approximate_tes
from rqc_modules import gene_methylation_analysis


def build_parser():
    parser = argparse.ArgumentParser(
        prog="rqc",
        description="RQC – a lightweight RNA-seq tool for transcript level analysis."
    )

    subparsers = parser.add_subparsers(
        title="commands",
        dest="command"
    )

    # ---- motif_finder command ----
    motif_finder_parser = subparsers.add_parser("motif_finder", help="find motifs in genome")
    motif_finder_parser.add_argument("-a", "--annotation", required=True, help="annotation file (GFF)")
    motif_finder_parser.add_argument("-g", "--genome", required=True, help="genome file (FASTA)")
    motif_finder_parser.add_argument("-o", "--output", required=True, help="output file (tsv)")
    motif_finder_parser.add_argument("-m", "--motif", required=True, help="motif to search for")
    motif_finder_parser.add_argument("-f", "--feature_filter", required=False, help="If provided, filter the annotation by feature type (e.g., 'exon', 'CDS', 'gene'). If not provided, no filtering is applied.")

    motif_finder_parser.set_defaults(func=motif_finder.motif_finder)

    # ---- calculate_offsets command ----
    calculate_offsets_parser = subparsers.add_parser("calculate_offsets", help="given files ")
    calculate_offsets_parser.add_argument("-r", "--reference", required=True, help="Bed file with list of coords used as offset point")
    calculate_offsets_parser.add_argument("-d", "--distance", required=True, type=int, help="distance (±) to calculate offsets for")
    calculate_offsets_parser.add_argument("-o", "--output", required=True, help="output file (tsv)")

    calculate_offsets_parser.add_argument("-i", "--inputs", required=True, nargs="*", help="input name and files (e.g., 'input1 file1.bed input2 file2.bed')")
    calculate_offsets_parser.add_argument("-s", "--same_strand", action="store_true", help="If provided, consider motifs on the same strands. If not provided, consider motifs on both strands.")

    calculate_offsets_parser.set_defaults(func=calculate_offsets.calculate_offsets)

    # ---- plot_relative_offsets command ----
    plot_relative_distance_parser = subparsers.add_parser("plot_relative_distance", help="plot relative offsets from a file")
    plot_relative_distance_parser.add_argument("-d", "--distance", required=True, type=int, help="distance (±) to plot offsets for")
    plot_relative_distance_parser.add_argument("-l", "--label", required=True, help="label for the plot")
    plot_relative_distance_parser.add_argument("-i", "--input", required=True, help="input file with offsets")
    plot_relative_distance_parser.add_argument("-o", "--output", required=False, help="output file suffix (e.g., 'plot.png', 'plot.pdf'). If not provided, will not save plot to file.")

    plot_relative_distance_parser.set_defaults(func=plot_relative_distance.plot_relative_distance)

    # ---- plot_coverage command ----
    plot_coverage_parser = subparsers.add_parser("plot_coverage", help="plot relative offsets from a file")
    plot_coverage_parser.add_argument("-i", "--input", required=True, help="input file listing coverage data")
    plot_coverage_parser.add_argument("-a", "--annotation", required=True, help="annotation file (GFF)")
    plot_coverage_parser.add_argument("-m", "--mode", required=True, help="input file listing coverage data")
    plot_coverage_parser.add_argument("-b", "--bins", required=True, type=int, help="input file listing coverage data")
    plot_coverage_parser.add_argument("-d", "--read_depth", required=True, type=int, help="input file listing coverage data")
    
    # at least one of these is required
    plot_coverage_parser.add_argument("--ids", required=False, nargs="*", help="input file listing coverage data")
    plot_coverage_parser.add_argument("--type", required=False, help="input file listing coverage data")

    # optional arguments
    plot_coverage_parser.add_argument("-v", "--verbose", action="store_true", help="verbose mode, benchmarking and printing additional information")
    plot_coverage_parser.add_argument("--plot_density", action="store_true", help="plot density of coverage (expensive)")
    plot_coverage_parser.add_argument("--separate_y_axes", action="store_true", help="plot density of coverage (expensive)")
    plot_coverage_parser.add_argument("--skip_malannotations", action="store_true", help="plot density of coverage (expensive)")

    plot_coverage_parser.add_argument("-c", "--coverage_method", required=False, default="max", help="How to calculate the coverage for each bin: sum, average, max")
    plot_coverage_parser.add_argument("-n", "--mod_normalisation", required=False, default="raw", help="input file listing coverage data")
    plot_coverage_parser.add_argument("-o", "--output", required=False, help="output file suffix (e.g., 'plot.png', 'plot.pdf'). If not provided, will not save plot to file.")
    plot_coverage_parser.add_argument("-p", "--padding", required=False, type=int, default=0, help="input file listing coverage data")
    plot_coverage_parser.add_argument("-r", "--padding_ratio", required=False, type=float, default=0.0, help="input file listing coverage data")
    plot_coverage_parser.add_argument("--line_width", required=False, type=int, default=1, help="input file listing coverage data")

    plot_coverage_parser.set_defaults(func=plot_coverage.plot_coverage)

    # ---- approximate_tes command ----
    approximate_tes_parser = subparsers.add_parser("approximate_tes", help="Based of RNA-seq data and a GFF file, approximate the transcription end sites (TES) of genes.")
    approximate_tes_parser.add_argument("-a", "--annotation", required=True, help="annotation file (GFF)")
    approximate_tes_parser.add_argument("-i", "--input", required=True, help="input file listing coverage data")

    approximate_tes_parser.add_argument("--adjust", required=False, type=int, default=0, help="adjust if the bed files are indexed funny, shifts indexes by this amount")
    approximate_tes_parser.add_argument("-o", "--output", required=False, help="output file suffix (e.g., 'plot.png', 'plot.pdf'). If not provided, will not save plot to file.")
    approximate_tes_parser.add_argument("-p", "--padding", required=False, type=int, default=0, help="input file listing coverage data")
    approximate_tes_parser.add_argument("-d", "--read_depth", required=False, type=int, default=20, help="input file listing coverage data")
    approximate_tes_parser.add_argument("--poly_a_filter", required=False, type=int, default=0, help="input file listing coverage data")
    approximate_tes_parser.add_argument("--filter_for_m6A", required=False, nargs="*", type=int, help="filter for only reads with m6A modification at these positions")
    approximate_tes_parser.add_argument("--filter_out_m6A", required=False, nargs="*", type=int, help="filter for only reads without m6A modification at these positions")
    approximate_tes_parser.add_argument("-v", "--verbose", action="store_true", help="verbose mode, benchmarking and printing additional information")
    approximate_tes_parser.add_argument("--feature_counts", action="store_true", help="if provided filter reads which featureCounts has mapped to a gene")
    
    approximate_tes_parser.add_argument("--exclude_contigs", required=False, nargs="*", help="input file listing coverage data.")
    approximate_tes_parser.add_argument("--compare_apa_between_treatments", required=False, nargs="*", help="filter for only reads without m6A modification at these positions")


    # at least one of these is required
    approximate_tes_parser.add_argument("--ids", required=False, nargs="*", help="input file listing coverage data")
    approximate_tes_parser.add_argument("--type", required=False, help="input file listing coverage data")

    approximate_tes_parser.set_defaults(func=approximate_tes.approximate_tes)

    # ---- gene_methylation_analysis command ----
    # My gut feeling is that mapping mod sites from a bedmethyl to a gff is not the best way, since gff annotation is poor and may not include readthrough sites.
    # I think it would be better to go through each gene, use my (nanopore specific) method of finding reads associated with a gene, then do mod analysis on those reads.    
    gene_methylation_analysis_parser = subparsers.add_parser("gene_methylation_analysis", help="Based of RNA-seq data and a GFF file, approximate the transcription end sites (TES) of genes.")
    gene_methylation_analysis_parser.add_argument("-i", "--input", required=True, help="input file listing coverage data.")
    gene_methylation_analysis_parser.add_argument("-o", "--output", required=False, help="output file suffix (e.g., 'plot.png', 'plot.pdf'). If not provided, will not save plot to file.")

    gene_methylation_analysis_parser.add_argument("--type", required=False, help="input file listing coverage data.")
    gene_methylation_analysis_parser.add_argument("--ids", required=False, nargs="*", help="input file listing coverage data.")
    gene_methylation_analysis_parser.add_argument("--exclude_contigs", required=False, nargs="*", help="input file listing coverage data.")

    gene_methylation_analysis_parser.add_argument("-a", "--annotation", required=True, help="annotation file (GFF)")
    gene_methylation_analysis_parser.add_argument("--poly_a_filter", required=False, type=int, default=0, help="input file listing coverage data")
    gene_methylation_analysis_parser.add_argument("-p", "--padding", required=False, type=int, default=0, help="input file listing coverage data")
    gene_methylation_analysis_parser.add_argument("-r", "--mod_ratio", required=False, type=float, default=0.5, help="input file listing coverage data")
    gene_methylation_analysis_parser.add_argument("-c", "--mod_prob_threshold", required=False, type=float, default=0.95, help="input file listing coverage data")

    gene_methylation_analysis_parser.add_argument("-d", "--read_depth", required=False, type=int, default=20, help="input file listing coverage data")
    gene_methylation_analysis_parser.add_argument("-v", "--verbose", action="store_true", help="verbose mode, benchmarking and printing additional information")
    gene_methylation_analysis_parser.add_argument("-g", "--genome", required=False, help="genome file (FASTA)")

    gene_methylation_analysis_parser.add_argument("--compare_methylation_between_treatments", required=False, nargs="*", help="filter for only reads without m6A modification at these positions")

    gene_methylation_analysis_parser.set_defaults(func=gene_methylation_analysis.gene_methylation_analysis)

    # ---- find_gene_neighbours command ----

    # ---- gene_neighbour_analysis command ----

    # ---- filter_bam_by_mod command ----

    # ---- m6A_specific_tes_analysis command ----

    # ---- logo command ----

    sequence_logo_parser = subparsers.add_parser("sequence_logo", help="find motifs in genome")
    sequence_logo_parser.add_argument("-l", "--length", required=True, type=int, help="distancce to search for motif (±)")
    sequence_logo_parser.add_argument("-o", "--output", required=False, help="output file to save image to")
    sequence_logo_parser.add_argument("-i", "--input", required=True, help="input bed file with genomic coords")
    sequence_logo_parser.add_argument("-a", "--annotation", required=True, help="annotation file (GFF)")
    sequence_logo_parser.add_argument("-g", "--genome", required=True, help="genome file (FASTA)")
    sequence_logo_parser.add_argument("--adjust", required=False, type=int, default=0, help="adjust if the bed files are indexed funny, shifts indexes by this amount")

    sequence_logo_parser.add_argument("-p", "--padding_distance", required=False, type=int, default=0, help="distance to pad the motif sequence (default: 2)")

    sequence_logo_parser.set_defaults(func=sequence_logo.sequence_logo)


    return parser

def main():
    parser = build_parser()
    args = parser.parse_args()

    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(1)

    args.func(args)

if __name__ == "__main__":
    main()


