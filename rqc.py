
import argparse
import numpy
import sys

numpy.set_printoptions(threshold=sys.maxsize)
numpy.seterr(divide='ignore', invalid='ignore')

from rqc_modules import motif_finder
from rqc_modules import calculate_offsets
from rqc_modules import plot_relative_distance


def build_parser():
    parser = argparse.ArgumentParser(
        prog="rqc",
        description="RQC – a RNA-seq tool for transcript level analysis."
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


