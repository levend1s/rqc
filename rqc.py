import pysam
import matplotlib.pyplot as plt
import argparse
import numpy

parser = argparse.ArgumentParser(description="filter bam file by qscores / mapqs")
parser.add_argument('function')
parser.add_argument('inputs', nargs='+')
parser.add_argument('-q', '--min_phred', type=int, default=0)
parser.add_argument('-m', '--min_mapq', type=int, default=0)
parser.add_argument('-o', '--outfile', type=str)

args = parser.parse_args()

INPUT = args.inputs
MIN_PHRED = args.min_phred
MIN_MAPQ = args.min_mapq
FUNCTION = args.function
OUTFILE = args.outfile

def plot_qscore_mapq_hists(qscores, mapqs):
    for key in qscores.keys():
        bins = numpy.arange(numpy.floor(qscores[key].min()), numpy.ceil(qscores[key].max()))
        qscore_hist, bin_edges = numpy.histogram(qscores[key], bins=bins, density=True)
        # bin_edges has len(qscore_hist) + 1. We could find the midpoint of each edge, but let's be lazy and take the leftmost edge
        plt.plot(bin_edges[:-1], qscore_hist, label=key)  # arguments are passed to np.histogram

    # plt.hist(mapqs, bins='auto')  # arguments are passed to np.histogram
    # plt.title("MAPQ")
    plt.title("PHRED")
    plt.legend(loc="upper right")
    if (OUTFILE):
        plt.savefig(OUTFILE)
    else:
        plt.show()

d_phred = {}
d_mapq = {}

for i in range(0, len(INPUT), 2):
    label = INPUT[i]
    filename = INPUT[i+1]

    phred_scores = []
    mapq_scores = []

    samfile = pysam.AlignmentFile(filename, 'rb')
    iter = samfile.fetch()

    for x in iter:
        # and not x.is_secondary ???
        if (x.mapping_quality > MIN_MAPQ and x.get_tag('qs') > MIN_PHRED):
            # write alignment to out file
            phred_scores.append(round(x.get_tag('qs'), 2))
            #mapq_scores.append(x.mapping_quality)

    d_phred[label] = numpy.array(phred_scores)
    d_mapq[label] = mapq_scores

    samfile.close()

if FUNCTION == "plot":
    plot_qscore_mapq_hists(d_phred, d_mapq)

