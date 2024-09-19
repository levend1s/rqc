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

# this calculates the NX for a reverse sorted list of read lengths
# You might use this to calculate the N50 or N90, to find the read 
# length which at least 50% or 90% of total nucleotides read belong to
def getNX(reverse_sorted_read_lengths, NX):
    csum = numpy.cumsum(reverse_sorted_read_lengths)
    total_bases_read = sum(reverse_sorted_read_lengths)
    NX_times_total_reads = int(total_bases_read * NX)

    idx = 0
    for i in reversed(range(len(csum))):
        if csum[i] > NX_times_total_reads:
            idx = i
            continue
        else:
            break

    return reverse_sorted_read_lengths[idx]

def print_tlen_distribution(tlens):
    d = []

    for key in tlens.keys():
        cur_tlens = tlens[key]
        mean = int(numpy.mean(cur_tlens))
        median = int(numpy.median(cur_tlens))
        q1, q3 = numpy.quantile(cur_tlens, [0.25, 0.75])
        min_read = int(min(cur_tlens))
        max_read = int(max(cur_tlens))

        # calculate n50
        sorted_tlens = sorted(cur_tlens, reverse=True)
        total_bases_read = sum(sorted_tlens)

        n10 = getNX(sorted_tlens, 0.1)
        n50 = getNX(sorted_tlens, 0.5)
        n90 = getNX(sorted_tlens, 0.9)


        a = [100, 70, 60, 50, 50, 40, 30]

        summary = {
            "sample": key,
            "count": len(cur_tlens),
            "total_bases_read": total_bases_read,
            "mean": mean,
            "median": median,
            "Q1": q1,
            "Q3": q3,
            "min": min_read,
            "max": max_read,
            "N10": n10,
            "N50": n50,
            "N90": n90
        }

        d.append(summary)

    # print summaries
    header = d[0].keys()
    rows = [r.values() for r in d]

    for x in header:
        print("{}\t".format(x), end="")

    print()

    for row in rows:
        for x in row:
            print("{}\t".format(x), end="")
        print()



def plot_qscore_hists(qscores):
    for key in qscores.keys():
        bins = numpy.arange(numpy.floor(qscores[key].min()), numpy.ceil(qscores[key].max()))
        qscore_hist, bin_edges = numpy.histogram(qscores[key], bins=bins, density=True)
        # bin_edges has len(qscore_hist) + 1. We could find the midpoint of each edge, but let's be lazy and take the leftmost edge
        plt.plot(bin_edges[:-1], qscore_hist, label=key)  # arguments are passed to np.histogram

    plt.title("PHRED")
    plt.legend(loc="upper right")
    if (OUTFILE):
        plt.savefig(OUTFILE)

def plot_mapq_hists(mapqs):
    plt.figure()
    plt.title("MAPQ")
    plt.hist(mapqs.values(), histtype="bar", label=mapqs.keys(), density=True)
    plt.legend(loc="upper right")


if FUNCTION == "plot":
    d_phred = {}
    d_mapq = {}
    d_tlen = {}


    for i in range(0, len(INPUT), 2):
        label = INPUT[i]
        filename = INPUT[i+1]

        phred_scores = []
        mapq_scores = []
        template_lengths = []


        samfile = pysam.AlignmentFile(filename, 'rb')
        iter = samfile.fetch()

        for x in iter:
            # and not x.is_secondary ???
            if (x.mapping_quality >= MIN_MAPQ and x.get_tag('qs') >= MIN_PHRED):
                phred_scores.append(round(x.get_tag('qs'), 2))
                mapq_scores.append(x.mapping_quality)
                # our dorado modbam files have TLEN as 0, why? Not all entries have an MN:i tag, why?
                if not x.is_secondary:
                    template_lengths.append(len(x.query_sequence))

                # if len(x.query_sequence) > 100000:
                #     print("caught a big one: {}".format(len(x.query_sequence)))
                #     print(x)

        d_phred[label] = numpy.array(phred_scores)
        d_mapq[label] = mapq_scores
        d_tlen[label] = template_lengths

        samfile.close()

    print_tlen_distribution(d_tlen)

    plot_qscore_hists(d_phred)
    plot_mapq_hists(d_mapq)

    try:
        plt.show()
    except:
        pass





