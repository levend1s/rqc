import pysam
import matplotlib.pyplot as plt
import argparse
import numpy
import pandas

parser = argparse.ArgumentParser(description="filter bam file by qscores / mapqs")
parser.add_argument('function')
parser.add_argument('inputs', nargs='+')
parser.add_argument('-q', '--min_phred', type=int, default=0)
parser.add_argument('-m', '--min_mapq', type=int, default=0)
parser.add_argument('-l', '--min_read_length', type=int, default=0)
parser.add_argument('-r', '--reverse_search', type=bool, default=False)
parser.add_argument('-n', '--num_results', type=int, default=0)
parser.add_argument('--sort_by', type=str, default="")
parser.add_argument('--check_duplicate_reads', type=bool, default=False)
parser.add_argument('-v', '--verbose', type=bool, default=False)
parser.add_argument('-o', '--outfile', type=str)

args = parser.parse_args()

INPUT = args.inputs
MIN_PHRED = args.min_phred
MIN_MAPQ = args.min_mapq
MIN_READ_LENGTH = args.min_read_length
FUNCTION = args.function
OUTFILE = args.outfile
CHECK_DUPLICATE_READS = args.check_duplicate_reads
VERBOSE = args.verbose
NUM_RESULTS = args.num_results
REVERSE_SEARCH = args.reverse_search
SORT_BY = args.sort_by

d_phred = {}
d_mapq = {}
d_tlen = {}
d_read_ids = {}

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

def calc_tlen_distribution(tlens):
    d = []

    for key in tlens.keys():
        if not tlens[key]:
            continue

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

        summary = {
            "sample": key,
            "count": len(cur_tlens),
            "total_bases_read": total_bases_read,
            "mean": mean,
            "median": median,
            "Q1": q1,
            "Q3": q3,
            "IQR": q3 - q1,
            "min": min_read,
            "max": max_read,
            "N10": n10,
            "N50": n50,
            "N90": n90
        }

        d.append(summary)

    df = pandas.DataFrame(d)
    return(df)


def plot_qscore_hists(qscores):
    for key in qscores.keys():
        bins = numpy.arange(numpy.floor(qscores[key].min()), numpy.ceil(qscores[key].max()))
        qscore_hist, bin_edges = numpy.histogram(qscores[key], bins=bins, density=True)
        # bin_edges has len(qscore_hist) + 1. We could find the midpoint of each edge, but let's be lazy and take the leftmost edge
        plt.plot(bin_edges[:-1], qscore_hist, label=key)  # arguments are passed to np.histogram

    plt.title("PHRED")
    plt.legend(loc="upper right")
    if (OUTFILE):
        plt.savefig("phred_{}".format(OUTFILE))

def plot_mapq_hists(mapqs):
    plt.figure()
    plt.title("MAPQ")
    plt.hist(mapqs.values(), histtype="bar", label=mapqs.keys(), density=True)
    plt.legend(loc="upper right")

    if (OUTFILE):
        plt.savefig("mapq_{}".format(OUTFILE))

def plot_read_distribution(tlens, read_summaries):
    plt.figure()
    plt.title("read lengths")

    print(read_summaries)

    # filter outliers from read_summaries
    outliers = {}
    filtered_for_outliers = {}
    for _, row in read_summaries.iterrows():
        label = row["sample"]
        upper = row["Q3"] + (row["IQR"] * 1.5)
        lower = row["Q1"] - (row["IQR"] * 1.5)

        outliers[label] = [x for x in tlens[label] if x > upper or x < lower]
        filtered_for_outliers[label] = [x for x in tlens[label] if x < upper and x > lower]

    labels = []
    for k in tlens.keys():
        labels.append(
            "{}" "\n"
            "n_reads={}" "\n"
            "n_outliers={}".format(k, len(filtered_for_outliers[k]), len(outliers[k]))
        )
    axes = plt.gca()
    axes.set_xticks(numpy.arange(1, len(labels) + 1), labels=labels)
    axes.set_ylabel("read length (nt)")

    plt.violinplot(filtered_for_outliers.values())

    if (OUTFILE):
        plt.savefig("readlengths_{}".format(OUTFILE))

    # ---------------- violin plot combining all samples ----------------
    plt.figure()
    plt.title("read lengths all samples")

    lengths_combined = {}
    lengths_combined_outliers = {}


    for k in filtered_for_outliers.keys():
        group = k.split("_")[1]
        if group not in lengths_combined:
            lengths_combined[group] = filtered_for_outliers[k]
            lengths_combined_outliers[group] = outliers[k]

        else:
            lengths_combined[group] += filtered_for_outliers[k]
            lengths_combined_outliers[group] += outliers[k]

    labels = []
    for k in lengths_combined.keys():
        labels.append(
            "{}" "\n"
            "n_reads={}" "\n"
            "n_outliers={}".format(k, len(lengths_combined[k]), len(lengths_combined_outliers[k]))
        )
    axes = plt.gca()
    axes.set_xticks(numpy.arange(1, len(labels) + 1), labels=labels)
    axes.set_ylabel("read length (nt)")

    plt.violinplot(lengths_combined.values())

    if (OUTFILE):
        plt.savefig("allreadlengths_{}".format(OUTFILE))



    

def find_multi_reference_alignments(d_reads):
    keys = list(d_reads.keys())
    d_intersects = {}

    if VERBOSE:
        print("starting multi ref check:")

    for i in list(range(len(keys))):
        for j in list(range(i+1, len(keys))):
            print("looking for read intersect in {} ({}) and {} ({})".format( 
                keys[i], len(d_reads[keys[i]]), keys[j], len(d_reads[keys[j]])
                ))

            intersect = [x for x in d_reads[keys[i]] if x in d_reads[keys[j]]]

            if VERBOSE:
                print("found {} common reads".format(len(intersect)))


            for r in intersect:
                if r in d_intersects:
                    if keys[j] not in d_intersects[r]:
                        d_intersects[r].append(keys[j])
                else:
                    d_intersects[r] = [keys[i], keys[j]]

    # print summaries
    print("FOUND {} READS COMMON TO SAMPLES:".format(len(d_intersects.keys())))

    for r in d_intersects.keys():
        print("{}\t".format(r), end="")

        for s in d_intersects[r]:
            print("{}\t".format(s), end="")
        print()

    # extra newline :_)
    print()

dataframes = {}

def process_bamfiles():
    for i in range(0, len(INPUT), 2):
        label = INPUT[i]
        filename = INPUT[i+1]

        phred_scores = []
        mapq_scores = []
        template_lengths = []
        read_ids = []


        samfile = pysam.AlignmentFile(filename, 'rb')
        iter = samfile.fetch()

        for x in iter:
            if (not x.is_secondary and 
                x.mapping_quality >= MIN_MAPQ and 
                x.get_tag('qs') >= MIN_PHRED and 
                len(x.query_sequence) >= MIN_READ_LENGTH):
                
                phred_scores.append(round(x.get_tag('qs'), 2))
                mapq_scores.append(x.mapping_quality)
                # our dorado modbam files have TLEN as 0, why? Not all entries have an MN:i tag, why?
                template_lengths.append(len(x.query_sequence))
                read_ids.append(x.query_name)

        d_phred[label] = numpy.array(phred_scores)
        d_mapq[label] = mapq_scores
        d_tlen[label] = template_lengths
        d_read_ids[label] = read_ids

        dataframes[label] = pandas.DataFrame({
            "read_id": read_ids,
            "phred_scores": phred_scores,
            "mapq_scores": mapq_scores,
            "template_lengths": template_lengths
        })

        samfile.close()

if FUNCTION == "search":
    print("searching...")
    process_bamfiles()

    for sample in dataframes.keys():
        print(sample)

        if (SORT_BY):
            dataframes[sample].sort_values(SORT_BY, inplace=True, ascending=(not REVERSE_SEARCH))

        if (NUM_RESULTS):
            print(dataframes[sample].head(NUM_RESULTS))
        else:
            print(dataframes[sample])

    if CHECK_DUPLICATE_READS:
        find_multi_reference_alignments(d_read_ids)


    df = calc_tlen_distribution(d_tlen)
    print(df)


if FUNCTION == "inspect":
    print("inspecting...")

    read_id = INPUT[0]
    alignments = []


    for i in range(1, len(INPUT), 2):
        label = INPUT[i]
        filename = INPUT[i+1]



        samfile = pysam.AlignmentFile(filename, 'rb')
        iter = samfile.fetch()

        for x in iter:
            if (x.query_name == read_id):
                alignments.append(x)


        samfile.close()

    for a in alignments:
        print(a)
    


if FUNCTION == "plot":
    print("plotting...")
    process_bamfiles()

    plot_qscore_hists(d_phred)
    plot_mapq_hists(d_mapq)

    read_summaries = calc_tlen_distribution(d_tlen)
    plot_read_distribution(d_tlen, read_summaries)


    try:
        plt.show()
    except:
        pass





