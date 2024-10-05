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
parser.add_argument('-f', '--feature', type=str)
parser.add_argument('-a', '--annotation_file', type=str)



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
FEATURE = args.feature
ANNOTATION_FILE= args.annotation_file


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
            "template_lengths": template_lengths,
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
            # print(len(dataframes[sample].cigar_tuples))
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
    
import gffpandas.gffpandas as gffpandas
import sys
import scipy

if FUNCTION == "base_coverage":
    # load annotation file
    feature_id = INPUT[0]
    print("summarising base coverage for {}...".format(feature_id))

    annotation_file = gffpandas.read_gff3(ANNOTATION_FILE)

    if feature_id == "chromosome":
        matches = pandas.DataFrame(columns = ["seq_id", "start", "end"])
        i = 0
        for line in annotation_file.header.splitlines():
            if "sequence-region" in line:
                s = line.split(" ")

                matches.loc[i] = [s[1]] + [int(s[2])] + [int(s[3])]
                i += 1
        num_matches = len(matches)


    else:
        matches = annotation_file.filter_feature_of_type([feature_id])
        if len(matches.df) == 0:
            print("WARNING: no matches of type {}".format(feature_id))
            matches = annotation_file.get_feature_by_attribute("ID", [feature_id])

            if len(matches.df) > 1:
                print("ERROR: multiple entries for {} found in gff file. exiting...".format(feature_id))
                sys.exit()

        num_matches = len(matches.df)
        matches = matches.df

    print("FOUND {} MATCHES FOR {}".format(num_matches, feature_id))


    for i in range(1, len(INPUT), 2):
        label = INPUT[i]
        filename = INPUT[i+1]
        samfile = pysam.AlignmentFile(filename, 'rb')

        output = pandas.DataFrame(columns = ["seq_id", "start", "end", "count_a", "count_t", "count_c", "count_g"])

        for index, row in matches.iterrows():
            # DEBUGGING
            print(index)
            # if i == num_matches:
            #     break

            # TODO: if two genes are close to each other, then this doesn't discern for only reads mapped to our gene of interest
            # so we can end up with weird lumps in the 5' end
            a, c, t, g = samfile.count_coverage(
                contig=row['seq_id'], 
                start=row['start'], 
                stop=row['end'],
                quality_threshold=0
            )
            # for column in samfile.pileup(
            #     contig=row['seq_id'], 
            #     start=row['start'], 
            #     stop=row['end'], 
            #     min_mapping_quality=MIN_MAPQ,
            #     truncate = True
            # ):
            #     print(column)

            sum_a = sum(a)
            sum_t = sum(t)
            sum_c = sum(c)
            sum_g = sum(g)

            output.loc[index] = [row['seq_id']] + [row['start']] + [row['end']] + [sum_a] + [sum_t] + [sum_c] + [sum_g]


        print(output)
        print("total {} {} {} {}".format(
            output['count_a'].sum(),
            output['count_t'].sum(),
            output['count_c'].sum(),
            output['count_g'].sum()
        ))
        
            


if FUNCTION == "plot_coverage":
    # load annotation file
    feature_id = INPUT[0]
    print("plotting gene coverage for {}...".format(feature_id))

    annotation_file = gffpandas.read_gff3(ANNOTATION_FILE)

    matches = annotation_file.filter_feature_of_type([feature_id])
    if len(matches.df) == 0:
        print("WARNING: no matches of type {}".format(feature_id))
        matches = annotation_file.get_feature_by_attribute("ID", [feature_id])

        if len(matches.df) > 1:
            print("ERROR: multiple entries for {} found in gff file. exiting...".format(feature_id))
            sys.exit()

    print("FOUND {} MATCHES FOR {}".format(len(matches.df), feature_id))

    num_matches = len(matches.df)
    coverage_lists = [None] * num_matches
    normalised_coverage_lists = [None] * num_matches

    num_bins = 100

    for i in range(1, len(INPUT), 2):
        label = INPUT[i]
        filename = INPUT[i+1]
        samfile = pysam.AlignmentFile(filename, 'rb')

        i = 0
        for index, row in matches.df.iterrows():
            # DEBUGGING
            print(i)
            # if i == num_matches:
            #     break

            # TODO: if two genes are close to each other, then this doesn't discern for only reads mapped to our gene of interest
            # so we can end up with weird lumps in the 5' end
            # a, c, t, g = samfile.count_coverage(
            #     contig=row['seq_id'], 
            #     start=row['start'], 
            #     stop=row['end'], 
            #     quality_threshold=MIN_PHRED
            # )
            # c = numpy.add(numpy.add(a, c), numpy.add(g, t))

            # This is the correct way to get gene coverage, since we can filter for reads of low mapping quality
            # It does not check the read is on the correct strand, and it does not check that the read mapped primarily
            # to the gene of interest (in the case two genes are close together, and some transcript overhand from one
            # gene bleeds into another gene coverage)
            c = [0] * (row['end'] - row['start'])
            j = 0
            for column in samfile.pileup(
                contig=row['seq_id'], 
                start=row['start'], 
                stop=row['end'], 
                min_mapping_quality=MIN_MAPQ,
                truncate = True
            ):
                # take the number of aligned reads at this column position (column.n) minus the number of aligned reads which have either a skip or delete base at this column position (r.query_position)
                read_depth = len(list(filter(None, column.get_query_sequences())))
                c[j] = read_depth
                j += 1
                

            if (row["strand"] == "-"):
                c = list(reversed(c))

            # resample coverage into an array of size num_bins
            x = numpy.arange(0, len(c), len(c) / num_bins)
            resampled_coverage = numpy.interp(x, range(len(c)), c)
            coverage_lists[i] = resampled_coverage

            if resampled_coverage.max() > 0:
                normalised_resampled_coverage = resampled_coverage * (1.0 / resampled_coverage.max())
                normalised_coverage_lists[i] = normalised_resampled_coverage
            else:
                normalised_coverage_lists[i] = resampled_coverage

            i += 1

        samfile.close()

    total_coverage = numpy.array([sum(i) for i in zip(*coverage_lists)])

    all_normalised_total_coverage = numpy.array([sum(i) for i in zip(*normalised_coverage_lists)])
    normalised_total_coverage = all_normalised_total_coverage * (1.0 / all_normalised_total_coverage.max())

    # this looks at coverage for each gene, resamples and normalises the coverage and adds it to a list
    # then takes the average of all those resampled and normalised coverages
    # this smooths out the cases where some genes might have read depth in the 1000's, and others in the 10's
    # so our data isn't skewed toward genes that are higher expressed 
    fig, axes = plt.subplots()
    plt.plot(total_coverage, color="red")
    axes.set_ylabel("total read depth (nt)", color="red")
    axes.set_ylim(ymin=0)

    # this looks at the coverage for each gene and resamples it, then takes the sum of all those resampled coverages and plots it
    # this can be skewed towards genes which have greater read depth
    axes_2 = axes.twinx()
    axes_2.set_ylabel("normalised read depth", color="blue")
    axes_2.plot(normalised_total_coverage, color="blue")
    axes_2.set_ylim(ymin=0)
    # axes_2.set_xlim(xmin=0, xmax=100)

    plt.title("read depth for {}".format(feature_id))
    fig.tight_layout()
    plt.show()

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



if FUNCTION == "de":
    de = pandas.read_csv(INPUT[0], sep='\t')
    # print(de.head())

    de_filtered = de[de["adj.P.Val"] < 0.05]
    print(de_filtered)

import seaborn

if FUNCTION == "plot_dmr":
    # plot transcript against dmr score

    d = {}

    heatmap_data = []

    for i in range(0, len(INPUT), 2):
        label = INPUT[i]
        filename = INPUT[i+1]

        bed = pandas.read_csv(filename, sep='\t')
        dmr_scores = bed.iloc[:,4]
        dmr_scores = dmr_scores.values.tolist()

        d[label] = dmr_scores
        heatmap_data.append(dmr_scores)

    dmr_cutoff = 1000
    high_dmr_scores = [x for x in dmr_scores if x > dmr_cutoff]

    # plt.violinplot(high_dmr_scores)
    # seaborn.stripplot(high_dmr_scores)
    a = numpy.random.random((16, 16))
    # plt.imshow([high_dmr_scores, high_dmr_scores], cmap='hot', interpolation='nearest')
    plt.plot(high_dmr_scores)
    plt.show()


if FUNCTION == "find_dmr":
    # we'll consider a region differentially methylated as a result of METTL3 knock sideways if
    # - the average scores across C1,2 vs K1,2 is above 1 likelihood
    # AND the average score across c1vsc2 and k1vsk2 is below 1 likelihood

    label = INPUT[0]
    filename = INPUT[1]
    bed = pandas.read_csv(INPUT[1], sep='\t')
    df = bed[bed.columns[3:5]]
    df.columns = ["region", label]

    for i in range(2, len(INPUT), 2):
        label = INPUT[i]
        filename = INPUT[i+1]

        this_bed = pandas.read_csv(filename, sep='\t')
        dmr_scores = this_bed[this_bed.columns[3:5]]
        dmr_scores.columns = ["region", label]

        # merge into df
        df = df.merge(dmr_scores, on="region", how="left")

    same_treatment_columns = df[df.columns[1:3]]
    diff_treatment_columns = df[df.columns[3:7]]

    df['same_treatment_average'] = same_treatment_columns.mean(numeric_only=True, axis=1)
    df['diff_treatment_average'] = diff_treatment_columns.mean(numeric_only=True, axis=1)

    SAME_TREATMENT_THRESHOLD = 0.1
    DIFF_TREATMENT_THRESHOLD = 10

    df['differentially_expressed'] = False
    df.loc[(df['same_treatment_average'] <= SAME_TREATMENT_THRESHOLD) & (df['diff_treatment_average'] >= DIFF_TREATMENT_THRESHOLD), 'differentially_expressed'] = True

    # write dataframe to file
    df.to_csv(OUTFILE, sep='\t')
        







