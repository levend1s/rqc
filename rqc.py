import pysam
import matplotlib.pyplot as plt
import argparse
import numpy
import pandas
import gffpandas.gffpandas as gffpandas
import sys
import scipy
import ast

numpy.seterr(divide='ignore', invalid='ignore')

parser = argparse.ArgumentParser(description="filter bam file by qscores / mapqs")
parser.add_argument('COMMAND')
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
parser.add_argument('--coverage_type', type=str, default = "gene")
parser.add_argument('--coverage_padding', type=int, default = 0)
parser.add_argument('--coverage_bins', type=int, default = 100)
parser.add_argument('--coverage_method', type=str, default = "max")

args = parser.parse_args()

INPUT = args.inputs
MIN_PHRED = args.min_phred
MIN_MAPQ = args.min_mapq
MIN_READ_LENGTH = args.min_read_length
COMMAND = args.COMMAND
OUTFILE = args.outfile
CHECK_DUPLICATE_READS = args.check_duplicate_reads
VERBOSE = args.verbose
NUM_RESULTS = args.num_results
REVERSE_SEARCH = args.reverse_search
SORT_BY = args.sort_by
FEATURE = args.feature
ANNOTATION_FILE= args.annotation_file
COVERAGE_TYPE= args.coverage_type
COVERAGE_PADDING=args.coverage_padding
COVERAGE_BINS=args.coverage_bins
COVERAGE_METHOD=args.coverage_method

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
        print("num bins {}".format(len(bins)))
        qscore_hist, bin_edges = numpy.histogram(qscores[key], bins=bins, density=True)
        # bin_edges has len(qscore_hist) + 1. We could find the midpoint of each edge, but let's be lazy and take the leftmost edge
        plt.plot(bin_edges[:-1], qscore_hist, label=key.replace("_pfal", ""))  # arguments are passed to np.histogram

    # plt.title("PHRED")
    plt.legend(loc="upper right")
    plt.xlabel("transcript quality (qscore)")
    plt.ylabel("density")

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
    # plt.title("read lengths")

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

    org_names = {
        '36C1_Pfal': 'P. falciparum',
        '36C1_Yeast': 'S. cerevisiae',
        '36C1_Human': 'H. sapien',

        'Pfal': 'P. falciparum',
        'Yeast': 'S. cerevisiae',
        'Human': 'H. sapien',

        'pfal': 'P. falciparum',
        'yeast': 'S. cerevisiae',
        'humans': 'H. sapien'
    }

    # labels = []
    # for k in tlens.keys():
    #     labels.append(
    #         "{}" "\n"
    #         "{:,} reads" "\n"
    #         "{:,} outliers".format(org_names[k], len(filtered_for_outliers[k]), len(outliers[k]))
    #     )
    # fig, axes = plt.subplots()
    # axes.set_xticks(numpy.arange(1, len(labels) + 1), labels=labels)
    # axes.set_ylabel("read length (nt)")
    # fig.tight_layout()
    # plt.violinplot(filtered_for_outliers.values())

    # if (OUTFILE):
    #     plt.savefig("readlengths_{}".format(OUTFILE))

    # ---------------- violin plot combining all samples ----------------
    plt.figure()
    # plt.title("read lengths all samples")

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
            "{:,} reads" "\n"
            "{:,} outliers".format(org_names[k], len(lengths_combined[k]), len(lengths_combined_outliers[k]))
        )
    axes = plt.gca()
    fig, axes = plt.subplots()
    axes.set_xticks(numpy.arange(1, len(labels) + 1), labels=labels)
    axes.set_ylabel("read length (nt)")
    fig.tight_layout()

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

if COMMAND == "search":
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


if COMMAND == "inspect":
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

if COMMAND == "base_coverage":
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

        output = pandas.DataFrame(columns = ["seq_id", "start", "end", "count_a", "count_c", "count_g", "count_t"])

        for index, row in matches.iterrows():
            # DEBUGGING
            print(index)
            # if i == num_matches:
            #     break

            # TODO: if two genes are close to each other, then this doesn't discern for only reads mapped to our gene of interest
            # so we can end up with weird lumps in the 5' end
            a, c, g, t = samfile.count_coverage(
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
            sum_c = sum(c)
            sum_g = sum(g)
            sum_t = sum(t)

            output.loc[index] = [row['seq_id']] + [row['start']] + [row['end']] + [sum_a] + [sum_c] + [sum_g] + [sum_t]


        print(output)
        print("total {} {} {} {}".format(
            output['count_a'].sum(),
            output['count_c'].sum(),
            output['count_g'].sum(),
            output['count_t'].sum()
        ))

# options: total, average, max
def resample_coverage(cov, bins, method):
    # handle the case where the subfeature (cov) is shorter than the number of bins
    if len(cov) < bins:
        temp_cov = []
        for e in cov:
            temp_cov += [e] * bins
        cov = temp_cov

    window_size = int(len(cov) / bins)

    resampled_coverage = numpy.zeros(bins)
    for window_index in range(bins):
        if method == "sum":
            resampled_coverage[window_index] = sum( cov[(window_size * window_index) : (window_size * (window_index+1))] )
        elif method == "max":
            resampled_coverage[window_index] = numpy.array(cov[(window_size * window_index) : (window_size * (window_index+1))]).max()
        elif method == "mean":
            resampled_coverage[window_index] = numpy.array(cov[(window_size * window_index) : (window_size * (window_index+1))]).mean()
        else:
            print("WARNING: resample coverage method not provided")

    return resampled_coverage

def normalise_coverage(cov):
    if cov.max() > 0:
        diff = cov.max() - cov.min()
        if diff == 0:
            normalised_resampled_coverage = cov * (1.0 / cov.max())
        else:
            normalised_resampled_coverage = (cov - cov.min()) * (1.0 / diff)

        return normalised_resampled_coverage
    else:
        return cov
    
# returns df of child features for tx_id
def get_feature_children(tx_id, gff_df):
    matches = gff_df.get_feature_by_attribute("Parent", tx_id)
    return matches

def plot_subfeature_coverage(coverages):
    fig, axes = plt.subplots()

    for label, cov in coverages['coverages'].items():
        plt.plot(cov, label=label)

    plt.legend(loc="upper right")
    axes.set_ylabel(coverages['y_label'], color="black")
    axes.set_ylim(ymin=0)
    plt.tick_params(
        axis='x',          
        which='both',
        bottom=False,
        top=False,
        labelbottom=False
    )

    textstr = "num matches: {}\nnum bins: {}".format(coverages['num_matches'], coverages['num_bins'])
    
    axes.text(
        0.05, 
        0.95, 
        textstr, 
        transform=axes.transAxes,
        fontsize=10,
        verticalalignment='top'
    )

    subfeature_names = coverages['subfeature_names']
    num_subfeatures = len(subfeature_names)
    subfeature_width = int(coverages['num_bins'] / num_subfeatures)

    # rename subfeature names to E1, E2, 3'UTR etc
    if num_subfeatures > 1:
        if 'UTR' in subfeature_names[0]:
            subfeature_names[0] = "5'UTR"
        if 'UTR' in subfeature_names[1]:
            subfeature_names[1] = "5'UTR"
        if 'UTR' in subfeature_names[-1]:
            subfeature_names[-1] = "3'UTR"
        if 'UTR' in subfeature_names[-2]:
            subfeature_names[-2] = "3'UTR"

    exon_idx = 1
    for i in range(num_subfeatures):
        if subfeature_names[i] == 'CDS':
            subfeature_names[i] = "E{}".format(exon_idx)
            exon_idx += 1

    label_rotation = 45
    if COVERAGE_TYPE == "gene":
        label_rotation = 0

    for l in range(num_subfeatures + 1):
        line_x_coord = subfeature_width * l
        label_x_coord = (line_x_coord + (subfeature_width / 2))
        plt.axvline(x= line_x_coord, color='darkgray', ls="--", linewidth=1.0)
        label_x_coord = line_x_coord + int(subfeature_width / 2)
        if l != num_subfeatures:
            axes.text(
                label_x_coord,
                -0.02, 
                subfeature_names[l], 
                fontsize=10,
                verticalalignment='top',
                horizontalalignment='center',
                rotation=label_rotation
            )

     # add vlines with intensity relative to value stored in sites_of_interest
    if 'sites_of_interest' in coverages and coverages['sites_of_interest'] is not None and coverages['sites_of_interest'].any():
        barcode_height = 0.05
        normalised_site_frequencies = normalise_coverage(coverages['sites_of_interest'])
        barcode_height_normalised_site_frequencies = normalised_site_frequencies * barcode_height

        for x in range(len(coverages['sites_of_interest'])):
            if normalised_site_frequencies[x]:
                plt.axvline(
                    x=x, 
                    ymin=1 - barcode_height, 
                    ymax=1 - barcode_height + barcode_height_normalised_site_frequencies[x], 
                    color='black', 
                    alpha=normalised_site_frequencies[x], 
                    ls="-",
                    linewidth=2.0
                )

    fig.tight_layout()

import time

gff_tree = {}
annotation_file = gffpandas.read_gff3(ANNOTATION_FILE)
gff_df = annotation_file.attributes_to_columns()

for row_index, row in gff_df.iterrows():
    if row['Parent'] in gff_tree:
        gff_tree[row['Parent']].append(row_index)
    else:
        gff_tree[row['Parent']] = [row_index]


def getSubfeatures(id, coverage_type, coverage_padding):

    if coverage_type == "subfeature":
        # plasmodium specific thing? drop exons, keep only CDS and UTR
        # exons seem to overlap with UTR regions in plasmodium gff

        # row_subfeatures = annotation_file.get_feature_by_attribute("Parent", [id])
        row_subfeatures = annotation_file.df.iloc[gff_tree[id]]
        row_subfeatures = row_subfeatures.sort_values(by=['start'])
        row_subfeatures = row_subfeatures[row_subfeatures.type != "exon"]

        # EDGE CASE: collapse multiple UTR's into a single UTR # eg utr_PF3D7_1105800.1_1
        if len(row_subfeatures[row_subfeatures.type == "five_prime_UTR"]) > 1:
            first_utr_index = row_subfeatures[row_subfeatures.type == "five_prime_UTR"].index[0]
            last_utr_index = row_subfeatures[row_subfeatures.type == "five_prime_UTR"].index[-1]
            row_subfeatures.at[first_utr_index, 'end'] = (row_subfeatures.loc[last_utr_index])['end']
            row_subfeatures.drop(index=last_utr_index, inplace=True)
        if len(row_subfeatures[row_subfeatures.type == "three_prime_UTR"]) > 1:
            first_utr_index = row_subfeatures[row_subfeatures.type == "three_prime_UTR"].index[0]
            last_utr_index = row_subfeatures[row_subfeatures.type == "three_prime_UTR"].index[-1]
            row_subfeatures.at[first_utr_index, 'end'] = (row_subfeatures.loc[last_utr_index])['end']
            row_subfeatures.drop(index=last_utr_index, inplace=True)

    elif coverage_type == "gene":
        # row_subfeatures = matches.iloc[index:index+1]
        row_subfeatures = gff_df[gff_df.ID == id]
        # drop everything after attributes
        row_subfeatures = row_subfeatures.drop(columns=row_subfeatures.columns[9:])
    else:
        print("WARNING: unknown coverage type: {}".format(coverage_type))

    # manually add upstream/downstream rows to the subfeature df
    # cannot assume 3' UTR occurs before the 5' UTR in the annotation file? But we just sorted by start so should be fine
    if coverage_padding > 0:
        first_feature = row_subfeatures.loc[row_subfeatures.index[0]]
        last_feature = row_subfeatures.loc[row_subfeatures.index[-1]]
        max_index = max(row_subfeatures.index)
        row_subfeatures.loc[max_index + 1] = [
            first_feature['seq_id'],
            "",
            "{}bp".format(coverage_padding),
            first_feature['start'] - coverage_padding,
            first_feature['start'],
            0,
            first_feature['strand'],
            "",
            ""
        ]
        row_subfeatures.loc[max_index + 2] = [
            last_feature['seq_id'],
            "",
            "{}bp".format(coverage_padding),
            last_feature['end'],
            last_feature['end'] + coverage_padding,
            0,
            last_feature['strand'],
            "",
            ""
        ]
        row_subfeatures = row_subfeatures.sort_values(by=['start'])

    return row_subfeatures


if COMMAND == "plot_coverage":
    # load annotation file
    feature_id = INPUT[0]

    # process input file. Each line contains a label, the type of file, and the filepath
    input_files = {}
    if len(INPUT[1:]) % 3 != 0:
        print("ERROR: not enough information in specified files. Check each input follows the format [LABEL] [TYPE] [PATH]")
        sys.exit()
    else:
        for in_index in range(1, len(INPUT), 3):
            if not INPUT[in_index].startswith("#"):
                input_files[INPUT[in_index]] = (INPUT[in_index+1], INPUT[in_index+2])

    annotation_file = gffpandas.read_gff3(ANNOTATION_FILE)

    # Try to find matches of provided type, if not, assume that input is a list of IDs
    matches = annotation_file.filter_feature_of_type([feature_id])
    if len(matches.df) == 0:
        evaluated_input = ast.literal_eval(feature_id)
        matches = annotation_file.get_feature_by_attribute("ID", evaluated_input)
        print("Looking for {} IDs, found {} matches. Plotting gene coverage for {}".format(len(evaluated_input) , len(matches.df), evaluated_input))
    else:
        print("Found {} matches for type {}. Plotting gene coverage...".format(len(matches.df), feature_id))

    num_matches = len(matches.df)
    PYSAM_PILEUP_MAX_DEPTH = 8000 # default
    subfeature_names = []
    matches = matches.attributes_to_columns()
    index = 0
    coverages = {}
    normalised_coverages = {}
    sites_of_interest = None

    for label, file in input_files.items():
        type, path = file
        feature_coverages = [None] * len(matches.index)
        normalised_feature_coverages = [None] * len(matches.index)

        if type == "bam":
            samfile = pysam.AlignmentFile(path, 'rb')
        elif type == "bedmethyl":
            modkit_bedmethyl_header = [
                "contig", "start", "end", "code", "score", "strand", 
                "start_2", "end_2", "color", "valid_cov", "percent_mod", "num_mod", 
                "num_canonical", "num_other_mod", "num_delete", "num_fail", "num_diff", "num_nocall"
            ]
            mods_file_df = pandas.read_csv(path, sep='\t', names=modkit_bedmethyl_header)
        elif type == "bed":
            modkit_bedmethyl_header = [
                "contig", "start", "end", "code", "score", "strand", 
                "start_2", "end_2", "color", "valid_cov", "percent_mod", "num_mod", 
                "num_canonical", "num_other_mod", "num_delete", "num_fail", "num_diff", "num_nocall"
            ]
            site_file_df = pandas.read_csv(path, sep='\t', names=modkit_bedmethyl_header)
        else:
            print("ERROR UNKNOWN FILE TYPE {}".format(type))

        # generate coverage for all matches in this bam file
        feature_index = 0
        for row_index, row in matches.iterrows():
            # find subfeatures
            row_subfeatures = getSubfeatures(row['ID'], COVERAGE_TYPE, COVERAGE_PADDING)

            if not subfeature_names:
                subfeature_names = row_subfeatures['type'].to_list()

            # gen coverage for each subfeature in a gene
            subfeature_index = 0
            num_subfeatures = len(row_subfeatures.index)
            subfeature_base_coverages = [None] * num_subfeatures

            for _, subfeature in row_subfeatures.iterrows():
                subfeature_length = subfeature['end'] - subfeature['start'] + 1
                subfeature_base_coverages[subfeature_index] = numpy.zeros(subfeature_length)

                if type == "bam":
                    # pysam indexes are zero indexed but gff are 1-indexed, so pysam index = gffindex-1
                    for column in samfile.pileup(
                        contig=subfeature['seq_id'], 
                        start=subfeature['start'] - 1, 
                        stop=subfeature['end'],
                        min_mapping_quality=MIN_MAPQ,
                        max_depth=PYSAM_PILEUP_MAX_DEPTH,
                        truncate = True
                    ):
                        # take the number of aligned reads at this column position (column.n) minus the number of aligned reads which have either a skip or delete base at this column position (r.query_position)
                        read_depth = len(list(filter(None, column.get_query_sequences())))
                        # reference pos is 0 indexed, gff (subfeature) is 1 indexed, add one to bring it back to zero

                        subfeature_base_coverages[subfeature_index][column.reference_pos - subfeature['start'] + 1] = read_depth

                elif type == "bedmethyl":
                    mod_matches = mods_file_df[
                        (mods_file_df.contig == subfeature['seq_id']) & 
                        (mods_file_df.start > subfeature['start']) & 
                        (mods_file_df.start < subfeature['end']) & 
                        (mods_file_df.strand == subfeature['strand'])
                    ]
                    # convert from genome to transcript space
                    subfeature_mod_positions = mod_matches['start'].to_numpy() - subfeature['start']
                    num_mods_at_pos = mod_matches['num_mod'].to_list()

                    for mod_pos_index in range(len(num_mods_at_pos)):
                        subfeature_base_coverages[subfeature_index][subfeature_mod_positions[mod_pos_index]] = num_mods_at_pos[mod_pos_index]

                elif type == "bed":
                    site_matches = site_file_df[
                        (site_file_df.contig == subfeature['seq_id']) & 
                        (site_file_df.start > subfeature['start']) & 
                        (site_file_df.start < subfeature['end']) & 
                        (site_file_df.strand == subfeature['strand'])
                    ]
                    # convert from genome to transcript space
                    subfeature_site_positions = site_matches['start'].to_numpy() - subfeature['start']

                    for site_pos_index in subfeature_site_positions:
                        subfeature_base_coverages[subfeature_index][site_pos_index] = 1
                else:
                    print("WARNING: unknown type: {}".format(type))

                subfeature_index += 1

            # resample coverage
            sf_bin_size = int(COVERAGE_BINS / num_subfeatures)
            sf_base_coverage_list = [None] * num_subfeatures

            # implicitly reset subfeature_index
            for subfeature_index in range(num_subfeatures):
                sf_resampled_coverage = resample_coverage(subfeature_base_coverages[subfeature_index], sf_bin_size, COVERAGE_METHOD)
                sf_base_coverage_list[subfeature_index] = sf_resampled_coverage

            # flatten resampled subfeature coverages into a single array
            resampled_base_coverage = numpy.concatenate(sf_base_coverage_list).ravel()
            # reverse coverages if necessary
            if (row["strand"] == "-"):
                resampled_base_coverage = numpy.flip(resampled_base_coverage)

            feature_coverages[feature_index] = resampled_base_coverage
            normalised_feature_coverages[feature_index] = normalise_coverage(resampled_base_coverage)

            feature_index += 1

            print("{}\t {}\t max coverage: {}\t tx length: {}".format(label, row['ID'], int(max(resampled_base_coverage)), row['end'] - row['start']))


        if type == "bam":
            samfile.close()

        # flatten down resampled coverage and store dict under label
        total_coverage = numpy.array([sum(i) for i in zip(*feature_coverages)])
        all_normalised_total_coverage = numpy.array([sum(i) for i in zip(*normalised_feature_coverages)])
        normalised_total_coverage = all_normalised_total_coverage * (1.0 / all_normalised_total_coverage.max())

        if type == "bed":
            sites_of_interest = normalised_total_coverage
        else:
            coverages[label] = total_coverage
            normalised_coverages[label] = normalised_total_coverage

    # plot coverages
    # this looks at coverage for each gene, resamples and normalises the coverage and adds it to a list
    # then takes the average of all those resampled and normalised coverages
    # this smooths out the cases where some genes might have read depth in the 1000's, and others in the 10's
    # so our data isn't skewed toward genes that are higher expressed
    # coverage: dict of read depths of transcript: eg {'sample1': [coverage...], 'sample2': [coverage...]}
    # mod_coverage: dict of mod coverages to plot: eg {'sample1m6As': [coverage...], 'sample2m6as': [coverage...], 'sample1pseU': [coverage...]}
    # sites_of_interest: dict of sites of interest to plot as vertical lines??? But how to do this for aggregate transcript searches?
    #       Maybe plot vlines for individual transcript plots and areas shaded with intensity according to how often motifs appear
    coverage_dict = {
        "coverages": coverages,
        "method": COVERAGE_METHOD,
        "num_matches": num_matches,
        "sites_of_interest": sites_of_interest,
        "num_bins": COVERAGE_BINS,
        "subfeature_names": subfeature_names,
        "y_label": "count (nt)"
    }

    plot_subfeature_coverage(coverage_dict)

    coverage_dict['coverages'] = normalised_coverages
    coverage_dict['y_label'] = "normalised coverage (au)"

    plot_subfeature_coverage(coverage_dict)

    plt.show()

    # if (OUTFILE):
    #     plt.savefig("coverage_{}".format(OUTFILE))




if COMMAND == "plot":
    print("plotting...")
    process_bamfiles()

    plot_qscore_hists(d_phred)
    #plot_mapq_hists(d_mapq)

    read_summaries = calc_tlen_distribution(d_tlen)
    plot_read_distribution(d_tlen, read_summaries)




    try:
        plt.show()
    except:
        pass



if COMMAND == "de":
    de = pandas.read_csv(INPUT[0], sep='\t')
    # print(de.head())

    de_filtered = de[de["adj.P.Val"] < 0.05]
    print(de_filtered)

if COMMAND == "plot_dmr":
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


if COMMAND == "find_dmr":
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
        
if COMMAND == "logo":
    import logomaker

    filename = INPUT[0]

    crp_df = -logomaker.get_example_matrix('crp_energy_matrix',
                                        print_description=False)
    
    f = open(filename, 'r')
    seqs = f.read().splitlines()
    counts_matrix = logomaker.alignment_to_matrix(seqs)
    print(counts_matrix)
    l = logomaker.Logo(counts_matrix)
    l.ax.set_xticklabels(["", "-2", "-1", "0", "+1", "+2"])
    l.ax.set_ylabel('count')
    plt.show()

if COMMAND == "plot_entropy":
    # read in an entropy regions bed file and plot a pdf histogram for mean entropy 
    column = INPUT[0]
    filename = INPUT[1]

    bed = pandas.read_csv(filename, sep='\t')
    entropy_values = bed[bed.columns[int(column)]].to_numpy()
    print(len(entropy_values))

    for i in range(2, len(INPUT)):
        filename = INPUT[i]
        b = pandas.read_csv(filename, sep='\t')
        e = b[b.columns[int(column)]].to_numpy()
        entropy_values = numpy.append(entropy_values, e)
        print(len(entropy_values))


    bins = numpy.arange(round(entropy_values.min(), 2), round(entropy_values.max(), 2), step=0.01)
    entropy_hist, bin_edges = numpy.histogram(entropy_values, bins=bins)
    # bin_edges has len(qscore_hist) + 1. We could find the midpoint of each edge, but let's be lazy and take the leftmost edge

    # print(qscore_hist)
    # print(bin_edges[:-1])

    log_entropy_hist = numpy.log10(entropy_hist)
    log_entropy_hist[log_entropy_hist == -numpy.inf] = 0

    xnew = numpy.linspace(bin_edges[:-1].min(), bin_edges[:-1].max(), 300) 
    spl = scipy.interpolate.make_interp_spline(bin_edges[:-1], log_entropy_hist, k=3)  # type: BSpline
    power_smooth = spl(xnew)

    axes = plt.gca()
    axes.set_ylim(ymin=0, ymax=4)
    plt.plot(xnew, power_smooth)

    plt.figure()
    plt.yscale('log')
    plt.plot(bin_edges[:-1], entropy_hist)  # arguments are passed to np.histogram


    plt.show()

    


