import os
import math
import numpy
import pandas
import matplotlib.pyplot as plt
import pysam
import gffpandas.gffpandas as gffpandas
import ast
import time

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


# exp function 
def exp_func(x, a, b, c):
    return a ** (x - b) + c

def power_func(x, a, b, c):
    return a * (x ** b) + c

def lin_func(x, a, b):
    return (a * x) + b

def normalise_numpy_array(a):
        return (a - numpy.min(a)) / (numpy.max(a) - numpy.min(a))

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

def process_bamfiles():
    for i in range(0, len(args.inputs), 2):
        label = args.inputs[i]
        filename = args.inputs[i+1]

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

# options: sum, average, max
def resample_coverage(cov, bins, method):

    window_size = len(cov)

    # increase size of cov so we can evenly bin it
    temp_cov = []
    for e in cov:
        temp_cov += [e] * bins
    cov = numpy.array(temp_cov)

    resampled_coverage = numpy.zeros(bins)
    for window_index in range(bins):

        if method == "sum":
            resampled_coverage[window_index] = math.ceil(sum( cov[(window_size * window_index) : (window_size * (window_index+1))] / bins ))

        elif method == "max":
            resampled_coverage[window_index] = numpy.array(cov[(window_size * window_index) : (window_size * (window_index+1))]).max()

        elif method == "mean":
            resampled_coverage[window_index] = numpy.array(cov[(window_size * window_index) : (window_size * (window_index+1))]).mean()
        else:
            print("WARNING: resample coverage method not provided")

    return resampled_coverage

def normalise_coverage(cov, min=0):
    if len(cov) > 0 and cov.max() > 0:
        diff = cov.max() - min
        if diff == 0:
            normalised_resampled_coverage = cov * (1.0 / cov.max())
        else:
            normalised_resampled_coverage = (cov - min) * (1.0 / diff)

        return normalised_resampled_coverage
    else:
        return cov
    
# returns df of child features for tx_id
def get_feature_children(tx_id, gff_df):
    matches = gff_df.get_feature_by_attribute("Parent", tx_id)
    return matches

def get_plot_color(l):
    this_color = "black"

    if 'read_depth' in l:
        this_color = 'skyblue'
    elif 'm6A' in l:
        this_color = 'indigo'
    elif 'm5C' in l:
        this_color = 'red'
    elif 'pseudouridine' in l:
        this_color = 'gold'
    elif 'inosine' in l:
        this_color = 'green'
    elif 'NGG' in l:
        this_color = 'green'
    elif 'TTTN' in l:
        this_color = 'blue'

    return this_color

def plot_subfeature_coverage(coverages, line_width, separate_y_axes, coverage_type):
    sample_labels = {}
    ymax = 0
    for label, cov in coverages['coverages'].items():
        sample_name = label.split("_")[0]

        if sample_name in sample_labels:
            sample_labels[sample_name].append(label)
        else:
            sample_labels[sample_name] = [label]
        
        if cov.max() > ymax:
            ymax = cov.max()


    num_samples = len(sample_labels.keys())
    height_ratios = [1 for x in range(num_samples)]
    if 'sites_of_interest' in coverages and coverages['sites_of_interest'] is not None and coverages['sites_of_interest'].any():
        barcode_ratio = 6
        height_ratios = [x*barcode_ratio for x in height_ratios]
        height_ratios.insert(0, 1)
        num_samples += 1

    fig, axes = plt.subplots(num_samples, gridspec_kw={'height_ratios': height_ratios})
    x_ticks = numpy.arange(coverages['num_bins'])

    plt_index = 0

    # add vlines with intensity relative to value stored in sites_of_interest
    if 'sites_of_interest' in coverages and coverages['sites_of_interest'] is not None and coverages['sites_of_interest'].any():
        this_axes = axes[plt_index]
        this_axes.bar(x_ticks, coverages['sites_of_interest'], color='black')
        this_axes.set_ylim(ymin=0, ymax=coverages['sites_of_interest'].max())
        this_axes.set_xlim(xmin=0, xmax=coverages['num_bins']-1)
        this_axes.set_yticks([0, coverages['sites_of_interest'].max()])
        this_axes.set_xticks([])
        this_axes.set_ylabel("DRACH\nfrequency", color="black")

        plt_index += 1

    for k, v in sample_labels.items():
        if num_samples > 1:
            this_axes = axes[plt_index]
        else:
            this_axes = axes

        if (len(v) == 2) and separate_y_axes:
            cov_1_label = v[0]
            cov_1 = coverages['coverages'][cov_1_label]
            cov_1_color = get_plot_color(cov_1_label)

            cov_2_label = v[1]
            cov_2 = coverages['coverages'][cov_2_label]
            cov_2_color = get_plot_color(cov_2_label)

            this_axes.plot(cov_1, label= ' '.join(cov_1_label.split("_")[1:]), color=cov_1_color, linewidth=line_width)
            this_axes.fill_between(x_ticks, cov_1, alpha=0.2, color=cov_1_color)
            this_axes.tick_params(axis='y', labelcolor=cov_1_color)
            this_axes.set_ylim(ymin=0)
            this_axes.set_xlim(xmin=0, xmax=coverages['num_bins']-1)
            this_axes.set_ylabel(' '.join(cov_1_label.split("_")[1:]), color=cov_1_color)

            this_axes_2 = this_axes.twinx()
            this_axes_2.plot(cov_2, label= ' '.join(cov_2_label.split("_")[1:]), color=cov_2_color, linewidth=line_width)
            this_axes_2.fill_between(x_ticks, cov_2, alpha=0.2, color=cov_2_color)
            this_axes_2.tick_params(axis='y', labelcolor=cov_2_color)
            this_axes_2.set_ylim(ymin=0)
            this_axes_2.set_ylabel(' '.join(cov_2_label.split("_")[1:]), color=cov_2_color)

        else:
            for label in v:
                cov = coverages['coverages'][label]
                this_color = get_plot_color(label)

                this_axes.plot(cov, label= ' '.join(label.split("_")[1:]), color=this_color, linewidth=line_width)
                this_axes.fill_between(x_ticks, cov, alpha=0.2, color=this_color)

            this_axes.legend(loc="upper left", title=k)
            this_axes.set_ylabel(coverages['y_label'], color="black")
            this_axes.set_ylim(ymin=0, ymax=ymax*1.1)
            this_axes.set_xlim(xmin=0, xmax=coverages['num_bins']-1)
            this_axes.set_yticks([0, ymax])
            this_axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        this_axes.tick_params(
            axis='x',          
            which='both',
            bottom=False,
            top=False,
            labelbottom=False
        )

        num_subfeatures = len(coverages['subfeature_names'])

        label_rotation = 45
        if coverage_type in ("gene"):
            label_rotation = 0

        curr_pos = 0

        for l in range(num_subfeatures):
            subfeature_width = coverages['subfeature_info'][coverages['subfeature_names'][l]]
            line_x_coord = curr_pos + subfeature_width

            this_axes.axvline(x= line_x_coord, color='darkgray', ls="--", linewidth=line_width)
            label_x_coord = curr_pos + int(subfeature_width / 2)

            if plt_index == (num_samples - 1):
                this_axes.text(
                    label_x_coord,
                    ymax*-0.02,
                    coverages['subfeature_names'][l], 
                    fontsize=10,
                    verticalalignment='top',
                    horizontalalignment='center',
                    rotation=label_rotation
                )

            curr_pos += subfeature_width

        plt_index += 1

def getSubfeatures(annotation, id, coverage_type, coverage_padding):

    if coverage_type in ["subfeature", "subfeature_cds"]:
        # plasmodium specific thing? drop exons, keep only CDS and UTR
        # exons seem to overlap with UTR regions in plasmodium gff

        # row_subfeatures = ANNOTATION_FILE.df.iloc[GFF_PARENT_TREE[id]]
        row_subfeatures = annotation[annotation['Parent'] == id]
        row_subfeatures = row_subfeatures.sort_values(by=['start'])
        row_subfeatures = row_subfeatures[row_subfeatures.type != "exon"]        

        # EDGE CASE: collapse multiple UTR's into a single UTR # eg utr_PF3D7_1105800.1_1
        # if len(row_subfeatures[row_subfeatures.type == "five_prime_UTR"]) > 1:
        #     first_utr_index = row_subfeatures[row_subfeatures.type == "five_prime_UTR"].index[0]
        #     last_utr_index = row_subfeatures[row_subfeatures.type == "five_prime_UTR"].index[-1]
        #     row_subfeatures.at[first_utr_index, 'end'] = (row_subfeatures.loc[last_utr_index])['end']
        #     row_subfeatures.drop(index=last_utr_index, inplace=True)
        # if len(row_subfeatures[row_subfeatures.type == "three_prime_UTR"]) > 1:
        #     first_utr_index = row_subfeatures[row_subfeatures.type == "three_prime_UTR"].index[0]
        #     last_utr_index = row_subfeatures[row_subfeatures.type == "three_prime_UTR"].index[-1]
        #     row_subfeatures.at[first_utr_index, 'end'] = (row_subfeatures.loc[last_utr_index])['end']
        #     row_subfeatures.drop(index=last_utr_index, inplace=True)

    elif coverage_type == "gene":
        # row_subfeatures = matches.iloc[index:index+1]
        row_subfeatures = GFF_DF[GFF_DF.ID == id]
        # drop everything after attributes
        row_subfeatures = row_subfeatures.drop(columns=row_subfeatures.columns[9:])
    else:
        print("WARNING: unknown coverage type: {}".format(coverage_type))

    # manually add upstream/downstream rows to the subfeature df
    # cannot assume 3' UTR occurs before the 5' UTR in the annotation file? But we just sorted by start so should be fine
    if coverage_padding > 0:
        first_feature = row_subfeatures.loc[row_subfeatures.index[0]]
        last_feature = row_subfeatures.loc[row_subfeatures.index[-1]]
        max_index = len(row_subfeatures)-1

        row_subfeatures = pandas.concat([row_subfeatures, row_subfeatures.iloc[[0, -1]]], ignore_index=True)

        row_subfeatures.loc[max_index + 1, 'start'] = first_feature['start'] - coverage_padding
        row_subfeatures.loc[max_index + 1, 'end'] = first_feature['start']
        row_subfeatures.loc[max_index + 1, 'type'] = "{}bp".format(coverage_padding)

        row_subfeatures.loc[max_index + 2, 'start'] = last_feature['end']
        row_subfeatures.loc[max_index + 2, 'end'] = last_feature['end'] + coverage_padding
        row_subfeatures.loc[max_index + 2, 'type'] = "{}bp".format(coverage_padding)

        row_subfeatures = row_subfeatures.sort_values(by=['start'])

    return row_subfeatures

def START_CLOCK(name):
    CLOCKS[name] = time.time()

def STOP_CLOCK(name, stop_name):
    clock_stop = time.time()
    diff = clock_stop - CLOCKS[name]
    if DEBUG:
        print("\tDEBUG: time between {} and {}: {}s".format(name, stop_name, diff))

def process_input_files(INPUT):
    input_files = {}

    with open(INPUT, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("#") or line.strip() == "":
                continue
            else:
                line = line.split()
                input_files[line[0]] = {
                    'group': line[1], 
                    'type': line[2], 
                    'path': line[3]
                }

    return input_files

def process_annotation_file(annotation_file_path=None):
    print("LOG - loading gff file...")
    # load annotation file and find indexes for all parent children
    ANNOTATION_FILE = gffpandas.read_gff3(annotation_file_path)

    gff_df = ANNOTATION_FILE.attributes_to_columns()
    gff_df['strand'] = gff_df['strand'].astype('category')
    gff_df['seq_id'] = gff_df['seq_id'].astype('category')
    gff_df['ID'] = gff_df['ID'].astype('category')
    gff_df['type'] = gff_df['type'].astype('category')

    # for row_index, row in GFF_DF.iterrows():
    #     if row['Parent'] in GFF_PARENT_TREE:
    #         GFF_PARENT_TREE[row['Parent']].append(row_index)
    #     else:
    #         GFF_PARENT_TREE[row['Parent']] = [row_index]

    return gff_df

def process_genome_file(genome_file_path=None, contig_lengths = None):
    print("LOG - loading fasta file...")

    FASTA_DICT = {}
    FASTA_LINE_LENGTH = 60

    # read two lines and set the line length as the second line
    if os.path.isfile(genome_file_path):
        with open(genome_file_path) as f:
            line = f.readline()
            line = f.readline()
            FASTA_LINE_LENGTH = len(line.strip())
            print("LOG - fasta line len = {}...".format(FASTA_LINE_LENGTH))
            print(contig_lengths)



    if os.path.isfile(genome_file_path):
        with open(genome_file_path) as f:
            line = f.readline()
            current_contig = None

            while line:
                # is line a header
                # PLASMODIUM SPECIFIC
                if line.startswith('>'):
                    if not contig_lengths:
                        header_attrs = line.split('|')
                        for h in header_attrs:
                            h = h.strip()
                            if h.startswith('>'):
                                contig = h[1:]
                                FASTA_DICT[contig] = {}
                                current_contig = contig
                            else:
                                k, v = h.split('=')
                                FASTA_DICT[current_contig][k] = v
                        FASTA_DICT[current_contig]['length'] = int(FASTA_DICT[current_contig]['length'])
                    else:
                        header_attrs = line.split(' ')
                        current_contig = header_attrs[0][1:]
                        FASTA_DICT[current_contig] = {}
                        FASTA_DICT[current_contig]['length'] = contig_lengths[current_contig]
                    
                    num_lines_seq_fasta = math.ceil(FASTA_DICT[current_contig]['length'] / FASTA_LINE_LENGTH) 
                    FASTA_DICT[current_contig]['sequence'] = [None] * num_lines_seq_fasta
                    i = 0
                else:
                    if line.startswith('>'):
                        current_contig = None

                    elif current_contig:
                        # FASTA_DICT[current_contig]['sequence'] += line.strip()
                        FASTA_DICT[current_contig]['sequence'][i] = line.strip()
                        i += 1

                line = f.readline()
    
    # verify each chromosome is the correct length
    for k, v, in FASTA_DICT.items():
        FASTA_DICT[k]['sequence'] = "".join(FASTA_DICT[k]['sequence'])
        print("{} - loaded seq: {}, header length: {}".format(k, len(v['sequence']), v['length']))

    return FASTA_DICT

def filter_gff_for_target_features(annotation_file):
    print("LOG - filtering gff file for target features...")

    # Try to find matches of provided type, if not, assume that input is a list of IDs
    feature_id = args.inputs[0]
    if os.path.isfile(feature_id):
        lines = []
        with open(feature_id) as f:
            lines = f.read().splitlines()

        matches = annotation_file.get_feature_by_attribute("ID", lines)

    else:
        matches = annotation_file.filter_feature_of_type([feature_id])
        if len(matches.df) == 0:
            evaluated_input = ast.literal_eval(feature_id)
            matches = annotation_file.get_feature_by_attribute("ID", evaluated_input)
            print("LOG - Looking for {} IDs, found {} matches: {}".format(len(evaluated_input) , len(matches.df), evaluated_input))
        else:
            print("LOG - Found {} matches for type {}...".format(len(matches.df), feature_id))

    matches = matches.attributes_to_columns()

    return matches

def find_canonical_mods(gff_panda_rows, input_files, bam_labels, mod_prop_threshold, read_depth_threshold, padding):
    print("LOG - finding canonical mods...")
    cannonical_mods_start_pos = {}
    for label in bam_labels:
        prefix = label.split("_")[0]
        mod_label = "{}_m6A_0.95".format(prefix)

        if DEBUG:
            print("mod_prop_threshold: {}".format(mod_prop_threshold))
            print("read_depth_threshold: {}".format(read_depth_threshold))

        mods_file_df = pandas.read_csv(input_files[mod_label]['path'], sep='\t', names=MODKIT_BEDMETHYL_HEADER)
        mods_file_df['strand'] = mods_file_df['strand'].astype('category')
        mods_file_df['contig'] = mods_file_df['contig'].astype('category')
        
        # TODO maybe don't need to filter this for read depth, just filter the gene for read depth
        mods_file_df = mods_file_df[
            (mods_file_df.percent_mod >= (mod_prop_threshold * 100)) & 
            (mods_file_df.valid_cov >= read_depth_threshold)
        ]

        for _, row in gff_panda_rows.iterrows():

            if row['ID'] not in cannonical_mods_start_pos:
                cannonical_mods_start_pos[row['ID']] = []

            row_mods = mods_file_df[
                (mods_file_df.start >= (row['start'] - padding)) &
                (mods_file_df.end <= (row['end'] + padding)) &
                (mods_file_df.strand == row['strand']) &
                (mods_file_df.contig == row['seq_id'])
            ]

            if DEBUG:
                print(row_mods)

            for mod_index, mod in row_mods.iterrows():
                if mod['start'] not in cannonical_mods_start_pos[row['ID']]:
                    cannonical_mods_start_pos[row['ID']].append(mod['start'])

    if DEBUG:
        print("cannonical_mods_start_pos: {}".format(cannonical_mods_start_pos[row['ID']]))

    return cannonical_mods_start_pos

def approximate_tes(gff_panda_rows, input_files, bam_labels):
    d_approximated_tes = {}

    for _, row in gff_panda_rows.iterrows():
        d_tes_vs_prop = {}
        d_readthrough_split_points = {}
        d_fitted_curve_r_squared = {}
        d_fitted_curve_coeff = {}

        # calculate readthrough proportions for each sample
        print("LOG - {} - approximating transcript end site...".format(row.ID))

        for label in bam_labels:
            samfile = pysam.AlignmentFile(input_files[label]['path'], 'rb')

            reads_in_region = samfile.fetch(
                contig=row.seq_id,
                start=row.start, 
                stop=row.end
            )

            reads_in_region = list(reads_in_region)

            samfile.close()

            tes = []

            for i in range(len(reads_in_region)):
                r = reads_in_region[i]

                # keep only reads in the same direction as this strand
                if (row.strand == "+" and r.is_forward) or (row.strand == "-" and r.is_reverse):
                    # if (row.strand == "-" and r.reference_start <= row.start) or (row.strand == "+" and r.reference_end >= row.end):
                    if row.strand == "-":
                        tes.append(r.reference_start)
                    else:
                        tes.append(r.reference_end)

            d_tes_vs_prop[label] = []

            for p in tes:
                if row['strand'] == "-":
                    num_read_throughs = len([x for x in tes if x < p])
                    num_normal = len([x for x in tes if x >= p])
                else:
                    num_read_throughs = len([x for x in tes if x > p])
                    num_normal = len([x for x in tes if x <= p])
                
                rt_prop = num_read_throughs / num_normal

                # only keep this if there are fewer readthroughs than normals
                # NOTE this assumption doesn't work if there is one clean cut site for all reads
                # or the majority of reads
                if rt_prop < 1:
                    d_tes_vs_prop[label].append((p, rt_prop))

            # if we couldn't add any TES which had fewer readthroughs than normals
            # Then this must be a clean cut gene... The final TES must be the splitpoint
            if len(d_tes_vs_prop[label]) == 1:
                d_readthrough_split_points[label] = d_tes_vs_prop[label][0][0]
                d_fitted_curve_r_squared[label] = numpy.inf
                d_fitted_curve_coeff[label] = numpy.inf
                # TODO add type info to reporting dict
                continue

            # sort by genomic position
            sorted_tes_prop = sorted(d_tes_vs_prop[label], key=lambda a: a[0], reverse=False)
            pos = [x[0] for x in sorted_tes_prop]
            prop = [x[1] for x in sorted_tes_prop]

            # fit curve to TES data that will be used to to find the kneedle of the curve
            # data will have bumps and many false knees/elbows that we want to smooth out
            # so the kneedle function finds the correct knee

            # interpolate points between the few points we have so we can better approximate a
            # curve that fits our data
            pos_normalised = normalise_numpy_array(numpy.array(pos))
            x_interp = numpy.linspace(0, 1, 100)
            pos_y_interp = scipy.interpolate.interp1d(pos_normalised, prop)
            prop_normalised_interpolated = [pos_y_interp(x) for x in x_interp]

            # if strand is positive direction in kneedle function are reversed
            if row['strand'] == "-":
                elbow_direction = "increasing"
                initial_guess = [1, 0.1, 0]
            else:
                elbow_direction = "decreasing"
                initial_guess = [-1, 0.1, 1]

            # len(prop) < 100
            abc, _ = scipy.optimize.curve_fit(
                power_func,
                x_interp,
                prop_normalised_interpolated,
                p0=initial_guess,
                maxfev=5000
            )

            y_fitted = power_func(x_interp, *abc)
            if DEBUG:
                print(abc)
                plt.scatter(x_interp, prop_normalised_interpolated)
                plt.scatter(x_interp, y_fitted)
                plt.show()
            # calculate R^2
            # https://stackoverflow.com/questions/19189362/getting-the-r-squared-value-using-curve-fit
            residuals = prop_normalised_interpolated - power_func(x_interp, *abc)
            residual_sum_squares = numpy.sum(residuals ** 2)
            total_sum_squares = numpy.sum((prop_normalised_interpolated - numpy.mean(prop_normalised_interpolated)) ** 2)
            r_squared = 1 - (residual_sum_squares / total_sum_squares)

            if DEBUG:
                print("r_squared: {}".format(r_squared))

            fitted_curve_coeff = abc[1]
            d_fitted_curve_coeff[label] = fitted_curve_coeff

            # if we estimated a concave curve, take the final TES as end site
            if fitted_curve_coeff <= 1 and row['strand'] == "-":
                d_readthrough_split_points[label] = d_tes_vs_prop[label][0][0]
                d_fitted_curve_r_squared[label] = numpy.inf
                continue

            if fitted_curve_coeff >= 1 and row['strand'] == "+":
                d_readthrough_split_points[label] = d_tes_vs_prop[label][-1][0]
                d_fitted_curve_r_squared[label] = numpy.inf
                continue

            # as the genomic position increases, the splitpoint readthrough proportion gets bigger
            kneedle = KneeLocator(x_interp, y_fitted, S=1.0, curve='convex', direction=elbow_direction)

            # convert normalised elbow point back to genomic space
            if kneedle.knee:
                normalised_elbow_pos = kneedle.knee
            else:
                print("WARNING: {} - couldn't determine knee, setting as max TES".format(row['ID']))
                if row['strand'] == "-":
                    normalised_elbow_pos = 0
                else:
                    normalised_elbow_pos = 1

            genomic_elbow_pos = pos[0] + ((max(pos) - min(pos)) * normalised_elbow_pos)
            d_readthrough_split_points[label] = genomic_elbow_pos
            d_fitted_curve_r_squared[label] = r_squared

        if DEBUG:
            print("{} - readthrough_split_points: {}".format(row['ID'], d_readthrough_split_points))

        # take the average of the control readthrough splitpoints and r^2
        readthrough_split_point = 0
        # average_r_squared = 0
        for label in bam_labels_control:
            readthrough_split_point += d_readthrough_split_points[label]
            # average_r_squared += d_fitted_curve_r_squared[label]

        readthrough_split_point = int(readthrough_split_point / len(bam_labels_control))
        d_approximated_tes[row.ID] = readthrough_split_point


    return d_approximated_tes

def get_filtered_reads_ids(gff_panda_rows, input_files, bam_labels, mod_prop_threshold, read_depth_threshold, padding, cannonical_mods_start_pos, approximated_tes = None):
    print("LOG - filtering reads...")

    print("" \
    "label\t" \
    "id\t" \
    "filtered reads\t" \
    "reads in region\t" \
    "filtered (strand)\t" \
    "filtered (fc)\t" \
    "filtered (3p)\t")
    
    PYSAM_MOD_THRESHOLD = int(256 * mod_prop_threshold) 

    d_sample_filtered_read_ids = {}

    for label in bam_labels:
        samfile = pysam.AlignmentFile(input_files[label]['path'], 'rb')
        
        # attempt to find the relevent featureCounts file in input_files
        if FILTER_FOR_FEATURE_COUNTS:
            feature_counts_sample_label = label.split("_")[0] + "_featureCounts"
            feature_counts_df = pandas.read_csv(input_files[feature_counts_sample_label]['path'], sep='\t', names=FEATURECOUNTS_HEADER)
            feature_counts_df['targets'] = feature_counts_df['targets'].astype('category')

        d_sample_filtered_read_ids[label] = {}

        # generate coverage for all matches in this bam file
        for _, row in gff_panda_rows.iterrows():
            reads_in_region = samfile.fetch(
                contig=row['seq_id'], 
                start=row['start'] - padding, 
                stop=row['end'] + padding
            )
            reads_in_region = list(reads_in_region)

            samfile.close()

            d_filtered_read_ids = {}
            d_filtered_read_ids['different_strand'] = []
            d_filtered_read_ids['3p'] = []
            d_filtered_read_ids['featureCounts'] = []
            d_filtered_read_ids['filtered_reads'] = []
            d_filtered_read_ids['tx_end_sites'] = {}

            d_filtered_read_ids['mods'] = {}

            d_filtered_read_ids['mods']['None'] = []
            d_filtered_read_ids['tx_end_sites']['all'] = []
            d_filtered_read_ids['tx_end_sites']['None'] = []

            for mod in cannonical_mods_start_pos[row.ID]:
                d_filtered_read_ids['mods'][mod] = []
                d_filtered_read_ids['tx_end_sites'][mod] = []

            if FILTER_FOR_FEATURE_COUNTS:
                gene_reads = feature_counts_df[feature_counts_df.targets == row['ID'].split(".")[0]]
                gene_read_ids_fc = set(gene_reads['read_id'])

            if approximated_tes:
                gene_3p_end = approximated_tes[row.ID]
                print("LOG - using approximated tes {}".format(gene_3p_end))
            else:
                if row.strand == "-":
                    gene_3p_end = row.start
                    gene_5p_end = row.end
                else:
                    gene_3p_end = row.end
                    gene_5p_end = row.start

            # NOTE: now the longest function in the TES analysis
            for i in range(len(reads_in_region)):
                r = reads_in_region[i]

                # keep only reads in the same direction as this strand
                if (row.strand == "+" and r.is_forward) or (row.strand == "-" and r.is_reverse):
                    if (row.strand == "-" and r.reference_start <= gene_3p_end and r.reference_end >= gene_3p_end) \
                    or (row.strand == "+" and r.reference_end >= gene_3p_end and r.reference_start <= gene_3p_end) \
                    or (row.strand == "-" and r.reference_start >= gene_3p_end and r.reference_end <= gene_5p_end) \
                    or (row.strand == "+" and r.reference_end <= gene_3p_end and r.reference_start >= gene_5p_end):
                        # get mod sites for this read. this is [(read index, 256 * mod_prob)...]
                        if FILTER_FOR_FEATURE_COUNTS and (r.qname not in gene_read_ids_fc):
                            d_filtered_read_ids['featureCounts'].append(i)
                        else:
                            d_filtered_read_ids['filtered_reads'].append(i)

                            if row.strand == "-":
                                mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES['m6A_rev'])
                                tx_end_site = r.reference_start
                            else:
                                mods_probs = r.modified_bases.get(PYSAM_MOD_TUPLES['m6A_for'])
                                tx_end_site = r.reference_end

                            d_filtered_read_ids['tx_end_sites']['all'].append(tx_end_site)

                            if mods_probs:
                                # keep only mod positions which are above mod prob threshold
                                ref_pos = numpy.array(r.get_reference_positions(full_length=True))
                                read_mod_positions = [x[0] for x in mods_probs if x[1] >= PYSAM_MOD_THRESHOLD]
                                
                                # read mod positions is the position from the start of the read
                                # aligned reads may contain indels, so we need to get reference index from get_reference_positions
                                canonical_mods_in_read = set(cannonical_mods_start_pos[row.ID]).intersection(ref_pos[read_mod_positions])
                                if len(canonical_mods_in_read) == 0:
                                    d_filtered_read_ids['mods']['None'].append(i)
                                    d_filtered_read_ids['tx_end_sites']['None'].append(tx_end_site)
                                else:
                                    for mod in canonical_mods_in_read:
                                        d_filtered_read_ids['mods'][mod].append(i)
                                        d_filtered_read_ids['tx_end_sites'][mod].append(tx_end_site)
                    else:
                        d_filtered_read_ids['3p'].append(i)
                else:
                    d_filtered_read_ids['different_strand'].append(i)

            print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                label, 
                row.ID, 
                len(d_filtered_read_ids['filtered_reads']), 
                len(reads_in_region), 
                len(d_filtered_read_ids['different_strand']), 
                len(d_filtered_read_ids['featureCounts']), 
                len(d_filtered_read_ids['3p'])
                ))
            
            d_sample_filtered_read_ids[label][row.ID] = d_filtered_read_ids

    return d_sample_filtered_read_ids

BASE_COMPLEMENTS = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    '[': ']',
    ']': '[',
    '^': '$',
    '$': '^',
    '(': ')',
    ')': '('
}

def reverse_complement(s):
    MOTIF_REVERSED = s[::-1]

    MOTIF_REVERSE_COMPLEMENT = ""
    for l in MOTIF_REVERSED:
        if l in BASE_COMPLEMENTS.keys():
            MOTIF_REVERSE_COMPLEMENT += BASE_COMPLEMENTS[l]
        else:
            MOTIF_REVERSE_COMPLEMENT += l

    return MOTIF_REVERSE_COMPLEMENT

