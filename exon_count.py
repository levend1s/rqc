import gffpandas.gffpandas as gffpandas
from pprint import pprint
from matplotlib import pyplot as plt
import json
import numpy
import pandas
import math

ANNOTATION_FILE = "/Users/joshualevendis/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"
OUTPUT_FILE = "/Users/joshualevendis/Documents/RNA/rqc/pfal_mRNA_exon_counts.json"

FEATURE_ID = "CDS"
PARENT_TYPE = "mRNA"

annotation_file = gffpandas.read_gff3(ANNOTATION_FILE)
gff_df = annotation_file.attributes_to_columns()

exons = annotation_file.filter_feature_of_type([FEATURE_ID])
three_prime_utrs = annotation_file.filter_feature_of_type(["three_prime_UTR"])

dict_exon_count = {}
feature_types = {}
dict_exon_count_counts = {}
dict_utr_count_counts = {}

gff_tree = {}

# return dict of attributes given separators
def split_attributes(attr, attr_sep, val_sep):
    d = {}
    for a in attr.split(attr_sep):
        k, v = a.split(val_sep)
        d[k] = v

    return d

# build dict of feature and type
for row_index, row in annotation_file.df.iterrows():
    attrs = split_attributes(row['attributes'], ';', '=')
    feature_types[attrs['ID']] = row.type

for row_index, row in gff_df.iterrows():
    if row['Parent'] in gff_tree:
        gff_tree[row['Parent']].append(row['ID'])
    else:
        gff_tree[row['Parent']] = [row['ID']]

# build dict of exon count
for row_index, row in exons.df.iterrows():
    attrs = split_attributes(row['attributes'], ';', '=')

    # if parent is mRNA add it to the dicts
    if feature_types[attrs['Parent']] == PARENT_TYPE:
        if attrs['Parent'] in dict_exon_count:
            dict_exon_count[attrs['Parent']] += 1
        else:
            dict_exon_count[attrs['Parent']] = 1

for k, v in dict_exon_count.items():
    if v not in dict_exon_count_counts:
        dict_exon_count_counts[v] = 1
    else:
        dict_exon_count_counts[v] += 1

def get_transcripts_with_n_exons(d, n):
    l = []
    for k, v in d.items():
        if v == n:
            l.append(k)

    return l

def split_transcripts_by_utrs(d):
    utrs = []
    no_utrs = []
    for k, v in d.items():
        break

    return utrs, no_utrs



# TODO: plot that shows isoforms
def get_transcripts_with_isoforms():
    return

def build_outdict():
    # exon_count_keys = sorted(dict_exon_count_counts.keys())
    # out_dict = {}
    # for c in exon_count_keys:
    #     print(c)
    #     out_dict[c] = get_transcripts_with_n_exons(dict_exon_count, c)
    #     print(get_transcripts_with_n_exons(c))
    out_dict = {}

    for k, v in gff_tree.items():
        if k and feature_types[k] == PARENT_TYPE:
            exon_count = len([x for x in v if "CDS" in x])
            utr_count = len([x for x in v if "utr" in x])

            if utr_count >= 2:
                label = "{}_UTRs".format(exon_count)
                if label in out_dict:
                     out_dict[label].append(k)
                else:
                     out_dict[label] = [k]
            elif utr_count == 1:
                print("{} has irregular UTR count, ignoring...".format(k))
            else:
                label = "{}_NO_UTRs".format(exon_count)
                if label in out_dict:
                     out_dict[label].append(k)
                else:
                     out_dict[label] = [k]

    # add mRNAs by chromosome
    chroms = set(gff_df['seq_id'].to_list())
    features = set(gff_df['type'].to_list())

    for c in chroms:
        all_items_in_chrom = gff_df[(gff_df['seq_id'] == c)]['ID'].to_list()
        out_dict[c] = all_items_in_chrom

        for f in features:
            items_in_chrom = gff_df[(gff_df['type'] == f) & (gff_df['seq_id'] == c)]['ID'].to_list()
            out_dict["{}_{}".format(c, f)] = items_in_chrom

    # add mRNAs by length
    # should first plot read length distribution
    transcript_lengths = [0, 500, 1000, 4000, 10000, 100000]

    for l in range(len(transcript_lengths) - 1):
        out_dict["mRNAs_{}-{}bp".format(transcript_lengths[l], transcript_lengths[l+1])] = []


    mRNAs = gff_df[(gff_df['type'] == "mRNA")]
    for m_index, m in mRNAs.iterrows():
        tx_length = m['end'] - m['start']
        
        for l in range(len(transcript_lengths)):
            if tx_length < transcript_lengths[l]:
                out_dict["mRNAs_{}-{}bp".format(transcript_lengths[l-1], transcript_lengths[l])].append(m['ID'])
                break

    return out_dict



# TODO write separate files of transcripts with exon counts to files
def to_file(out_dict, filename):
    with open(filename, 'w') as f:
        json.dump(out_dict, f)

# plot that shows exon counts as bar graphs
def plot_exon_counts(o):
    sort_strat = lambda e: int(e.split('_')[0])
    utr_keys = sorted([k for k in o.keys() if "NO_UTRs" in k], key = sort_strat)
    no_utr_keys = sorted([k for k in o.keys() if "NO_UTRs" not in k], key = sort_strat)

    max_exon_count = max([int(x.split('_')[0]) for x in utr_keys + no_utr_keys])

    utr_exon_counts = {}
    no_utr_exon_counts = {}

    for c in range(max_exon_count + 1):
        no_utr_label = "{}_NO_UTRs".format(c)
        utr_label = "{}_UTRs".format(c)

        if no_utr_label in o.keys():
            no_utr_exon_counts[c] = len(o[no_utr_label])
        else:
            no_utr_exon_counts[c] = 0

        if utr_label in o.keys():
            utr_exon_counts[c] = len(o[utr_label])
        else:
            utr_exon_counts[c] = 0

    d = {
        "mRNA - annotated UTR": utr_exon_counts.values(),
        "mRNA - no annotated UTR": no_utr_exon_counts.values()
    }
    df = pandas.DataFrame(d)
    tick_size = 5
    c = (math.ceil(max_exon_count / tick_size) + 1) * tick_size

    df.plot.bar(xticks=range(0, c, tick_size), width=0.8)

    for i in range(max_exon_count):
        print("{}: {}".format(i, utr_exon_counts[i] + no_utr_exon_counts[i]))

    # plt.yscale('log',base=10)
    plt.ylabel("number of transcripts")
    plt.xlabel("exon count")
    plt.legend(loc="upper right")

    plt.show()

o = build_outdict()
to_file(o, OUTPUT_FILE)

# plot_exon_counts(o)