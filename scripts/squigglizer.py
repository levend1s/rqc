from ont_fast5_api.fast5_interface import get_fast5_file
import gzip
import numpy

import h5py
import pysam
from matplotlib import pyplot as plt

fast5_filepath = "/Users/joshualevendis/FAQ59548_pass_2fc928ab_0.fast5" # This can be a single- or multi-read file
dorado_stdout = "/Users/joshualevendis/dorado_stdout.bam"
nanopolish_ks1_polya_results = "/Users/joshualevendis/ks1_polya_results.tsv"

hdf = h5py.File(fast5_filepath, 'r')

dorado_stdout_bamfile = pysam.AlignmentFile(dorado_stdout, "rb")
dorado_read = None
for x in dorado_stdout_bamfile.fetch(until_eof=True):
    dorado_read = x

dorado_read_tags = dict(dorado_read.tags)
dorado_move_table = dorado_read_tags['mv'] # array

def move_table_to_signal_space(table, signal_start):
    stride = table[0]
    events = []
    head = signal_start
    events.append(head)

    for s in table[1:]:
        head += stride

        if s == 1:
            events.append(head)
            
    return events

#read_id = "00088e28-7499-4a85-afaf-5c7ae916cedc"
#read_id = "000e19b4-e5de-430e-9363-62b0e50d89bf"

# readname	                            contig	    position	leader_start	adapter_start	polya_start	transcript_start	read_rate	polya_length	qc_tag
# 0068d8ef-a397-4e80-89be-520bb339472e	Pf3D7_05_v3	125081	    2.0	            3.0	            5096.0	    5434.0	            130.96	    9.65	        PASS
#read_id = "0068d8ef-a397-4e80-89be-520bb339472e"

# ks1 nanopolish polya:
# 02c953bc-207b-4102-9ba6-4bc7fcc81330	Pf3D7_04_v3	668314	    2.0	            3.0	            7040.0	    10654.0	            100.40	    115.43	        PASS
# ks1 bam file:
# 02c953bc-207b-4102-9ba6-4bc7fcc81330    0       Pf3D7_04_v3     668315  60      8M1I16M1D38M2D11M1D34M1D31M4D15M1D21M3D17M82S       *       0       0       TTTAAATATTTTTAATTTGAATAATAAACATTTAAATTTGTTATAAAGAGGTGTTACAAGATCAAAAAAATAATAAAAAAAGAATTAAATAAAAAAAATTTTAAAATTAAAATATTAAATATATATATATATATATATAATATCCATATTGACATTTTTAAAATGCAATTATTTAAAAAGAATGACAAAAAACAGAATACCCCATCCCTACACATAACCCTATCATTTATATCCATATATATCTTACATATATGCCCTCCTACCATCCTCTAAT      -5:;=?<4;>AC?983.<4:;2711:<<(%&+9:;;58;.79=??:;=3551-+%#'*72==68==56:(8;;<973497A7=EA87($7988<;>;*7?B@8?@C:=8<:=73?<:8%5'8*87=$0'1'35,.70.76,?>58.7&37212019EG?7;@/=226<216%4370//723<56::64:4-.++$$%%#(*'%$$')+$$''%$%%%&,/*1.)))&&%'%)+12)&'#&,*.-*$%$%)33+$%(+%*/&$+$&%'($)',-(      NM:i:20 ms:i:146        AS:i:143    nn:i:0  ts:A:+  tp:A:P  cm:i:19 s1:i:127        s2:i:0  de:f:0.0704     rl:i:42
read_id = "02c953bc-207b-4102-9ba6-4bc7fcc81330"
leader_start = 2.0
adapter_start = 3.0
poly_start = 7040.0
transcript_start = 10654.0
read_rate = 100.40
polya_length = 115.43

# https://hasindu2008.github.io/slow5specs/fast5_demystified.pdf
# http://simpsonlab.github.io/2015/04/08/eventalign/
# https://bioinformatics.stackexchange.com/questions/18931/how-can-i-get-the-duration-of-the-events-in-nanopore
# https://samtools.github.io/hts-specs/SAMv1.pdf
# https://www.biorxiv.org/content/10.1101/2024.02.19.581111v2.full.pdf

# https://timd.one/blog/genomics/file_formats.php

# dorado output 
# https://github.com/nanoporetech/dorado/blob/36218004593b0ebba37ff870eb3d92c8238fb0c6/documentation/SAM.md

# CIGAR string
# https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/
# is 0-based
# M	Match 		Exact match of x positions
# N	Alignment gap 	Next x positions on ref don’t match
# D	Deletion 	Next x positions on ref don’t match
# I	Insertion 	Next x positions on query don’t match


# TODO draw event lines on graph

# NOTE: the fastq sequence is the the bam file complement sequence in reverse. Why? Is this true for all?

base_colors = {
    "A": "green",
    "U": "red",
    "G": "orange",
    "C": "blue"
}

def view_read_h5py(read_id):
    read = hdf['read_{}'.format(read_id)]

    raw_fastq = read['Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]
    fastq = raw_fastq.decode('UTF-8').split()
    basecalled_string = fastq[-3] # magic number
    basecall_qual = fastq[-1] # magic number


    basecall_summary = read['Analyses/Basecall_1D_000/Summary/basecall_1d_template/']
    basecall_location = basecall_summary.attrs['basecall_location']
    mean_qscore = basecall_summary.attrs['mean_qscore']
    num_events = basecall_summary.attrs['num_events']
    sequence_length = basecall_summary.attrs['sequence_length']
    block_stride = basecall_summary.attrs['block_stride']
    
    segmentation_summary = read['Analyses/Segmentation_000/Summary/segmentation/']
    duration_template = segmentation_summary.attrs['duration_template']
    first_sample_template = segmentation_summary.attrs['first_sample_template']
    raw_signal = read['Raw/Signal'][()]


    digitisation = read['channel_id'].attrs['digitisation']
    offset = read['channel_id'].attrs['offset']
    read_range = float("{0:.2f}".format(read['channel_id'].attrs['range']))
    sampling_rate = read['channel_id'].attrs['sampling_rate'] # number of samples read per second

    scaled_signal = (raw_signal + offset) * (read_range / digitisation)
    signal_min = min(scaled_signal)
    signal_max = max(scaled_signal)

    fig, axes = plt.subplots(nrows=3, ncols=1)
    axes[0].plot(scaled_signal)
    axes[0].add_patch(plt.Rectangle(xy=(first_sample_template,signal_min), width=duration_template, height=signal_max, color='green', alpha=0.3))
    axes[0].axvline(x = poly_start, color = "black", linestyle = '-')
    axes[0].axvline(x = leader_start, color = "black", linestyle = '-')
    axes[0].axvline(x = adapter_start, color = "black", linestyle = '-')
    axes[0].axvline(x = transcript_start, color = "black", linestyle = '-')

    
    qual = list(map(ord, list(basecall_qual)))
    #qual score are stored ASCII 33-126 (ie 33 indicates quality score of 0)
    qual = [q-33 for q in qual]

    # "a Phred quality score is a measure of confidence based on error rate, it is stored alongside fastq
    # phred_qscore = -10 x log(Pe), where Pe is the estimated probability of error" 
    #   - https://help.nanoporetech.com/en/articles/6629615-where-can-i-find-out-more-about-quality-scores
    # To calculate mean q score, you must convert the Phred back to probabilities, find average, then convert
    # back to phred - https://www.jmferguson.science/post/split-nanopore-reads-by-qscore 
    # mean_err = numpy.exp(numpy.array(qual) * (-numpy.log(10) / 10.)).mean()
    # score = -10 * numpy.log10(max(mean_err, 1e-4))
    # score should be close enough to mean_qscore
    # Taking the average of Phred scores will yield a score skewed higher than mean_qscore
    # The lower the error rate, the higher the Phred score.

    axes[1].plot(qual)
    print(mean_qscore)
    print(sequence_length)
    print(len(basecalled_string))
    axes[1].set_xticks(list(range(0, (len(basecalled_string)))))
    axes[1].xaxis.set_ticks_position('none')
    axes[1].set_xticklabels(list(basecalled_string), fontsize=5)

    axes[1].axhline(y = 7, color = "red", linestyle = '-')
    axes[1].axhline(y = float(mean_qscore), color = "green", linestyle = '-')

    # add event lines for dorado basecalled information
    dorado_events = move_table_to_signal_space(dorado_read_tags['mv'], 0)
    for e in dorado_events:
        axes[1].axvline(x = e, color = "black", linestyle = '-')


    errors = numpy.exp(numpy.array(qual) * (-numpy.log(10) / 10.))
    axes[2].plot(errors)
    axes[2].set_xticks(list(range(0, (len(basecalled_string)))))
    axes[2].xaxis.set_ticks_position('none')
    axes[2].set_xticklabels(list(basecalled_string), fontsize=5)

    for l in axes[1].get_xticklabels():
        l.set_color(base_colors[l._text])
    for l in axes[2].get_xticklabels():
        l.set_color(base_colors[l._text])

    plt.show()


view_read_h5py(read_id)
