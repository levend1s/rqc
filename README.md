## RQC (RNA Quality Control)

RQC helps you analyse the quality of your BAM files. Input BAM files must be sorted and indexed. RQC can do the following:

* Summarise read PHRED scores for one or multiple BAM files
* Summarise read MAPQ scores for one or multiple BAM files
* Generate transcript coverage plots
    * Accepts bedmethyl and BAM files for coverage plots
    * Can specify padding
    * Can show coverage for: 
        * Genomic regions
        * Genomic regions split by UTR and CDS
        * Subfeatures (Exon and UTR boundaries) NOTE all input genes must have the same number of subfeatures for this to work accurately.
    * Can show frequncy plot of DRACH sites (locations specified by bed file)
    * Can generate normalised coverage plots, which average normalised individual gene coverage
    * Can show raw proportion of mod coverage against read depth, or relative abundance
    * Can calculate locations of 'mod peaks' if given a read depth and mod proportion threshold. If mod abundance exceeds this threshold it's considered a peak, or 'cannonical modification'.
    * TODO: separate_plots_by_label_prefix
    * TODO: specify colours in samples.txt
    * TODO: coverage smoothing


* Mapping % bar graphs against multiple genomes
* Sequence logos
* Violin plots for read lengths


## SETUP
```
# set up and enter virtual environment
$ python3 -m venv .venv
$ source .venv/bin/activate
$ pip install -r requirements.txt

or

$ virtualenv env
$ source env/bin/activate
$ pip install -r requirements.txt

# run rqc
$ python rqc.py -q 7 -m 5 plot "36C1 Pfal" /Users/joshualevendis/Downloads/36C1_to_pfal.sorted.bam

# exit the virtual environment
$ deactivate
```
If you have a file that looks like this:
samples.txt
```
36C1_Pfal /Users/joshualevendis/Downloads/36C1_to_pfal.sorted.bam 
36C1_Yeast /Users/joshualevendis/Downloads/36C1_to_yeast_sorted.bam 
36C1_Human /Users/joshualevendis/Downloads/36C1_to_humans.sorted.bam
```
You can run something like this:
```
$ python rqc.py -q 7 -m 5 plot $(cat samples.txt)
```
![Screenshot 2024-09-21 at 10 36 34 PM](https://github.com/user-attachments/assets/9eab357d-b738-48aa-b24b-7f32687c2180)


After playing around with quality filters, you can then print out the read ids of reads which passed the filter:
```
$ python rqc.py --check_duplicate_reads True -n 10 --sort_by phred_scores -q 7 -m 5 -l 100000 search $(cat samples.txt)

searching...
36C1_Pfal
Empty DataFrame
Columns: [read_id, phred_scores, mapq_scores, template_lengths]
Index: []
36C1_Yeast
                                read_id  phred_scores  mapq_scores  template_lengths
0  ed3f5afc-5154-42aa-b71c-ebfaeb361bbc          7.08            8            232367
36C1_Human
Empty DataFrame
Columns: [read_id, phred_scores, mapq_scores, template_lengths]
Index: []
looking for read intersect in 36C1_Pfal (0) and 36C1_Yeast (1)
looking for read intersect in 36C1_Pfal (0) and 36C1_Human (0)
looking for read intersect in 36C1_Yeast (1) and 36C1_Human (0)
FOUND 0 READS COMMON TO SAMPLES:

sample	count	total_bases_read	mean	median	Q1	Q3	min	max	N10	N50	N90	
36C1_Yeast	1	232367	232367	232367	232367.0	232367.0	232367	232367	232367	232367	232367	
```


You might then like to find information on a specific read
```
$ python rqc.py inspect ed3f5afc-5154-42aa-b71c-ebfaeb361bbc $(cat samples.txt)
```

## Usage

# plot_coverage

The plot_coverage command takes in an annotation file, a bam file and a list of gene/transcript ids to generate coverage plots for. It generates a plot with total coverage and normalised coverage against % through gene, 0% 5' -> 3' 100%. Normalised coverage means the sequencing depth is scaled from a true number (like 1000 bases) to a number between 0 and 1. This reduces the impact a highly expressed gene has on the curve and gives a less bias comparison of coverage across the selected genes.

```
# By type...
python rqc.py -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff plot_coverage tRNA $(head -n 1 samples.txt)

# By list of gene ids...
python rqc.py -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff plot_coverage "['PF3D7_0209800']" $(head -n 1 samples.txt)

# if you have a json file you can extract an entry (must be list of gene/transcript ids) with jq in a subshell
python rqc.py -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff plot_coverage "$(cat /Users/joshualevendis/Documents/RNA/rqc/pfal_mRNA_exon_counts.json | jq -r '."25"')" $(head -n 1 samples.txt)
```
