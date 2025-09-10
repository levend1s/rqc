# SETUP
```
# set up and enter virtual environment
$ python3 -m venv .venv
$ source .venv/bin/activate
$ pip install -r requirements.txt

or

$ python3 -m virtualenv env
$ source env/bin/activate
$ pip install -r requirements.txt

# run rqc
$ python rqc.py -q 7 -m 5 plot "36C1 Pfal" /Users/joshualevendis/Downloads/36C1_to_pfal.sorted.bam

# exit the virtual environment
$ deactivate
```

# Usage

## plot_coverage

The plot_coverage command takes in an annotation file, a bam file and a list of gene/transcript ids to generate coverage plots for. It generates a plot with total coverage and normalised coverage against % through gene, 0% 5' -> 3' 100%. Normalised coverage means the sequencing depth is scaled from a true number (like 1000 bases) to a number between 0 and 1. This reduces the impact a highly expressed gene has on the curve and gives a less bias comparison of coverage across the selected genes.

```
python rqc.py plot_coverage -m subfeature -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 1000 --read_depth 0 --ids PF3D7_1123900.1 -i ./samples/laptop_all_mods_samples.txt --separate_y_axes

python rqc.py plot_coverage -m subfeature_cds -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --padding 100 --padding_ratio 0.1 --bins 100 --read_depth 0 --ids PF3D7_1123900.1 -i ./samples/laptop_all_mods_samples.txt
```

The input file laptop_all_mods_samples.txt should look something like this:
```
inosine_read_depth control bam /Users/joshlevendis/Downloads/bams/28C1_to_pfal.50MAPQ.sorted.bam
inosine_inosine control bedmethyl /Users/joshlevendis/Downloads/bedmethyls/28C1_to_pfal_0.65.0p.17596.bed

m6A_read_depth control bam /Users/joshlevendis/Downloads/bams/28C1_to_pfal.50MAPQ.sorted.bam
m6A_m6A control bedmethyl /Users/joshlevendis/Downloads/bedmethyls/28C1_to_pfal_0.65.0p.a.bed

pseudouridine_read_depth control bam /Users/joshlevendis/Downloads/bams/28C1_to_pfal.50MAPQ.sorted.bam
pseudouridine_pseudouridine control bedmethyl /Users/joshlevendis/Downloads/bedmethyls/28C1_to_pfal_0.65.0p.17802.bed

m5C_read_depth control bam /Users/joshlevendis/Downloads/bams/28C1_to_pfal.50MAPQ.sorted.bam
m5C_m5C control bedmethyl /Users/joshlevendis/Downloads/bedmethyls/28C1_to_pfal_0.65.0p.m.bed
```


## motif_finder

```
# find all motifs of TTTN or NGG, where N=anything
python rqc.py motif_finder -a CpBGF_genome_V1.gff -g CpBGF_genome_v1.fasta -o TTTN_crypto_bgf.tsv -m "TTT."
python rqc.py motif_finder -a CpBGF_genome_V1.gff -g CpBGF_genome_v1.fasta -o TTTN_crypto_bgf.tsv -m ".GG"

# find all stop codons within CDS regions (stop codon is TAG, TAA, TGA, must be at the end of the region ($)) or start codons
python rqc.py motif_finder -a CpBGF_genome_V1.gff -g CpBGF_genome_v1.fasta -o stop_codons_crypto_bgf.tsv -m "(TA[AG]|TGA)$" -f CDS
python rqc.py motif_finder -a CpBGF_genome_V1.gff -g CpBGF_genome_v1.fasta -o start_codons_crypto_bgf.tsv -m "^ATG" -f CDS


# get the list of gene ids from the crypto gff file (parsing is annotation specific)
grep protein_coding CpBGF_genome_V1.gff | awk '$3=="gene" {split ($9,x,/[-;]/); print x[2]}' > pcg_list_crypto_bgf.tsv

# based on the 7th column, and sorted by version, keep only the first entry for start codons (sort -V) or last entry for stop codons (sort -rV). These should be the true start and stop codons.
awk '{split ($7,x,/[-]/); print $0"\t"x[1]}' start_codons_crypto_bgf.tsv | sort -V -k7 - | awk '!seen[$8]++' > start_codons_crypto_bgf.filtered.tsv
awk '{split ($7,x,/[-]/); print $0"\t"x[1]}' stop_codons_crypto_bgf.tsv | sort -rV -k7 - | awk '!seen[$8]++' > stop_codons_crypto_bgf.filtered.tsv

# print any gene ids missing from the list of start codons 
while read line; do grep -q "$line" start_codons_crypto_bgf.filtered.tsv || echo $line; done < pcg_list_crypto_bgf.tsv
while read line; do grep -q "$line" stop_codons_crypto_bgf.filtered.tsv || echo $line; done < pcg_list_crypto_bgf.tsv
```

## calculate_offsets

```
python rqc.py calculate_offsets -d 100 -r pam_analysis/start_codons_crypto_bgf.filtered.tsv -o pam_analysis/crypto_bgf_start_offsets.tsv -i NGG pam_analysis/NGG_crypto_bgf.tsv TTTN pam_analysis/TTTN_crypto_bgf.tsv
```

## plot_relative_distance
```
python rqc.py plot_relative_distance -l "start codon" -d 50 -i pam_analysis/crypto_bgf_start_offsets.tsv -o crypto_start.eps
```

## sequence_logo
```
python rqc.py sequence_logo -a PlasmoDB-67_Pfalciparum3D7.gff -g PlasmoDB-67_Pfalciparum3D7_Genome.fasta -l 1 -p 2 -i drach_m6A_sites_28u32u36.bed

# if the bed file is 1 indexed, can adjust
python rqc.py sequence_logo -a ~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff -g ~/Downloads/Pfalciparum3D7/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta -l 3 -p 0 --adjust -1 -i pam_analysis/start_codons_plasmo.tsv
```

## approximate_tes
```
python rqc.py approximate_tes -a ~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff -i samples/samples.txt --ids PF3D7_1123900.1

python rqc.py approximate_tes -a ~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff -i samples/samples.txt --ids $(jq -r '.Pf3D7_01_v3_mRNA[]' ./pfal_mRNA_exon_counts.json)
```

## gene methylation analysis

This calculates average methylation (m6A) for genes. It calculates separates statistic

```
python rqc.py gene_methylation_analysis -a ~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff -i samples/samples.txt --ids PF3D7_1123900.1 -d 20 -r 0.7
python rqc.py gene_methylation_analysis -a ~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff -i samples/samples.txt --ids PF3D7_1338200.1 -d 20 -r 0.7
time python rqc.py gene_methylation_analysis -a ~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff -i samples/samples.txt --ids $(jq -r '.Pf3D7_01_v3_mRNA[]' ./pfal_mRNA_exon_counts.json) -d 10 -r 0.5
```


python rqc.py gene_methylation_analysis -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff -i laptop_samples.txt --type mRNA -d 10 -r 0.5 --compare_methylation_between_treatments 28C 28K -p 500 -o mrna_gene_methylation_analysis.tsv --exclude_contigs Pf3D7_API_v3 Pf3D7_MIT_v3

## filter_bam_by_mod
```
python rqc.py filter_bam_by_mod --include 213060 -s + -c Pf3D7_02_v3 -i laptop_samples.txt
```
