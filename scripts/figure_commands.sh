# Basecalling POD5 + uBAM alignment

# I merged my previous two slurms into a single super slurm. With dorado 5.1.0 you can't call m6A_DRACH and inosine_m6A at the same time since they have the same cannonical base, but calling inosine_m6A does both (although this may not be the case for DNA?). I then submitted:

# TODO: compartmentalise scripts into STEPs and folders

#!/bin/bash

#SBATCH --nodes=1
#SBATCH -p gpu-h100
#SBATCH --gres=gpu:1
#SBATCH --time=0-12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=1

SAMPLE="28C1"
DORADO="/data/gpfs/projects/punim0875/josh/dorado-0.8.3-linux-x64/bin/dorado"
RNASEQ_PROJ_DIR="/data/gpfs/projects/punim0875/josh/projects/proj-5210_ralphlab-1128.4.524/rnaseq/nanopore/METTL3_KS_TIME_SERIES"
READS="${RNASEQ_PROJ_DIR}/${SAMPLE}/20240826_1458_MN24403_FAZ60852_d9c3cda9/pod5/"
OUTDIR="/data/gpfs/projects/punim0875/josh/dorado_0.8.3_basecalled/"
UBAM="${OUTDIR}/${SAMPLE}.bam"

PFAL_FASTA="${PROJ_DIR}/references/Pfalciparum3D7/fasta/data/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
HUMAN_FASTA="${PROJ_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
YEAST_FASTA="${PROJ_DIR}/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"

$DORADO basecaller --emit-moves --estimate-poly-a --verbose sup,pseU,m5C,inosine_m6A $READS > $UBAM

$DORADO aligner $PFAL_FASTA $UBAM > "${OUTDIR}/${SAMPLE}_to_pfal.bam"
$DORADO aligner $HUMAN_FASTA $UBAM > "${OUTDIR}/${SAMPLE}_to_humans.bam"
$DORADO aligner $YEAST_FASTA $UBAM > "${OUTDIR}/${SAMPLE}_to_yeast.bam"
# And monitored it with:
# sbatch nanopore_basecall_and_align.slurm
# squeue --me


# Filtering BAMs:
SAMPLE=36K1 FILEPATH=/home/ubuntu/mediaflux_mnt/rnaseq/nanopore/METTL3_KS_TIME_SERIES_DORADO_8.3/MODBAMS; samtools view -bq 50 ${FILEPATH}/${SAMPLE}_to_pfal.bam > ${FILEPATH}/${SAMPLE}_to_pfal.50MAPQ.bam && samtools sort ${FILEPATH}/${SAMPLE}_to_pfal.50MAPQ.bam > ${FILEPATH}/${SAMPLE}_to_pfal.50MAPQ.sorted.bam && samtools index ${FILEPATH}/${SAMPLE}_to_pfal.50MAPQ.sorted.bam && samtools quickcheck ${FILEPATH}/${SAMPLE}_to_pfal.50MAPQ.sorted.bam && echo "filtered all ok"


# Generating bedmethyls:
SAMPLE=32C1 DIR=/home/ubuntu/mediaflux_mnt/rnaseq/nanopore/METTL3_KS_TIME_SERIES_DORADO_8.3/MODBAMS THRESHOLD=0.95; ~/mount_rna/dist_modkit_v0.4.2_10d99bc/modkit pileup --filter-threshold ${THRESHOLD} ${DIR}/${SAMPLE}_to_pfal.50MAPQ.sorted.bam ${SAMPLE}_to_pfal_${THRESHOLD}.bed --log-filepath ${SAMPLE}_to_pfal_${THRESHOLD}.log
TYPE="a" PERCENT_CUTOFF=0 SAMPLE=32C1; cat ${SAMPLE}_to_pfal_0.95.bed | awk -v TYPE=$TYPE -v PERCENT_CUTOFF=$PERCENT_CUTOFF '{if ($4 == TYPE && $11 > PERCENT_CUTOFF) print $0 }' > ${SAMPLE}_to_pfal_0.95.${PERCENT_CUTOFF}p.${TYPE}.bed


# Counting features:
DIR=/home/ubuntu/mediaflux_mnt/rnaseq/nanopore/METTL3_KS_TIME_SERIES_DORADO_8.3/MODBAMS; ~/mount_rna/subread-2.0.6-Linux-x86_64/bin/featureCounts -s 1 -L -R CORE -a ~/mount_rna/Pfalciparum3D7/gff/data/PlasmoDB-68_Pfalciparum3D7.gff -o ./8.3_featureCounts ${DIR}/28C1_to_pfal.bam ${DIR}/28C2_to_pfal.bam ${DIR}/28K1_to_pfal.50MAPQ.sorted.bam ${DIR}/28K2_to_pfal.bam ${DIR}/32C1_to_pfal.merged.bam ${DIR}/32C2_to_pfal.bam ${DIR}/32K1_to_pfal.merged.bam ${DIR}/32K2_to_pfal.bam ${DIR}/36C1_to_pfal.bam ${DIR}/36C2_to_pfal.bam ${DIR}/36K1_to_pfal.bam ${DIR}/36K2_to_pfal.bam

# TES change analysis (WAM change analysis)
python rqc.py approximate_tes -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --poly_a_filter 10 -i laptop_samples.txt --type mRNA -d 10 --compare_apa_between_treatments 36C 36K -p 500 -o mrna_tes_analysis_36hpi.tsv --exclude_contigs Pf3D7_API_v3 Pf3D7_MIT_v3
python rqc.py gene_methylation_analysis -a $ANNOTATION -i laptop_samples.txt --type mRNA -d 10 -r 0.5 --compare_methylation_between_treatments 36C 36K -p 500 --poly_a_filter 10 -o output/mRNA_methylation_analysis_36hpi.tsv --exclude_contigs Pf3D7_API_v3 Pf3D7_MIT_v3

# ------------ Fig 2 ------------ #

# A: Mod coverage
python rqc.py plot_coverage -m subfeature -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 1000 --read_depth 0 --ids PF3D7_1123900.1 -i ./samples/laptop_samples_all_mods_coverage.txt --separate_y_axes -o PF3D7_1123900_mod_compare.eps  --plot_type bar

# generate file of nuclear mRNA regions only
awk '($3 == "mRNA" && $1 != "Pf3D7_MIT_v3" && $1 != "Pf3D7_API_v3") {print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' ${ANNOTATION} > nuclear_mRNA_regions.bed 
modkit summary --filter-threshold 0.95 --include-bed nuclear_mRNA_regions.bed --no-sampling ~/Downloads/bams/28C1_to_pfal.50MAPQ.sorted.bam
modkit summary --filter-threshold 0.65 --sampling-frac 0.1 --include-bed nuclear_mRNA_regions.bed --interval-size 100000 ~/Downloads/bams/28C1_to_pfal.50MAPQ.sorted.bam


# B: IGV screenshot of PF3D7_1123900.1

# C: PCA (transcript count)
# D: PCA top 50 loadings (transcript count)
wam_pca.R

# E: Differential abundance (transcript count)
differential_abundance.R

# ------------ Fig 3 ------------ #

# A: Gene body coverage (28C1, 28K1)
python rqc.py plot_coverage -m gene -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 1000 --read_depth 0 --padding 0 --padding_ratio 0.1 --ids PF3D7_1123900.1 -i ./samples/laptop_samples.txt -o PF3D7_1338200.1.eps --mod_normalisation other
# IGV screenshot, hide small indels (<100bp), show m6A > 0.95p

# B: coverage (top PC1 contributors)
# the trick for padding here is padding_ratio = padding / (gene length + padding)

python rqc.py plot_coverage -m gene -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 1000 --read_depth 0 --padding 0 --padding_ratio 0.1 --ids PF3D7_1237700.1 -i ./samples/laptop_samples.txt -o PF3D7_1237700.1.eps --mod_normalisation other
python rqc.py plot_coverage -m gene -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 1000 --read_depth 0 --padding 300 --padding_ratio 0.12 --ids PF3D7_0309600.1 -i ./samples/laptop_samples_coverage_plot.txt --mod_normalisation other -o PF3D7_0309600.1.eps
python rqc.py plot_coverage -m gene -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 1000 --read_depth 0 --padding 300 --padding_ratio 0.18 --ids PF3D7_1128200.1 -i ./samples/laptop_samples_coverage_plot.txt --mod_normalisation other -o PF3D7_1128200.1.eps

# C: Differential methylation of individual sites
diff_meth_site_advanced.R
# ONLY 28 HPI SAMPLES IN FIGURE
# ONLY GENES WITH NO OVERLAPS (0BP)
# FILTER ALL READS >= 10
# FILTER ALL SITES WITH >= 20 M6A IN ALL CONTROL SAMPLES 
# FILTER OUT DE SITES (NB GLM)
# DM EDGER GROUP * METH_STATUS (NB GLM)
# FDR SIGNIFICANCY CUTOFF 0.2

# D: WAM change
wam_change_t_test.R

# ------------ Fig 4 ------------ #
# A: cartoon

# FINDING m6A COUNTS
bash glori_x_nanopore_analysis.sh

# B: venn diagrams
# glori_nanopore_venn.R

# C: sequence logos
python rqc.py sequence_logo -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff -g ~/Documents/RNA/honours/Pfalciparum3D7/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta -l 1 -p 2 -i glori_x_nanopore_analysis/nanopore_28u32u36.bed -o ont_union_15d25p.eps
python rqc.py sequence_logo -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff -g ~/Documents/RNA/honours/Pfalciparum3D7/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta -l 1 -p 2 -i glori_x_nanopore_analysis/glori_12u24u48.bed -o glori_union_15d25p.eps

# D: IGV screenshot of PF3D7_0102300, loaded 28C1, glori_12i24i48.bed, nanopore_28u32u36.bed and DRACH.bed

# E: metagene coverage

# input samples file looks like this:
# ONT_total_m6A control bed ./glori_x_nanopore_analysis/nanopore_28u32u36.bed
# GLORI-seq_total_m6A control bed ./glori_x_nanopore_analysis/glori_12u24u48.bed

python rqc.py plot_coverage -m subfeature_cds -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 100 --read_depth 0 --padding 500 --padding_ratio 0.1 --type mRNA -i ./samples/samples_m6A_sites.txt --mod_normalisation other --skip_malannotations --plot_density -o nanoporeXglori15d25p.eps

# F: nanopore offset from canonical PAS
awk '$6 != 0 {print $1"\t"$2"\t"$6"\t"$6+1"\t"$3"\t"0"\t"$4}' '/Users/joshualevendis/rqc/output/mRNA_tes_analysis_28hpi_compare.tsv' > 28hpi_canonical_pas.tsv
awk '$6 != 0 {print $1"\t"$2"\t"$5"\t"$5+1"\t"$3"\t"0"\t"$4}' '/Users/joshualevendis/rqc/output/mRNA_tes_analysis_28hpi_compare.tsv' > 28hpi_3p_annotation.tsv
# must manually add headers to above output files
python rqc.py calculate_offsets -d 1000 -r ./28hpi_canonical_pas.tsv -o 3p_m6a_offsets_from_canonical_pas.tsv -i m6a /Users/joshualevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/rqc_output/15depth25prop/28u32u36.bed annotation_three_p ./28hpi_3p_annotation.tsv
python rqc.py plot_relative_distance -l "canonical PA" -d 1000 -i 3p_m6a_offsets_from_canonical_pas.tsv -o 3p_m6a_offsets.eps

# G: m6A specific analysis
python rqc.py plot_coverage -m gene -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 1000 --read_depth 0 --padding 100 --padding_ratio 0.1 --ids PF3D7_0402000.1 -i ./samples/laptop_samples_coverage_plot.txt -o PF3D7_0402000.1.eps --mod_normalisation other --plot_type bar --alpha 1
python rqc.py m6A_specific_tes_analysis -a $ANNOTATION -i samples/laptop_samples.txt --ids PF3D7_0402000.1 --poly_a_filter 10 -d 10 --separate_mod_tracks --offset_padding 1000 --offset 114415

# ------------ Fig 5 ------------ #
# A:
# IGV 28C1, 28K1 UBC E2 (PF3D7_1203900)

# B: cartoon

# C: 
differential_abundance.R

# D:
run_on_t_test.R

# E:
wam_pca.R

# F: 
wam_pca.R

# G: neighbour gene coverage
python rqc.py plot_coverage -m gene -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 1000 --read_depth 0 --padding 1823 --padding_ratio 0.25 --ids PF3D7_1439700 -i ./samples/laptop_samples_coverage_plot.txt -o PF3D7_1439700.1.eps --mod_normalisation other --plot_type bar --alpha 1
python rqc.py plot_coverage -m gene -a ~/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff --bins 1000 --read_depth 0 --padding 2106 --padding_ratio 0.33 --ids PF3D7_0518300 -i ./samples/laptop_samples_coverage_plot.txt -o PF3D7_0518300.1.eps --mod_normalisation other --plot_type bar --alpha 1

# ------------ Supplementary figure 1 ------------ #

# RNA-seq library sizes (steve)
# Read length distributions (steve)


# PCA (canonical wam %)
# Heatmap (canonical wam %)
wam_pca.R

# ------------ Fig 6 ------------ #

# From each featureCounts file the third column is number of assignments. We can filter for just entries with multiple assignments (>1) and get a list of gene ids that are involved in overlapping assignments. The below gives us a list of genes with the number of overlapping read counts:

SAMPLE=28C1; awk '$3 > 1 {print $4}'  ~/Downloads/featureCountsStrandedOverlapMAPQFilter/${SAMPLE}_to_pfal.50MAPQ.sorted.bam.featureCounts | tr ',' '\n' > ${SAMPLE}_overlaps.txt
SAMPLE=28C1; sort ${SAMPLE}_overlaps.txt | uniq -c > ${SAMPLE}_overlaps_counts.txt
# There are 3 parent types of features in the Plasmodium off file: protein_coding_gene, ncRNA_gene, pseudogene. I'll create a list of all gene IDs with the following commands:

cat ~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff | awk '$3=="protein_coding_gene" {split ($9,x,/[=;]/); print x[2]}' > pcg_pseudogene_ncRNA_list.tsv
cat ~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff | awk '$3=="ncRNA_gene" {split ($9,x,/[=;]/); print x[2]}' >> pcg_pseudogene_ncRNA_list.tsv
cat ~/Downloads/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff | awk '$3=="pseudogene" {split ($9,x,/[=;]/); print x[2]}' >> pcg_pseudogene_ncRNA_list.tsv
wc -l pcg_pseudogene_ncRNA_list.tsv (=5720)
# Now I'll combine my overlap count files into a single file:

awk 'NR==FNR{a[$2]=$1;next} {if ($1 in a) {print $1"\t" a[$1]} else {print $1"\t"0}}' 28C1_overlaps_counts.txt pcg_pseudogene_ncRNA_list.tsv > 28hpi_overlap_counts.txt
SAMPLE=28KC2; awk 'NR==FNR{a[$2]=$1;next} {if ($1 in a) {print $0"\t" a[$1]} else {print $0"\t"0}}' ${SAMPLE}_overlaps_counts.txt 28hpi_overlap_counts.txt > temp.txt && mv temp.txt 28hpi_overlap_counts.txt