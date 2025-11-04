#!/usr/bin/env bash

# Define your input files
FILES=(
    "/Users/joshualevendis/Downloads/a/28C1_to_pfal_0.95.0p.a.bed" 
    "/Users/joshualevendis/Downloads/a/28C2_to_pfal_0.95.0p.a.bed" 
    "/Users/joshualevendis/Downloads/a/28K1_to_pfal_0.95.0p.a.bed" 
    "/Users/joshualevendis/Downloads/a/28K2_to_pfal_0.95.0p.a.bed"
)

ANNOTATION="/Users/joshualevendis/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff"
GENOME_FAI="/Users/joshualevendis/Documents/RNA/honours/Pfalciparum3D7/fasta/data/PlasmoDB-67_Pfalciparum3D7_Genome.fasta.fai"

# Output directory
OUTDIR="feature_mapped_bedmethyls"
mkdir -p "$OUTDIR"

# awk '{print $3}' /Users/joshualevendis/Documents/RNA/honours/Pfalciparum3D7/gff/data/PlasmoDB-67_Pfalciparum3D7.gff | sort | uniq

# Extract relevant features from GFF
awk '$3 == "five_prime_UTR" {print $0}' ${ANNOTATION} > ${OUTDIR}/pfal_5p.gff
awk '$3 == "three_prime_UTR" {print $0}' ${ANNOTATION} > ${OUTDIR}/pfal_3p.gff
awk '$3 == "CDS" {print $0}' ${ANNOTATION} > ${OUTDIR}/pfal_cds.gff
awk '$3 == "ncRNA" {print $0}' ${ANNOTATION} > ${OUTDIR}/pfal_ncRNA.gff
awk '$3 == "rRNA" {print $0}' ${ANNOTATION} > ${OUTDIR}/pfal_rRNA.gff
awk '$3 == "tRNA" {print $0}' ${ANNOTATION} > ${OUTDIR}/pfal_tRNA.gff
awk '$3 == "snRNA" {print $0}' ${ANNOTATION} > ${OUTDIR}/pfal_snRNA.gff
awk '$3 == "snoRNA" {print $0}' ${ANNOTATION} > ${OUTDIR}/pfal_snoRNA.gff

awk '$3 == "mRNA" || $3 == "ncRNA_gene" {print $0}' ${ANNOTATION} > ${OUTDIR}/pfal_transcripts.gff

bedtools slop -i ${OUTDIR}/pfal_transcripts.gff -g ${GENOME_FAI} -b 500 > ${OUTDIR}/ex500_pfal_transcripts.gff
bedtools slop -i ${OUTDIR}/pfal_transcripts.gff -g ${GENOME_FAI} -b 1000 > ${OUTDIR}/ex1000_pfal_transcripts.gff

# Get list of mRNAs that overlap with another mRNA
bedtools intersect -a ${OUTDIR}/pfal_transcripts.gff -b ${OUTDIR}/pfal_transcripts.gff -wa -wb | awk -F '\t' '$4 != $13' > ${OUTDIR}/pfal_transcripts_overlaps.gff
bedtools intersect -a ${OUTDIR}/ex500_pfal_transcripts.gff -b ${OUTDIR}/ex500_pfal_transcripts.gff -wa -wb | awk -F '\t' '$4 != $13' > ${OUTDIR}/ex500_pfal_transcripts_overlaps.gff
bedtools intersect -a ${OUTDIR}/ex1000_pfal_transcripts.gff -b ${OUTDIR}/ex1000_pfal_transcripts.gff -wa -wb | awk -F '\t' '$4 != $13' > ${OUTDIR}/ex1000_pfal_transcripts_overlaps.gff


# Loop over each file
for f in "${FILES[@]}"; do
    echo "Processing $f ..."
    BASENAME=$(basename "$f")          # 28C1_to_pfal_0.95.0p.a.bed
    PREFIX=$(echo "$BASENAME" | cut -d'_' -f1)          # 28C1
    OUTFILE="${OUTDIR}/${PREFIX}_to_pfal_0.95.0p.a.feature_mapped.bed"
    
    # Perform the intersection with each feature file
    bedtools intersect -a "$f" -b "${OUTDIR}/pfal_5p.gff" -wa -c -s | \
    bedtools intersect -a - -b "${OUTDIR}/pfal_3p.gff" -wa -c -s | \
    bedtools intersect -a - -b "${OUTDIR}/pfal_cds.gff" -wa -c -s | \
    bedtools intersect -a - -b "${OUTDIR}/pfal_ncRNA.gff" -wa -c -s | \
    bedtools intersect -a - -b "${OUTDIR}/pfal_rRNA.gff" -wa -c -s | \
    bedtools intersect -a - -b "${OUTDIR}/pfal_tRNA.gff" -wa -c -s | \
    bedtools intersect -a - -b "${OUTDIR}/pfal_snRNA.gff" -wa -c -s | \
    bedtools intersect -a - -b "${OUTDIR}/pfal_snoRNA.gff" -wa -c -s > "$OUTFILE"

    echo "Output written to $OUTFILE"

    # Keep only sites which are not in overlapping mRNAs
    # ie genes which are far enough away from their neighbours
    bedtools intersect -a "$OUTFILE" -b "${OUTDIR}/pfal_transcripts_overlaps.gff" -v > "${OUTDIR}/${PREFIX}_to_pfal_0.95.0p.a.feature_mapped.nonoverlap.bed"
    bedtools intersect -a "$OUTFILE" -b "${OUTDIR}/ex500_pfal_transcripts_overlaps.gff" -v > "${OUTDIR}/${PREFIX}_to_pfal_0.95.0p.a.feature_mapped.nonoverlap500.bed"
    bedtools intersect -a "$OUTFILE" -b "${OUTDIR}/ex1000_pfal_transcripts_overlaps.gff" -v > "${OUTDIR}/${PREFIX}_to_pfal_0.95.0p.a.feature_mapped.nonoverlap1000.bed"

done

echo "All files processed. Results in $OUTDIR/"
