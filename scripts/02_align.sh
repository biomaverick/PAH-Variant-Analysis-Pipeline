#!/usr/bin/env bash
# 02_align.sh — align trimmed reads to reference, sort, mark duplicates
# Usage: ./02_align.sh <sample_name> <reference.fa>
# Expects trimmed reads at results/trimmed/<sample>_R1.trimmed.fastq.gz etc.
set -euo pipefail

SAMPLE=$1
REF=$2

mkdir -p results/aligned

# Index reference once (skip if already indexed)
if [ ! -f "${REF}.bwt" ]; then
  echo ">> Indexing reference $REF"
  bwa index "$REF"
  samtools faidx "$REF"
fi

echo ">> Aligning $SAMPLE with BWA-MEM"
bwa mem -t 4 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    "$REF" \
    "results/trimmed/${SAMPLE}_R1.trimmed.fastq.gz" \
    "results/trimmed/${SAMPLE}_R2.trimmed.fastq.gz" \
  | samtools sort -@ 4 -o "results/aligned/${SAMPLE}.sorted.bam"

samtools index "results/aligned/${SAMPLE}.sorted.bam"

echo ">> Marking duplicates for $SAMPLE"
gatk MarkDuplicates \
  -I "results/aligned/${SAMPLE}.sorted.bam" \
  -O "results/aligned/${SAMPLE}.dedup.bam" \
  -M "results/aligned/${SAMPLE}.dup_metrics.txt"

samtools index "results/aligned/${SAMPLE}.dedup.bam"

echo ">> Done. Final BAM: results/aligned/${SAMPLE}.dedup.bam"
echo ">> Sanity check: open this BAM in IGV at the PAH locus (chr12:102,836,882-102,958,435, GRCh38)"
