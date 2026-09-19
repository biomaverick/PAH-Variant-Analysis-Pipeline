#!/usr/bin/env bash
# 01_qc.sh — raw read QC + adapter trimming
# Usage: ./01_qc.sh <sample_name> <R1.fastq.gz> <R2.fastq.gz>
set -euo pipefail

SAMPLE=$1
R1=$2
R2=$3

mkdir -p results/fastqc results/trimmed results/fastp_reports

echo ">> Running FastQC on raw reads for $SAMPLE"
fastqc "$R1" "$R2" -o results/fastqc/

echo ">> Trimming adapters with fastp for $SAMPLE"
fastp \
  -i "$R1" -I "$R2" \
  -o "results/trimmed/${SAMPLE}_R1.trimmed.fastq.gz" \
  -O "results/trimmed/${SAMPLE}_R2.trimmed.fastq.gz" \
  --detect_adapter_for_pe \
  --json "results/fastp_reports/${SAMPLE}.fastp.json" \
  --html "results/fastp_reports/${SAMPLE}.fastp.html"

echo ">> Done. Trimmed reads in results/trimmed/"
echo ">> Run 'multiqc results/' after processing all samples for a combined report."
