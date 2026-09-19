#!/usr/bin/env bash
# 03_call_variants.sh — per-sample GVCF calling, joint genotyping, hard filtering
# Usage: ./03_call_variants.sh <reference.fa> <region> <sample1> [sample2 ...]
# region example: chr12:102836882-102958435
set -euo pipefail

REF=$1
REGION=$2
shift 2
SAMPLES=("$@")

mkdir -p results/variants

GVCF_ARGS=()
for SAMPLE in "${SAMPLES[@]}"; do
  echo ">> HaplotypeCaller (GVCF mode) for $SAMPLE"
  gatk HaplotypeCaller \
    -R "$REF" \
    -I "results/aligned/${SAMPLE}.dedup.bam" \
    -O "results/variants/${SAMPLE}.g.vcf.gz" \
    -ERC GVCF \
    -L "$REGION"
  GVCF_ARGS+=(-V "results/variants/${SAMPLE}.g.vcf.gz")
done

echo ">> Combining GVCFs into GenomicsDB"
rm -rf results/variants/pah_db
gatk GenomicsDBImport \
  "${GVCF_ARGS[@]}" \
  --genomicsdb-workspace-path results/variants/pah_db \
  -L "$REGION"

echo ">> Joint genotyping"
gatk GenotypeGVCFs \
  -R "$REF" \
  -V gendb://results/variants/pah_db \
  -O results/variants/pah_cohort.vcf.gz

echo ">> Hard-filtering (GATK best-practice thresholds for small cohorts)"
gatk VariantFiltration \
  -R "$REF" \
  -V results/variants/pah_cohort.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  -O results/variants/pah_cohort.filtered.vcf.gz

echo ">> Done. Filtered VCF: results/variants/pah_cohort.filtered.vcf.gz"
echo ">> If using NA12878, validate against the GIAB truth VCF now with bcftools isec or hap.py"
