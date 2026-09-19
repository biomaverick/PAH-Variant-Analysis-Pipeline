#!/usr/bin/env bash
# 04_annotate.sh — VEP annotation of filtered variants
# Usage: ./04_annotate.sh
# Requires VEP cache installed for GRCh38 (run once: vep_install -a cf -s homo_sapiens -y GRCh38)
set -euo pipefail

mkdir -p results/annotation

echo ">> Running VEP"
vep --input_file results/variants/pah_cohort.filtered.vcf.gz \
    --output_file results/annotation/pah_annotated.vcf \
    --vcf \
    --cache --offline --assembly GRCh38 \
    --sift b --polyphen b \
    --force_overwrite

echo ">> Done. Annotated file: results/annotation/pah_annotated.vcf"
echo ">> Next: filter to missense/nonsense variants in PAH, cross-reference ClinVar and PAHvdb by hand"
echo ">> (PAHvdb: http://www.biopku.org/pah/ — search each candidate variant's HGVS notation)"
