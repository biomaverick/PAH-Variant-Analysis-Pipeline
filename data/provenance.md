# Data Provenance

## Reference genome
- Build: hg19 / GRCh37
- Source: UCSC Genome Browser
- URL: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
- Chromosome naming convention: `chr`-prefixed (e.g. `chr12`)
- PAH coordinates (hg19/GRCh37): chr12:103,232,104-103,311,244

## Sample sequencing data

| Sample ID | Source | Dataset | Platform | Files | Notes |
|---|---|---|---|---|---|
| NA12878 (HG001) | GIAB | Garvan_NA12878_HG001_HiSeq_Exome | Illumina HiSeq2500, Nextera Rapid Capture Exome | NIST7035, lanes L001+L002, R1+R2 (concatenated) | Raw, unaligned FASTQ — real "off the sequencer" starting point, not pre-aligned |

**Source directory:** https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/

**Files used:**
- `NIST7035_TAAGGCGA_L001_R1_001.fastq.gz`
- `NIST7035_TAAGGCGA_L001_R2_001.fastq.gz`
- `NIST7035_TAAGGCGA_L002_R1_001.fastq.gz`
- `NIST7035_TAAGGCGA_L002_R2_001.fastq.gz`

Lanes were concatenated per read direction (same library sequenced across two physical flow-cell lanes) into `NIST7035_R1.fastq.gz` / `NIST7035_R2.fastq.gz` (40,796,412 read pairs, matched exactly between R1/R2).

**Second replicate available but not used in this pilot:** NIST7086 (same source directory) — flagged as a future-work reproducibility check.

## Exome capture regions
- File: `nexterarapidcapture_expandedexome_targetedregions.bed`
- Source: same Garvan_NA12878_HG001_HiSeq_Exome directory as above
- 201,071 targeted regions genome-wide; 13 regions confirmed overlapping the PAH locus
- Format: 4-column BED (chrom, start, end, gene-region label)

## Truth/validation data (GIAB high-confidence calls)
- Source: NIST/GIAB, NISTv3.3.2 release, GRCh37
- VCF: `HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz`
- URL: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv3.3.2/GRCh37/
- Total variant records (whole file): 3,775,119
- **Chromosome naming note:** this file uses plain numeric naming (`12`, not `chr12`) — required renaming via `bcftools annotate --rename-chrs` to match the hg19.fa/GATK-called VCF convention before comparison

## Annotation
- Tool: Ensembl VEP, GRCh37 assembly
- Run via the **web interface** (https://www.ensembl.org/Tools/VEP), not command-line, due to local Perl/DBI dependency conflicts (see main README's "Known issues" section)
- Annotations included: consequence type, canonical transcript flag, AlphaMissense classification/score, ClinPred score

## Pipeline results summary

- **Alignment:** 99.97% mapped, 99.41% properly paired, 7.67% duplication rate (BWA-MEM + GATK MarkDuplicates)
- **Depth:** ~14.3x mean (chr1 sample region), ~9.0x mean at PAH locus specifically
- **Variants called (PAH region):** 11, all PASS filter
- **Validation vs. GIAB truth (region-restricted to PAH):** TP=8, FP=3, FN=47 → Precision=72.7%, Recall=14.5%
- **Missense variants:** 1 — chr12:103,234,252 T>C, p.Tyr414Cys (Y414C), heterozygous (0/1), rs5030860
  - ClinGen PAH Expert Panel classification: pathogenic (mild PKU-associated)
  - AlphaMissense: 0.8118 (pathogenic class); ClinPred: 0.9387
  - gnomAD v2.1.1 frequency: ~0.00037
