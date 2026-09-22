# PAH Variant Discovery & Structural Impact Pipeline

## Overview

An end-to-end NGS variant-calling pipeline built on real, unaligned exome sequencing data (NA12878/HG001, NIST7035 replicate), targeting the *PAH* (phenylalanine hydroxylase) locus. The pipeline goes from raw FASTQ through alignment, variant calling, validation against the GIAB gold-standard truth set, and functional annotation. It identified a well-documented, clinically classified pathogenic *PAH* variant (Y414C) in a heterozygous state — consistent with the sample's known healthy carrier phenotype — and cross-validated the finding using two independent computational pathogenicity predictors (AlphaMissense, ClinPred) against expert clinical curation (ClinGen).

This project was built to extend a structural/computational biology background into NGS variant-calling — connecting a called variant back to its known clinical and functional significance rather than stopping at "variant called."

## Scope

- **Sample:** NA12878 (HG001), replicate NIST7035, Garvan Institute HiSeq2500 exome dataset (GIAB)
- **Region analyzed:** *PAH* locus, chr12:103,232,104-103,311,244 (hg19/GRCh37)
- **Reference build:** hg19/GRCh37 (chosen to match the exome capture kit design and avoid a liftover step)
- **Sequencing depth:** ~14.3x mean genome-wide (sampled region), ~9.0x mean at the PAH locus specifically — both consistent with this being an older (2013-era) benchmark exome dataset, well below modern clinical exome standards (100x+)
- **This is a targeted, single-gene, single-sample pilot analysis** — not a full WGS or exome-wide variant survey, and not intended to be. Scope was deliberately kept narrow to allow full validation and interpretation of every call.

## Pipeline

```
Raw FASTQ (NIST7035, 2 lanes, concatenated)
   │  01_qc.sh — FastQC, fastp (adapter trimming, quality filtering)
   ▼
Trimmed reads (38.5M read pairs retained, 94.5% pass rate)
   │  02_align.sh — BWA-MEM alignment, MarkDuplicates
   ▼
Deduplicated BAM (99.97% mapped, 99.41% properly paired, 7.67% duplication rate)
   │  03_call_variants.sh — GATK4 HaplotypeCaller → GenotypeGVCFs → hard filtering
   ▼
Filtered VCF (11 variants in PAH region, all PASS)
   │  Validation — bcftools isec vs. GIAB HG001 GRCh37 truth VCF (region-restricted)
   ▼
Precision/recall against gold-standard truth set
   │  04_annotate.sh — Ensembl VEP (GRCh37), AlphaMissense + ClinPred pathogenicity scores
   ▼
Annotated variants, consequence classification, pathogenicity prediction
   │  Manual literature/database cross-reference (ClinVar, LOVD/PAH database)
   ▼
Clinical interpretation of the one missense variant identified (Y414C)
```

## Key results

**Variant calling:** 11 variants called in the PAH region, all passing GATK hard filters. Of these, only 1 was a missense (protein-altering) variant — the remaining calls were synonymous (2) or intronic (8), consistent with the small target region and modest sequencing depth.

**Validation against GIAB truth set** (region-restricted to PAH, hg19):

| Metric | Value |
|---|---|
| True positives | 8 |
| False positives | 3 |
| False negatives (missed) | 47 |
| **Precision** | **72.7%** |
| **Recall** | **14.5%** |

Recall is low, and that is expected and explainable: at ~9x mean depth over PAH, HaplotypeCaller lacks sufficient read support to confidently call many true variants, particularly heterozygous ones, and conservatively drops them rather than over-calling. Precision is more robust at this depth, since sites that do get called tend to have reasonably strong support. This depth/recall relationship is a legitimate and expected outcome of using a lower-coverage historical benchmark dataset, not a pipeline defect — a production pipeline targeting clinical-grade sensitivity would require 50-100x+ depth.

**Headline finding — Y414C (rs5030860, c.1241A>G):**

The single missense variant identified, chr12:103,234,252 T>C (p.Tyr414Cys), is a well-documented, previously characterized *PAH* variant:
- Classified **pathogenic for phenylketonuria** by the ClinGen PAH Variant Curation Expert Panel (ACMG/AMP criteria: PM5, PP3, PS3, PM3_VeryStrong, PP4_Moderate)
- Independently classified pathogenic by multiple clinical diagnostic laboratories (GeneDx, Eurofins)
- Described in the literature as **the most common mild-PKU-associated mutation in Europe**, retaining significant residual enzyme activity and typically producing mild PKU/hyperphenylalaninemia phenotypes rather than classic severe PKU
- gnomAD population frequency ≈ 0.00037 (rare, consistent with heterozygous carrier frequency)

Computational pathogenicity predictors from this pipeline's own annotation step were concordant with the clinical classification: **AlphaMissense** score 0.81 (pathogenic class) and **ClinPred** score 0.94.

**Structural stability prediction (DynaMut2, PDB 6HYC):** four independent structure-based methods all predict Y414C to be destabilizing:

| Method | ΔΔG (kcal/mol) | Verdict |
|---|---|---|
| DynaMut2 (consensus) | -1.627 | Destabilizing |
| mCSM | -1.817 | Destabilizing |
| DUET | -1.734 | Destabilizing |
| SDM | -0.350 | Destabilizing |
| NMA/ENCoM | -0.559 | Destabilizing |

Vibrational entropy analysis (ΔΔSvib = +0.698 kcal·mol⁻¹·K⁻¹) shows a localized increase in flexibility concentrated in the C-terminal helical region surrounding the mutation site, and interatomic contact mapping shows the wild-type Tyr414 forms hydrogen-bonding interactions with neighboring residues that are lost in the Cys414 mutant — consistent with the destabilization signal and providing a concrete structural mechanism (loss of an aromatic, H-bond-capable side chain) rather than just a statistical prediction. See `results/structural/dynamut2/` for full output (normal-mode eigenvector files, PyMOL session comparing wild-type/mutant, per-residue vibrational entropy differences).

This gives Y414C four independent, converging lines of evidence: clinical curation (ClinGen), clinical lab concordance (GeneDx, Eurofins), sequence-based ML prediction (AlphaMissense, ClinPred), and structure-based biophysical prediction (DynaMut2/mCSM/DUET/SDM) — all pointing to the same pathogenic, destabilizing classification.

**Resolving the apparent contradiction:** NA12878 is a healthy reference individual with no PKU/HPA phenotype. This is fully consistent with the finding — the variant was called in the **heterozygous state (0/1)**, and PAH deficiency is autosomal recessive. A single copy of even a well-established pathogenic PAH allele does not cause disease. This distinction (checking zygosity before interpreting a pathogenicity score) was confirmed directly from the VCF genotype field rather than assumed.

## Reproduce

```bash
conda env create -f environment.yml
conda activate ngs-pah

./scripts/01_qc.sh NIST7035 data/raw/NIST7035_R1.fastq.gz data/raw/NIST7035_R2.fastq.gz
./scripts/02_align.sh NIST7035 hg19.fa
./scripts/03_call_variants.sh hg19.fa chr12:103232104-103311244 NIST7035
# Annotation was run via the Ensembl VEP web tool (GRCh37) rather than command-line,
# due to local Perl/DBI dependency conflicts on this system — see Known issues below.
```

## Data sources
See [`data/provenance.md`](data/provenance.md) for full accession details, download commands, and dataset notes.

## Known issues / environment notes

Several non-trivial environment issues came up building this on a Fedora system with conda-managed bioinformatics tools — documented here in case they help anyone reproducing this pipeline on a similar setup:

- **FastQC/Perl `libnsl.so.1` not found:** conda's Perl build expects the legacy `libnsl.so.1`, but conda-forge's `libnsl` package ships `.so.3`. Fixed with a symlink: `ln -s $CONDA_PREFIX/lib/libnsl.so.3.0.0 $CONDA_PREFIX/lib/libnsl.so.1`
- **Perl `libcrypt.so.1` GLIBC version mismatch on Fedora:** Fedora dropped the legacy `libcrypt` ABI in favor of `libxcrypt`. A plain symlink fails here (symbol versioning mismatch, not just a missing file). Fixed with Fedora's official compatibility package: `sudo dnf install libxcrypt-compat`
- **VEP command-line `--database` mode fails** with a missing `DBI.pm` error even after installing `perl-dbi`/`perl-dbd-mysql` via conda, due to a Perl version mismatch with this build of ensembl-vep (v88, from 2017). Worked around by using the **Ensembl VEP web tool** instead of the command-line client for annotation.
- **VEP requires exact chromosome-naming convention matching (`chr12` vs `12`)** between input VCF and the reference/annotation source — a silent zero-result failure mode, not an error message. Same issue arose comparing the GATK-called VCF (`chr12`) against the GIAB truth VCF (`12`); fixed both cases with `bcftools annotate --rename-chrs`.
- **VEP offline cache download (~1GB+) was skipped** in favor of the web tool, since annotating only 11 variants did not justify the download for this pilot scope.
- **FoldX/PyMOL local installation was skipped** in favor of the **DynaMut2 web server** (https://biosig.lab.uq.edu.au/dynamut2/) for structural stability prediction — avoided licensing/build overhead for local FoldX while still obtaining consensus ΔΔG predictions across four independent structure-based methods (DynaMut2, mCSM, DUET, SDM) plus normal-mode/vibrational entropy analysis. Submission required the wild-type residue, position, mutant residue, and chain matched exactly to the PDB's native numbering (verified beforehand via `grep`/`awk` on the PDB file directly, since a mismatch here silently fails).

## Tools used
FastQC, fastp, BWA-MEM, GATK4 (HaplotypeCaller, GenomicsDBImport, GenotypeGVCFs, VariantFiltration), bcftools, samtools, Ensembl VEP (AlphaMissense, ClinPred annotations), DynaMut2 (mCSM, DUET, SDM, ENCoM/NMA)

## Future work
- Extend to the second NA12878 replicate (NIST7086) as a reproducibility check
- Re-run at higher-depth exome/WGS data to assess recall improvement
- Wrap the pipeline in Snakemake/Nextflow for full workflow reproducibility
- Extend structural modeling (DynaMut2/FoldX) to any additional missense variants identified in future runs
