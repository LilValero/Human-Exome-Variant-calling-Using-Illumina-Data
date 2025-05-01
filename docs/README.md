# Human Exome Variant Calling Pipeline (Illumina, Single-End)

This project demonstrates a full variant calling pipeline for human exome data, executed in a local Conda environment using Illumina FASTQ files from NCBI SRA. It includes read QC, alignment, variant calling, and annotation using well-established bioinformatics tools.

> **Project Goal**: Detect and annotate functional variants in a human exome dataset (chr22) and visualize high-impact mutations.

---

## 🔬 Summary for Recruiters & Interviewers

This project was developed to showcase bioinformatics skills in real-world genomic data analysis. It reflects fluency with standard pipelines, data hygiene, reproducibility practices, and interpretation of high-throughput sequencing results. It is fully documented and reproducible on macOS.

- 🧪 **Tools**: FastQC, Trimmomatic, BWA, SAMtools, BCFtools, SnpEff, MultiQC
- 🧬 **Pipeline**: Illumina FASTQ → QC → Trimming → Alignment → Variant Calling → Annotation → Visualization
- 💻 **Stack**: Bash, Conda, macOS, HTML, GitHub Pages
- 🧠 **Skills**: NGS analysis, reproducible workflows, VCF interpretation, shell scripting, documentation

---

## 🚀 Pipeline Overview

### 1. 📥 Data Acquisition
- Dataset: [SRR2138889](https://www.ncbi.nlm.nih.gov/sra/SRR2138889)
- Source: NCBI SRA
- Platform: Illumina, single-end

### 2. ✅ Quality Control
- Tool: FastQC + MultiQC
- Outcome: Per-base quality, overrepresented sequences, adapter contamination

### 3. ✂️ Trimming
- Tool: Trimmomatic
- Settings: Adapter removal + quality filtering (`SLIDINGWINDOW`, `MINLEN`)

### 4. 🧬 Reference Genome
- Source: Ensembl GRCh38 release 110 (chr22)
- Indexing: `bwa index`

### 5. 🧷 Alignment
- Tool: BWA MEM
- Output: `aligned.sam`, converted/sorted/indexed to BAM

### 6. 🔎 Variant Calling
- Tool: `bcftools mpileup` + `bcftools call`
- Output: `variants.vcf` (8,719 raw variants)

### 7. 🧠 Annotation
- Tool: SnpEff with GRCh38.86 database
- Output: `annotated_variants.vcf`
- Summary: `snpeff_summary.html` (HTML report)

### 8. 🎯 High-Impact Extraction
- `grep "HIGH"` to isolate functionally significant variants → `high_impact_variants.vcf`

---

## 📂 Repository Structure

```
exome_project/
├── data/                          # Raw FASTQ and reference
├── qc/                            # FastQC + MultiQC outputs
├── results/                       # BAM, VCF, snpEff output
├── scripts/                       # (optional) pipeline.sh
├── snpeff_summary.html            # HTML annotation summary
├── README.md                      # This file
```

---

## 📈 GitHub Pages Report

[Click to View Variant Annotation Report](https://yourusername.github.io/exome_variant_pipeline/snpeff_summary.html)

Includes:
- Distribution of variant impacts
- Top affected genes
- Effect categories (e.g., missense, nonsense)

---

## 🧪 Technologies Used
- `conda` for reproducible environment management
- `bash` for scripting and pipeline chaining
- `bioconda` for easy installation of genomics tools
- `GitHub Pages` for hosting visual outputs

---

## 🧬 Keywords / Tags (for discoverability)
```
bioinformatics, variant calling, snpeff, exome, FASTQ, VCF, bwa, samtools, bcftools, Illumina, human genome, GRCh38, next-generation sequencing, cancer genomics, data science, computational biology
```

---

## ✅ Resume-Friendly Project Line
> Built a complete Illumina exome variant-calling pipeline from raw FASTQ to annotated VCF using BWA, SAMtools, bcftools, and snpEff. Identified 8,700+ variants, extracted high-impact mutations, and generated a public HTML annotation summary using GitHub Pages.


