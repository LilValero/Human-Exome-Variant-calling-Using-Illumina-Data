# Human Exome Variant Calling Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description <a name="description"></a>

This repository contains a robust bioinformatics pipeline designed for processing Illumina human exome sequencing data to identify genetic variants (SNPs and InDels). The workflow utilizes a suite of standard, widely-adopted bioinformatics tools within a defined Conda environment to ensure reproducibility.

This project demonstrates an end-to-end solution for NGS data analysis, emphasizing automation, quality control, and clear documentation practices suitable for research and development environments.

## Features <a name="features"></a>

✨ Key capabilities of this pipeline include:

* 🔬 **Quality Control:** Initial assessment of raw sequencing read quality using FastQC.
* ✂️ **Read Trimming:** Adapter removal and quality filtering of reads via Trimmomatic.
* 🧬 **Alignment:** Mapping processed reads to a reference genome utilizing BWA-MEM.
* ⚙️ **Post-Alignment Processing:** Efficient BAM file sorting and indexing with Samtools.
* 🎯 **Variant Calling:** Accurate identification of SNPs and short InDels using the bcftools `mpileup` and `call` workflow.
* 🏷️ **Variant Annotation:** Functional consequence prediction and annotation of identified variants using SnpEff.
* 📊 **Aggregate Reporting:** Consolidation of QC metrics across the workflow into a comprehensive report using MultiQC.
* 🚀 **Automation:** The entire pipeline is orchestrated via a single, configurable shell script (`pipeline.sh`).
* 📦 **Reproducibility:** A defined Conda environment (`environment.yml`) ensures consistent results by managing software versions and dependencies.

## Project Structure <a name="structure"></a>

The repository is organized following standard practices for clarity and maintainability:

```text
.
+-- alignment/            # Output: Sorted BAM files and indices
+-- data/                 # Input: FASTQ reads, reference genome, adapter sequences
|   +-- adapters/
|   +-- reference/
+-- docs/                 # Output: HTML reports for GitHub Pages (SnpEff, MultiQC)
+-- logs/                 # Output: Log files from pipeline execution
+-- qc/                   # Output: FastQC reports and MultiQC results directory
+-- trimmed/              # Output: Trimmed FASTQ files
+-- variants/             # Output: VCF files and SnpEff annotation directory
|   +-- snpEff_annotation/
+-- .gitignore            # Configuration: Specifies intentionally untracked files
+-- environment.yml       # Configuration: Conda environment definition
+-- LICENSE               # Documentation: Project license details
+-- pipeline.sh           # Core: Main pipeline execution script
+-- README.md             # Documentation: This file

Setup and Installation <a name="setup"></a>
Follow these steps to set up the project environment locally.
Prerequisites:
 * Git
 * Miniconda or Anaconda
 * Clone Repository:
   git clone [https://github.com/LilValero/Human-Exome-Variant-calling-Using-Illumina-Data.git](https://github.com/LilValero/Human-Exome-Variant-calling-Using-Illumina-Data.git)
cd Human-Exome-Variant-calling-Using-Illumina-Data

 * Create Conda Environment:
   This command uses the provided environment.yml file to install all necessary tools with the correct versions.
   conda env create -f environment.yml
conda activate variant_calling_env # Activate the created environment

 * Prepare Input Data:
   * Place raw paired-end FASTQ files (*.fastq.gz) in the data/ directory.
   * Download & index the full GRCh38 primary assembly (all chromosomes) on‐the‐fly:
  ```bash
  REF=Homo_sapiens.GRCh38.dna.primary_assembly.fa
  mkdir -p data
  curl -O ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/${REF}.gz
  gunzip ${REF}.gz
  bwa index data/${REF}
   * Place the Trimmomatic adapter sequence file (.fa) in data/adapters/.
   * Ensure the required SnpEff database (e.g., GRCh38.p13.RefSeq) is available to SnpEff. This may involve running snpEff download <database_name> if it's the first time using it.
Usage <a name="usage"></a>
 * Configure Pipeline:
   * Open pipeline.sh in a text editor.
   * Modify the variables within the --- Configuration --- section to match your input file names/paths (FASTQ_R1, FASTQ_R2, SAMPLE_NAME, REF_GENOME, ADAPTERS) and SnpEff settings (SNPEFF_DB, SNPEFF_CONFIG).
 * Make Executable:
   (Only needed once) Grant execute permissions to the script:
   chmod +x pipeline.sh

 * Run Pipeline:
   Execute the script from the project's root directory:
   ./pipeline.sh

   Pipeline progress and tool-specific output will be directed to log files within the logs/ directory.
Results <a name="results"></a>
Upon successful completion, the pipeline generates several key outputs:
 * 📄 SnpEff Annotation Summary: An interactive HTML report detailing variant annotations.
   * View SnpEff Report (Live Demo)
 * 📈 MultiQC Aggregate Report: A comprehensive HTML report summarizing QC metrics from various pipeline stages.
   * View MultiQC Report (Live Demo)
 * 🧬 Annotated Variants: The final annotated variant calls are available in VCF format at: variants/snpEff_annotation/<SAMPLE_NAME>.annotated.vcf.gz
Intermediate files (trimmed reads, alignments, raw VCFs) are stored in their respective directories (trimmed/, alignment/, variants/).
Contributing <a name="contributing"></a>
Contributions, issues, and feature requests are welcome. Please refer to the issues page.
License <a name="license"></a>
This project is licensed under the MIT License. See the LICENSE file for full details.