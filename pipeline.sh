#!/bin/bash
set -e -o pipefail # Exit on error, treat errors in pipes as failure

# --- Configuration ---
# Adjust these paths and settings according to your project setup

# Input/Reference Files
FASTQ_R1="data/SRRXXXXXX_1.fastq.gz" # CHANGE to your Read 1 FASTQ file
FASTQ_R2="data/SRRXXXXXX_2.fastq.gz" # CHANGE to your Read 2 FASTQ file
SAMPLE_NAME="SRRXXXXXX"             # CHANGE to your sample identifier
REF_GENOME="data/Homo_sapiens.GRCh38.dna.primary_assembly.fa" # CHANGE to your reference genome path
ADAPTERS="data/adapters/TruSeq3-PE-2.fa"      # CHANGE to your Trimmomatic adapter file path

# SnpEff Settings
SNPEFF_DB="GRCh38.p13.RefSeq" # CHANGE to your SnpEff database name (e.g., hg38, GRCh38.XYZ)
SNPEFF_CONFIG="$HOME/snpEff/snpEff.config" # CHANGE path if snpEff.config is not default

# Resources
THREADS=8 # Number of threads to use

# Tool Paths (Assumes tools are in the Conda environment's PATH)
# If not, specify full paths here, e.g., TRIMMOMATIC_JAR="/path/to/trimmomatic.jar"
TRIMMOMATIC_JAR=$(which trimmomatic.jar) # Assumes 'trimmomatic.jar' script is in PATH or provide path directly

# Output Directories (Relative to project root)
QC_DIR="qc"
TRIM_DIR="trimmed"
ALIGN_DIR="alignment"
VAR_DIR="variants"
SNPEFF_DIR="variants/snpEff_annotation"
LOG_DIR="logs"
MULTIQC_DIR="qc/multiqc_report"

# --- Pipeline ---

echo "Starting Variant Calling Pipeline for Sample: ${SAMPLE_NAME}"
date

# Create output directories if they don't exist
mkdir -p "${QC_DIR}" "${TRIM_DIR}" "${ALIGN_DIR}" "${VAR_DIR}" "${SNPEFF_DIR}" "${LOG_DIR}" "${MULTIQC_DIR}"

# --- Step 1: Initial Quality Control (FastQC) ---
echo "[$(date +%T)] Running FastQC on raw reads..."
fastqc -t "${THREADS}" -o "${QC_DIR}" "${FASTQ_R1}" "${FASTQ_R2}" &> "${LOG_DIR}/${SAMPLE_NAME}_fastqc_raw.log"
echo "[$(date +%T)] FastQC on raw reads complete."

# --- Step 2: Adapter Trimming (Trimmomatic) ---
echo "[$(date +%T)] Running Trimmomatic..."
TRIM_R1_PAIRED="${TRIM_DIR}/${SAMPLE_NAME}_R1_paired.fastq.gz"
TRIM_R1_UNPAIRED="${TRIM_DIR}/${SAMPLE_NAME}_R1_unpaired.fastq.gz"
TRIM_R2_PAIRED="${TRIM_DIR}/${SAMPLE_NAME}_R2_paired.fastq.gz"
TRIM_R2_UNPAIRED="${TRIM_DIR}/${SAMPLE_NAME}_R2_unpaired.fastq.gz"

java -jar "${TRIMMOMATIC_JAR}" PE -threads "${THREADS}" \
    "${FASTQ_R1}" "${FASTQ_R2}" \
    "${TRIM_R1_PAIRED}" "${TRIM_R1_UNPAIRED}" \
    "${TRIM_R2_PAIRED}" "${TRIM_R2_UNPAIRED}" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10:8:true" \
    "LEADING:3" "TRAILING:3" "SLIDINGWINDOW:4:15" "MINLEN:36" \
    &> "${LOG_DIR}/${SAMPLE_NAME}_trimmomatic.log"
echo "[$(date +%T)] Trimmomatic complete."

# --- Step 3: Post-Trimming Quality Control (FastQC) ---
echo "[$(date +%T)] Running FastQC on trimmed reads..."
fastqc -t "${THREADS}" -o "${QC_DIR}" "${TRIM_R1_PAIRED}" "${TRIM_R2_PAIRED}" &> "${LOG_DIR}/${SAMPLE_NAME}_fastqc_trimmed.log"
echo "[$(date +%T)] FastQC on trimmed reads complete."

# --- Step 4: Download & index full GRCh38 primary assembly ---
if [[ ! -f "$REF_GENOME" ]]; then
  mkdir -p data
  cd data
  curl -O ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  cd ..
fi
echo "[$(date +%T)] Indexing full GRCh38 with BWA…"
bwa index "$REF_GENOME"
echo "[$(date +%T)] BWA indexing complete."

# --- Step 5: Alignment (BWA-MEM) and Sorting (SAMtools) ---
echo "[$(date +%T)] Aligning reads with BWA-MEM and sorting..."
SORTED_BAM="${ALIGN_DIR}/${SAMPLE_NAME}.sorted.bam"
# Define Read Group: crucial for downstream tools like GATK, bcftools
RG_ID="group_${SAMPLE_NAME}"
RG_SM="${SAMPLE_NAME}"
RG_PL="ILLUMINA"
RG_LB="lib_${SAMPLE_NAME}"
READ_GROUP_INFO="@RG\\tID:${RG_ID}\\tSM:${RG_SM}\\tPL:${RG_PL}\\tLB:${RG_LB}"

bwa mem -t "${THREADS}" -R "${READ_GROUP_INFO}" "${REF_GENOME}" "${TRIM_R1_PAIRED}" "${TRIM_R2_PAIRED}" 2> "${LOG_DIR}/${SAMPLE_NAME}_bwa_mem.log" | \
    samtools view -@ "${THREADS}" -bhS - > "${ALIGN_DIR}/${SAMPLE_NAME}.bam"

# Sort the BAM file
samtools sort -@ "${THREADS}" -o "${SORTED_BAM}" "${ALIGN_DIR}/${SAMPLE_NAME}.bam"
rm "${ALIGN_DIR}/${SAMPLE_NAME}.bam" # Remove unsorted BAM
echo "[$(date +%T)] Alignment and sorting complete."

# --- Step 6: Indexing BAM File (SAMtools) ---
echo "[$(date +%T)] Indexing sorted BAM file..."
samtools index -@ "${THREADS}" "${SORTED_BAM}" &> "${LOG_DIR}/${SAMPLE_NAME}_samtools_index.log"
echo "[$(date +%T)] BAM indexing complete."

# --- Step 7: Variant Calling (bcftools) ---
echo "[$(date +%T)] Calling variants with bcftools..."
RAW_VCF="${VAR_DIR}/${SAMPLE_NAME}.raw.vcf.gz"
FINAL_VCF="${VAR_DIR}/${SAMPLE_NAME}.final.vcf.gz" # Using final name here, filtering is optional step

bcftools mpileup -Ou -f "${REF_GENOME}" "${SORTED_BAM}" --threads "${THREADS}" 2> "${LOG_DIR}/${SAMPLE_NAME}_mpileup.log" | \
    bcftools call -mv -Ov --threads "${THREADS}" -o "${VAR_DIR}/${SAMPLE_NAME}.call.vcf" - 2> "${LOG_DIR}/${SAMPLE_NAME}_bcftools_call.log"

# Optional: Add filtering step here using bcftools filter if desired
# Example: bcftools filter -Oz -o "${FINAL_VCF}" -i 'QUAL>30 && DP>10' "${VAR_DIR}/${SAMPLE_NAME}.call.vcf"
# If filtering, use RAW_VCF for the output of 'call' and FINAL_VCF for the output of 'filter'
# For now, we just compress and index the called VCF
bgzip -c "${VAR_DIR}/${SAMPLE_NAME}.call.vcf" > "${FINAL_VCF}"
rm "${VAR_DIR}/${SAMPLE_NAME}.call.vcf" # Remove uncompressed VCF
tabix -p vcf "${FINAL_VCF}" &> "${LOG_DIR}/${SAMPLE_NAME}_tabix.log"

echo "[$(date +%T)] Variant calling complete: ${FINAL_VCF}"

# --- Step 8: Variant Annotation (SnpEff) ---
echo "[$(date +%T)] Annotating variants with SnpEff..."
ANNOTATED_VCF="${SNPEFF_DIR}/${SAMPLE_NAME}.annotated.vcf.gz"
SNPEFF_HTML_REPORT="${SNPEFF_DIR}/snpeff_summary.html" # Note: Path inside variants/snpEff_annotation

# Ensure the SnpEff output directory exists
mkdir -p "$(dirname ${ANNOTATED_VCF})"

snpEff -Xmx8g -v "${SNPEFF_DB}" "${FINAL_VCF}" \
    -stats "${SNPEFF_HTML_REPORT}" \
    -config "${SNPEFF_CONFIG}" \
    -gz \
    > "$(echo ${ANNOTATED_VCF} | sed 's/\.gz$//')" 2> "${LOG_DIR}/${SAMPLE_NAME}_snpEff.log"

# SnpEff might not output gzipped directly depending on version/config, gzip if needed
if [[ -f "$(echo ${ANNOTATED_VCF} | sed 's/\.gz$//')" && ! -f "${ANNOTATED_VCF}" ]]; then
    bgzip "$(echo ${ANNOTATED_VCF} | sed 's/\.gz$//')"
fi

# Move the HTML report to the docs folder for GitHub Pages
# Adjust source path if needed
cp "${SNPEFF_HTML_REPORT}" "docs/snpeff_summary.html"
echo "[$(date +%T)] SnpEff annotation complete. Report copied to docs/."

# --- Step 9: Aggregate Quality Control (MultiQC) ---
echo "[$(date +%T)] Running MultiQC to aggregate logs..."
multiqc . -f -o "${MULTIQC_DIR}" -n multiqc_report &> "${LOG_DIR}/multiqc.log"
# Copy MultiQC report to docs as well if desired
cp "${MULTIQC_DIR}/multiqc_report.html" "docs/multiqc_report.html"
echo "[$(date +%T)] MultiQC complete. Report copied to docs/."

# --- Pipeline End ---
echo "[$(date +%T)] Variant Calling Pipeline Finished Successfully!"
date
