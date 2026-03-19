#!/bin/bash

# =================================================================
# GATK Variant Calling Pipeline (Optimized for GitHub)
# =================================================================

# --- Parameter Definitions (Adjust relative paths as needed) ---
# It is recommended to place reference and raw data in the 'data' directory
BASE_DIR=$(pwd)
REFERENCE="${BASE_DIR}/reference/genomic.fna"
SAMPLES_DIR="${BASE_DIR}/raw_reads/"
OUTPUT_DIR="${BASE_DIR}/results"
THREADS=48

# --- Create Output Directory Structure ---
echo "Initializing directory structure..."
DIRS=(
    "aligned_reads" "dedup_sorted" "metrics" "Raw_variants_Round1"
    "Raw_snps_Round1" "Raw_indels_Round1" "Filtered_snps_round1" "Filtered_indels_round1"
    "bqsr_snps" "bqsr_indels" "recal_data" "recal_reads" "recalibration_plots"
    "post_recal_data" "raw_variants_recal" "Raw_snps_recal" "Raw_indels_recal"
    "filtered_snps_final" "filtered_indels_final" "filtered_snps_final_ann"
)

for dir in "${DIRS[@]}"; do
    mkdir -p "${OUTPUT_DIR}/${dir}"
done

# --- Check Reference Genome Index ---
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "Reference index not found. Generating index..."
    bwa index "$REFERENCE"
fi

# --- Batch Processing Samples ---
# Assuming file suffix is _R1.trim.fq.gz
for R1_FILE in ${SAMPLES_DIR}/*_R1.trim.fq.gz; do
    # Extract sample name (get filename from path and remove suffix)
    FILENAME=$(basename "$R1_FILE")
    SAMPLE="${FILENAME%_R1.trim.fq.gz}" 
    
    R2_FILE="${SAMPLES_DIR}/${SAMPLE}_R2.trim.fq.gz"
    BASENAME="$SAMPLE"
    
    # Define output file paths for each step
    SAM_FILE="${OUTPUT_DIR}/aligned_reads/${BASENAME}.sam"
    DEDUP_BAM_FILE="${OUTPUT_DIR}/dedup_sorted/${BASENAME}_sorted_dedup.bam"
    METRICS_FILE="${OUTPUT_DIR}/metrics/${BASENAME}_dedup_metrics.txt"
    ALIGNMENT_METRICS="${OUTPUT_DIR}/metrics/${BASENAME}_alignment_metrics.txt"
    INSERT_METRICS="${OUTPUT_DIR}/metrics/${BASENAME}_insert_metrics.txt"
    INSERT_HISTOGRAM="${OUTPUT_DIR}/metrics/${BASENAME}_insert_size_histogram.pdf"
    DEPTH_FILE="${OUTPUT_DIR}/metrics/${BASENAME}_depth_out.txt"
    RAW_VARIANTS_Round1="${OUTPUT_DIR}/Raw_variants_Round1/${BASENAME}_raw_variants_Round1.vcf"
    RAW_SNPS_Round1="${OUTPUT_DIR}/Raw_snps_Round1/${BASENAME}_raw_snps_Round1.vcf"
    RAW_INDELS_Round1="${OUTPUT_DIR}/Raw_indels_Round1/${BASENAME}_raw_indels_Round1.vcf"
    Filtered_SNPS_Round1="${OUTPUT_DIR}/Filtered_snps_round1/${BASENAME}_filtered_snps_Round1.vcf"
    Filtered_INDELS_Round1="${OUTPUT_DIR}/Filtered_indels_round1/${BASENAME}_filtered_indels_Round1.vcf"
    BQSR_SNPS="${OUTPUT_DIR}/bqsr_snps/${BASENAME}_bqsr_snps.vcf"
    BQSR_INDELS="${OUTPUT_DIR}/bqsr_indels/${BASENAME}_bqsr_indels.vcf"
    Recal_data="${OUTPUT_DIR}/recal_data/${BASENAME}_recal_data.table"
    Recal_reads="${OUTPUT_DIR}/recal_reads/${BASENAME}_recal_reads.bam"
    Post_recal_data="${OUTPUT_DIR}/post_recal_data/${BASENAME}_post_recal_data.table"
    Recalibration_plots="${OUTPUT_DIR}/recalibration_plots/${BASENAME}_recalibration_plots.pdf"
    Raw_variants_recal="${OUTPUT_DIR}/raw_variants_recal/${BASENAME}_raw_variants_recal.vcf"
    RAW_SNPS_recal="${OUTPUT_DIR}/Raw_snps_recal/${BASENAME}_raw_snps_recal.vcf"
    RAW_INDELS_recal="${OUTPUT_DIR}/Raw_indels_recal/${BASENAME}_raw_indels_recal.vcf"
    Filtered_SNPS_final="${OUTPUT_DIR}/filtered_snps_final/${BASENAME}_filtered_snps_final.vcf"
    Filtered_INDELS_final="${OUTPUT_DIR}/filtered_indels_final/${BASENAME}_filtered_indels_final.vcf"
    
    echo "----------------------------------------------------"
    echo "Processing Sample: $BASENAME"
    echo "----------------------------------------------------"
    
    # Step 1: BWA Alignment
    echo "Step 1: Running BWA MEM..."
    bwa mem -K 100000000 -Y \
        -R "@RG\tID:${BASENAME}\tLB:${BASENAME}\tPL:ILLUMINA\tPM:HISEQ\tSM:${BASENAME}" \
        -t $THREADS "$REFERENCE" "$R1_FILE" "$R2_FILE" > "$SAM_FILE"
    
    # Step 2: Mark Duplicates
    echo "Step 2: Running MarkDuplicatesSpark..."
    gatk MarkDuplicatesSpark -I "$SAM_FILE" -M "$METRICS_FILE" -O "$DEDUP_BAM_FILE"  
    
    # Step 3 & 4: Picard Metrics
    echo "Step 3-4: Collecting Quality Metrics (Picard)..."
    picard CollectAlignmentSummaryMetrics R="$REFERENCE" I="$DEDUP_BAM_FILE" O="$ALIGNMENT_METRICS"
    picard CollectInsertSizeMetrics INPUT="$DEDUP_BAM_FILE" OUTPUT="$INSERT_METRICS" HISTOGRAM_FILE="$INSERT_HISTOGRAM" 
    
    # Step 5: Samtools depth
    echo "Step 5: Calculating Coverage Depth (Samtools)..."
    samtools depth -a "$DEDUP_BAM_FILE" > "$DEPTH_FILE"
    
    # Step 6: Call Variants Round 1
    echo "Step 6: Running HaplotypeCaller Round 1..."
    gatk HaplotypeCaller --sample-ploidy 50 -R "$REFERENCE" -I "$DEDUP_BAM_FILE" -O "$RAW_VARIANTS_Round1"
    
    # Step 7 & 8: Extract SNPs & INDELs
    echo "Step 7-8: Extracting SNPs and INDELs..."
    gatk SelectVariants -R "$REFERENCE" -V "$RAW_VARIANTS_Round1" -select-type SNP -O "$RAW_SNPS_Round1"
    gatk SelectVariants -R "$REFERENCE" -V "$RAW_VARIANTS_Round1" -select-type INDEL -O "$RAW_INDELS_Round1"
       
    # Step 9 & 10: Variant Filtration
    echo "Step 9-10: Applying Variant Hard Filtration..."
    gatk VariantFiltration -R "$REFERENCE" -V "$RAW_SNPS_Round1" -O "$Filtered_SNPS_Round1" \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

    gatk VariantFiltration -R "$REFERENCE" -V "$RAW_INDELS_Round1" -O "$Filtered_INDELS_Round1" \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0"
    
    # Step 11: Exclude Filtered Variants
    echo "Step 11: Excluding filtered variants for BQSR..."
    gatk SelectVariants --exclude-filtered -V "$Filtered_SNPS_Round1" -O "$BQSR_SNPS"
    gatk SelectVariants --exclude-filtered -V "$Filtered_INDELS_Round1" -O "$BQSR_INDELS"
    
    # Step 12-15: BQSR Workflow
    echo "Step 12-15: Performing Base Quality Score Recalibration (BQSR)..."
    gatk BaseRecalibrator -R "$REFERENCE" -I "$DEDUP_BAM_FILE" --known-sites "$BQSR_SNPS" --known-sites "$BQSR_INDELS" -O "$Recal_data"
    gatk ApplyBQSR -R "$REFERENCE" -I "$DEDUP_BAM_FILE" -bqsr "$Recal_data" -O "$Recal_reads"
    gatk BaseRecalibrator -R "$REFERENCE" -I "$Recal_reads" --known-sites "$BQSR_SNPS" --known-sites "$BQSR_INDELS" -O "$Post_recal_data"
    gatk AnalyzeCovariates -before "$Recal_data" -after "$Post_recal_data" -plots "$Recalibration_plots"
     
    # Step 16-19: Final Call & Filter
    echo "Step 16-19: Final Variant Calling and Filtration..."
    gatk HaplotypeCaller --sample-ploidy 50 -R "$REFERENCE" -I "$Recal_reads" -O "$Raw_variants_recal"
    
    gatk SelectVariants -R "$REFERENCE" -V "$Raw_variants_recal" -select-type SNP -O "$RAW_SNPS_recal"
    gatk SelectVariants -R "$REFERENCE" -V "$Raw_variants_recal" -select-type INDEL -O "$RAW_INDELS_recal"
    
    gatk VariantFiltration -R "$REFERENCE" -V "$RAW_SNPS_recal" -O "$Filtered_SNPS_final" \
        -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

    gatk VariantFiltration -R "$REFERENCE" -V "$RAW_INDELS_recal" -O "$Filtered_INDELS_final" \
        -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"
    
    # Step 20: Annotation
    echo "Step 20: Functional Annotation with snpEff..."
    snpEff -v GCF_000146045.2 "$Filtered_SNPS_final" \
        -s "${OUTPUT_DIR}/filtered_snps_final_ann/${BASENAME}_snpEff_summary.html" \
        > "${OUTPUT_DIR}/filtered_snps_final_ann/${BASENAME}_filtered_snps_final.ann.vcf"
    
    # Cleanup Intermediate Files
    echo "Cleaning up intermediate file: $SAM_FILE"
    rm "$SAM_FILE"

done

echo "All samples processed successfully! Results stored in: $OUTPUT_DIR"