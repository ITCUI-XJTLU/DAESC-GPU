#!/bin/bash

#SAMPLE_IDS=("717_718" "718_719" "720_721" "721_722" "722_723" "736_737" "737_738" "738_739" "739_740" "740_741" "743_744" "968_969" "970_971" "971_972") # sample5
SAMPLE_IDS=("970_971" "971_972")

# Define chromosomes
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
# CHROMOSOMES=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

# Define reference and output directories
REFERENCE_DIR="/projects/gqilab/DAESC_GPU/data/reference/refdata-gex-GRCh38-2020-A"
GATK_DIR="/projects/gqilab/DAESC_GPU/data/reference/gatk"
PHASING_REF="/projects/gqilab/DAESC_GPU/data/reference/phasing/biallelic_SNV"
OUTPUT_DIR_BASE="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample5"
LOG_DIR="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample5/step3_scripts"

BARCODE_DIR="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample5/Sample5_cellbarcode"

# Create log directory if it doesn't exist
mkdir -p ${LOG_DIR}


# Loop through each sample
for SAMPLE_ID in "${SAMPLE_IDS[@]}"; do

    # Define the input BAM file for the current sample
    INPUT_BAM="${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/merged_${SAMPLE_ID}.bam"

    # Loop over each chromosome
    for CHR in "${CHROMOSOMES[@]}"; do
        echo "Processing chromosome ${CHR} for sample ${SAMPLE_ID}"

        # Define log file

        LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_chr${CHR}.log"

        # Step 1: Run GATK Genotype and log output
        #bash /projects/gqilab/DAESC_GPU/data/SALSA/step1_gatk_genotype.sh \
        #-i ${INPUT_BAM} \
        #-n ${SAMPLE_ID} \
        #-r ${REFERENCE_DIR} \
        #-g ${GATK_DIR} \
        #-d ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step1 \
        #-o com_${SAMPLE_ID}.rna.chr${CHR}.vcf.gz \
        #-l chr${CHR} \
        #-m rna \
        #-threads 6 >> ${LOG_FILE} 2>&1
        #echo "Finish step 1, chromosome ${CHR} for sample ${SAMPLE_ID}"

        # Step 2: Filter VCF file using bcftools and log output
        #bcftools view -i 'FORMAT/GQ>=30 && strlen(REF)=1 && strlen(ALT)=1' \
        #${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step1/com_${SAMPLE_ID}.rna.chr${CHR}.vcf.gz \
        #-Oz -o ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step1/filtered_com_${SAMPLE_ID}.rna.chr${CHR}.vcf.gz >> ${LOG_FILE} 2>&1

        # Index the filtered VCF file using bcftools
        #bcftools index -t ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step1/filtered_com_${SAMPLE_ID}.rna.chr${CHR}.vcf.gz >> ${LOG_FILE} 2>&1
        #echo "Finish step 2, chromosome ${CHR} for sample ${SAMPLE_ID}"

        # Step 4: Annotate VCF with GATK Funcotator
        bash /projects/gqilab/DAESC_GPU/data/SALSA/step4_gatk_anno_vcf.sh \
        --library_id pbmc_1k \
        --inputvcf ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step1/filtered_com_${SAMPLE_ID}.rna.chr${CHR}.vcf.gz \
        --outputdir ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step4 \
        --outputvcf com_${SAMPLE_ID}.pass.rna.chr${CHR}hcphase.funco.vcf.gz \
        --reference ${REFERENCE_DIR} \
        --funcotation /projects/gqilab/DAESC_GPU/data/reference/funcotator_dataSources.v1.6.20190124g \
        --output_table com_${SAMPLE_ID}.pass.rna.chr${CHR}hcphase.formatted.csv \
        --modality rna \
        --threads 6 >> ${LOG_FILE} 2>&1
        echo "Finish step 4, chromosome ${CHR} for sample ${SAMPLE_ID}"
        
        # Step 5: Filter BAM file with cell barcodes
        bash /projects/gqilab/DAESC_GPU/data/SALSA/step5_filterbam.sh \
        --library_id pbmc_1k \
        --validate \
        --inputbam ${INPUT_BAM} \
        --modality rna \
        --interval chr${CHR} \
        --barcodes ${BARCODE_DIR}/individual_${SAMPLE_ID}.csv \
        --outputdir ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step5 \
        --outputbam pbmc_1k_${SAMPLE_ID}.bcfilter.chr${CHR}.bam \
        --threads 6 >> ${LOG_FILE} 2>&1
        echo "Finish step 5, chromosome ${CHR} for sample ${SAMPLE_ID}"
        
        # Step 6: Perform variant-aware realignment with WASP
        bash /projects/gqilab/DAESC_GPU/data/SALSA/step6_wasp.sh \
        --inputvcf ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step4/com_${SAMPLE_ID}.pass.rna.chr${CHR}hcphase.funco.vcf.gz \
        --inputbam ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step5/pbmc_1k_${SAMPLE_ID}.bcfilter.chr${CHR}.bam \
        --outputdir ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step6 \
        --outputbam pbmc_1k_${SAMPLE_ID}.hcphase.chr${CHR}wasp.bam \
        --genotype rna \
        --stargenome ${REFERENCE_DIR}/star \
        --library_id pbmc_1k \
        --modality rna \
        --isphased \
        --interval chr${CHR} \
        --threads 6 >> ${LOG_FILE} 2>&1
        echo "Finish step 6, chromosome ${CHR} for sample ${SAMPLE_ID}"
        
        # Step 7: Get allele-specific read counts with SALSA
        bash /projects/gqilab/DAESC_GPU/data/SALSA/step7_gatk_alleleCount.sh \
        --inputvcf ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step4/com_${SAMPLE_ID}.pass.rna.chr${CHR}hcphase.funco.vcf.gz \
        --inputbam ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step6/pbmc_1k_${SAMPLE_ID}.hcphase.chr${CHR}wasp.bam \
        --outputdir ${OUTPUT_DIR_BASE}/com_${SAMPLE_ID}/com_${SAMPLE_ID}_step7 \
        --barcodes ${BARCODE_DIR}/individual_${SAMPLE_ID}.csv \
        --genotype rna_genotype \
        --library_id pbmc_1k \
        --modality rna \
        --reference ${REFERENCE_DIR} \
        --pseudobulk_counts \
        --single_cell_counts \
        --celltype_counts \
        --interval chr${CHR} \
        --isphased \
        --threads 6 >> ${LOG_FILE} 2>&1
        echo "Finish step 7, chromosome ${CHR} for sample ${SAMPLE_ID}"
        
    done

done
