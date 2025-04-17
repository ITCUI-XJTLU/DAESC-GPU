#!/bin/bash
set -euo pipefail

Individual_ID=$1
REFERENCE_DIR=$2
OUTPUT_DIR_BASE=$3
LOG_DIR=$4
BARCODE_DIR=$5
SCRATCH1=$6
CHR_STRING="$7"

echo "🚀 Start pipeline for $1"

IFS=' ' read -r -a CHR_LIST <<< "$CHR_STRING"

INPUT_BAM="${OUTPUT_DIR_BASE}/Ind${Individual_ID}/outs/possorted_genome_bam.bam"
BARCODE_FILE="${BARCODE_DIR}/${Individual_ID}.csv"
mkdir -p "$LOG_DIR"

for CHR in "${CHR_LIST[@]}"; do
    echo "Processing chr${CHR} for individual ${Individual_ID}"
    LOG_FILE="${LOG_DIR}/${Individual_ID}_chr${CHR}.log"

    # Step 1: Genotype
    bash ./salsa1_gatk_genotype.sh \
        -i ${INPUT_BAM} \
        -n ${Individual_ID} \
        -r ${REFERENCE_DIR}/refdata-gex-GRCh38-2020-A \
        -g ${REFERENCE_DIR}/gatk \
        -d ${OUTPUT_DIR_BASE}/${Individual_ID}_step1 \
        -o ${Individual_ID}.rna.chr${CHR}.vcf.gz \
        -l chr${CHR} \
       -m rna \
        -threads 5 >> ${LOG_FILE} 2>&1

    # Step 1: Filter VCF
    bcftools view -i 'FORMAT/GQ>=30 && strlen(REF)=1 && strlen(ALT)=1' \
        ${OUTPUT_DIR_BASE}/${Individual_ID}_step1/${Individual_ID}.rna.chr${CHR}.vcf.gz \
        -Oz -o ${OUTPUT_DIR_BASE}/${Individual_ID}_step1/filtered_${Individual_ID}.rna.chr${CHR}.vcf.gz >> ${LOG_FILE} 2>&1

    bcftools index -t ${OUTPUT_DIR_BASE}/${Individual_ID}_step1/filtered_${Individual_ID}.rna.chr${CHR}.vcf.gz >> ${LOG_FILE} 2>&1

    # Step 2: Phasing
    bash ./salsa2_phase_vcf_shapeit5.sh \
        --library_id pbmc_1k \
        --inputvcf ${OUTPUT_DIR_BASE}/${Individual_ID}_step1/filtered_${Individual_ID}.rna.chr${CHR}.vcf.gz \
        --outputdir ${OUTPUT_DIR_BASE}/${Individual_ID}_step3 \
        --outputvcf pipe3_${Individual_ID}.pass.rna.chr${CHR}hcphase.vcf.gz \
        --phasingref ${REFERENCE_DIR} \
        --interval chr${CHR} \
        --hcphase \
        --snvonly \
        --verbose \
        --threads 5 >> ${LOG_FILE} 2>&1

    # Step 3: Annotation
    bash ./salsa3_gatk_anno_vcf.sh \
        --library_id pbmc_1k \
        --inputvcf ${OUTPUT_DIR_BASE}/${Individual_ID}_step3/pipe3_${Individual_ID}.pass.rna.chr${CHR}hcphase.vcf.gz \
        --outputdir ${OUTPUT_DIR_BASE}/${Individual_ID}_step4 \
        --outputvcf pipe3_${Individual_ID}.pass.rna.chr${CHR}hcphase.funco.vcf.gz \
        --reference ${REFERENCE_DIR}/refdata-gex-GRCh38-2020-A \
        --funcotation ${REFERENCE_DIR}/funcotator_dataSources.v1.6.20190124g \
        --output_table pipe3_${Individual_ID}.pass.rna.chr${CHR}hcphase.formatted.csv \
        --modality rna \
        --threads 5 >> ${LOG_FILE} 2>&1

    # Step 4: Barcode Filter
    bash ./salsa4_filterbam.sh \
        --library_id pbmc_1k \
        --validate \
        --inputbam ${INPUT_BAM} \
        --modality rna \
        --interval chr${CHR} \
        --barcodes ${BARCODE_FILE} \
        --outputdir ${OUTPUT_DIR_BASE}/${Individual_ID}_step5 \
        --outputbam pbmc_1k_${Individual_ID}.bcfilter.chr${CHR}.bam \
        --threads 5 >> ${LOG_FILE} 2>&1

    # Step 5: WASP
    bash ./salsa5_wasp.sh \
        --inputvcf ${OUTPUT_DIR_BASE}/${Individual_ID}_step4/pipe3_${Individual_ID}.pass.rna.chr${CHR}hcphase.funco.vcf.gz \
        --inputbam ${OUTPUT_DIR_BASE}/${Individual_ID}_step5/pbmc_1k_${Individual_ID}.bcfilter.chr${CHR}.bam \
        --outputdir ${OUTPUT_DIR_BASE}/${Individual_ID}_step6 \
        --outputbam pbmc_1k_${Individual_ID}.hcphase.chr${CHR}wasp.bam \
        --genotype rna \
        --stargenome ${REFERENCE_DIR}/refdata-gex-GRCh38-2020-A/star \
        --chromInfo ${REFERENCE_DIR}/hg38_chromInfo.txt.gz \
        --library_id pbmc_1k \
        --modality rna \
        --isphased \
        --interval chr${CHR} \
        --threads 5 >> ${LOG_FILE} 2>&1

    # Step 6: Allele-specific count
    bash ./salsa6_gatk_alleleCount.sh \
        --inputvcf ${OUTPUT_DIR_BASE}/${Individual_ID}_step4/pipe3_${Individual_ID}.pass.rna.chr${CHR}hcphase.funco.vcf.gz \
        --inputbam ${OUTPUT_DIR_BASE}/${Individual_ID}_step6/pbmc_1k_${Individual_ID}.hcphase.chr${CHR}wasp.bam \
        --outputdir ${OUTPUT_DIR_BASE}/${Individual_ID}_step7 \
        --barcodes ${BARCODE_FILE} \
        --genotype rna_genotype \
        --library_id pbmc_1k \
        --modality rna \
        --reference ${REFERENCE_DIR}/refdata-gex-GRCh38-2020-A \
        --pseudobulk_counts \
        --single_cell_counts \
        --celltype_counts \
        --interval chr${CHR} \
        --isphased \
        --threads 5 >> ${LOG_FILE} 2>&1

done

