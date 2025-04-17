#!/bin/bash

txt_1="./srr_list.txt"
n_lines=$(wc -l < "$txt_1")

root_path="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline"
barcodes_dir="${root_path}/barcodes"

# 注意这里的单引号和双引号配合使用
sbatch --job-name=splitfq \
       --cpus-per-task=3 \
       --mem=50G \
       --time=05:00:00 \
       --partition=gqi-32c256g,all-12c128g \
       --array=1-${n_lines} \
       --output=logs/splitfq_%A_%a.log \
       --error=logs/splitfq_%A_%a.log \
       --wrap="bash pipeline4_step1_split_fastq.sh root_path=${root_path} txt_1=${txt_1} barcodes_dir=${barcodes_dir} line_number=\${SLURM_ARRAY_TASK_ID}"

