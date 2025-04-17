#!/bin/bash

# ---------------- 参数定义 ----------------
individual_list="./list_individuals.txt"
root_path="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline"

# 自动获取行数
n_lines=$(wc -l < "$individual_list")

# ---------------- 提交 SLURM array job ----------------
sbatch --job-name=mergefq \
       --cpus-per-task=1 \
       --mem=10G \
       --time=02:00:00 \
       --partition=gqi-32c256g,all-12c128g \
       --array=1-${n_lines} \
       --output=logs/mergefq_%A_%a.log \
       --error=logs/mergefq_%A_%a.log \
       --wrap="bash merge_fastq_one_individual.sh root_path=${root_path} individual_list=${individual_list} line_number=\${SLURM_ARRAY_TASK_ID}"

