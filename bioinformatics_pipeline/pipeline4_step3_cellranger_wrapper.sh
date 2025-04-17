#!/bin/bash

# ---------------- 参数定义 ----------------
individual_list="./list_individuals.txt"
root_path="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline"
reference_path="/projects/gqilab/references"  # 替换成你的参考路径

# 自动获取行数
n_lines=$(wc -l < "$individual_list")

# ---------------- 提交 SLURM array job ----------------
sbatch --job-name=cellranger \
       --cpus-per-task=3 \
       --mem=60G \
       --time=2-00:00:00 \
       --partition=gqi-32c256g,all-12c128g \
       --array=1-${n_lines} \
       --output=logs/cellranger_%A_%a.log \
       --error=logs/cellranger_%A_%a.log \
       --wrap="bash run_cellranger_one_individual.sh root_path=${root_path} individual_list=${individual_list} reference_path=${reference_path} line_number=\${SLURM_ARRAY_TASK_ID}"

