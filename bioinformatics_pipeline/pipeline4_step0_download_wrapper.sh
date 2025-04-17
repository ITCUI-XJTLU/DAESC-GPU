#!/bin/bash
#SBATCH --job-name=dl_sra_array
#SBATCH --output=logs/dl_sra_%A_%a.out
#SBATCH --error=logs/dl_sra_%A_%a.err
#SBATCH --array=1-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=gqi-32c256g

# 参数设置
root_path="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline"
txt_1="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline/list_test_srr.txt"
conda_activate_cmd="source /home/users/tfcui23/Stat_Gene/GATK/miniconda/etc/profile.d/conda.sh"
line_number=${SLURM_ARRAY_TASK_ID}

# 运行主脚本
bash pipeline4_step0_download.sh root_path="$root_path" txt_1="$txt_1" line_number="$line_number" conda_activate_cmd="$conda_activate_cmd"

