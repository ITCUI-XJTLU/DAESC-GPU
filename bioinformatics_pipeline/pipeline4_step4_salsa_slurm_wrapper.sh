#!/bin/bash
#SBATCH --job-name=submit_pipeline4
#SBATCH --output=./logs/submit_pipeline4.log
#SBATCH --error=./logs/submit_pipeline4.err
#SBATCH --partition=gqi-32c256g,all-12c128g
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

# ====== 用户参数 ======
ROOT_PATH="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline"
REFERENCE_PATH="/projects/gqilab/DAESC_GPU/data/reference"
INPUT_LIST="./list_individuals.txt"
SIF_PATH="/home/users/tfcui23/Stat_Gene/My_SALSA/salsa_latest.sif"
SUBMIT_SCRIPT="pipeline4_step4_salsa_slurm.sh"

# ====== SLURM 资源配置 ======
PARTITION="gqi-32c256g,all-12c128g"
NTASKS=1
CPUS=5
MEMORY="120G"
TIME="7-00:00:00"

# ====== 自动计算数组范围 ======
NUM_LINES=$(wc -l < "$INPUT_LIST")
ARRAY_RANGE="0-$(($NUM_LINES - 1))"

echo "📦 样本总数: $NUM_LINES"
echo "🚀 提交数组任务: $ARRAY_RANGE"
echo "🎯 使用容器: $SIF_PATH"

# ====== 提交任务 ======
sbatch \
  --array=0-1 \
  --job-name=pipeline4_array \
  --output=./logs/pipeline4_array_%A_%a.log \
  --error=./logs/pipeline4_array_%A_%a.log \
  --partition=$PARTITION \
  --ntasks=$NTASKS \
  --cpus-per-task=$CPUS \
  --mem=$MEMORY \
  --time=$TIME \
  $SUBMIT_SCRIPT \
    root_path=$ROOT_PATH \
    reference_path=$REFERENCE_PATH \
    input_list=$INPUT_LIST \
    sif_path=$SIF_PATH

##############################################################################
# Tengfei's playground
sbatch \
  --array=0-1 \
  --job-name=pipeline4_array \
  --output=./logs/pipeline4_array_%A_%a.log \
  --error=./logs/pipeline4_array_%A_%a.log \
  --partition=gqi-32c256g,all-12c128g \
  --ntasks=1 \
  --cpus-per-task=5 \
  --mem=250G \
  --time=10:00:00 \
  pipeline4_step4_salsa_slurm.sh \
    root_path=/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline \
    reference_path=/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline/reference_test \
    input_list=/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline/list_test_individuals.txt \
    sif_path=/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/mypipeline/reference_test/salsa_latest.sif















