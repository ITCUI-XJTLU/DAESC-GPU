#!/bin/bash
#SBATCH --job-name=st3_s3
#SBATCH --output=sample3_SALSA.log
#SBATCH --error=sample3_SALSA.log
#SBATCH --partition=gqi-32c256g
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00           # 作业运行的最大时间
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=16G

# 运行 Singularity 容器并执行脚本
singularity exec \
--bind $HOME:/projects/gqilab/DAESC_GPU \
--bind $project:/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample3 \
--bind $reference:/projects/gqilab/DAESC_GPU/data/reference \
--bind $SALSA:/projects/gqilab/DAESC_GPU/data/SALSA \
--bind $Barcodes:/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample3/Sample3_cellbarcode \
--bind $SCRATCH1:/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample3/step3_scratch \
--env SCRATCH1="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample3/step3_scratch" \
/home/users/tfcui23/Stat_Gene/My_SALSA/salsa_latest.sif \
bash /projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample3/step3_scripts/step3_sample3.sh >> sample3_SALSA.log 2>&1








