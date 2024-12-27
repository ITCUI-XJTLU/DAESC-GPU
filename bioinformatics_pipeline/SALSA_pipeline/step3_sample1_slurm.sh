#!/bin/bash
#SBATCH --job-name=st3_s1
#SBATCH --output=sample1_SALSA.log
#SBATCH --error=sample1_SALSA.log
#SBATCH --partition=all-12c128g
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00 
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=16G


# 运行 Singularity 容器并执行脚本, 原本的servers：gqi-32c256g,all-12c128g
singularity exec \
--bind $HOME:/projects/gqilab/DAESC_GPU \
--bind $project:/projects/gqilab/DAESC_GPU/data/Test_run_pool1 \
--bind $reference:/projects/gqilab/DAESC_GPU/data/reference \
--bind $SALSA:/projects/gqilab/DAESC_GPU/data/SALSA \
--bind $Barcodes:/projects/gqilab/DAESC_GPU/data/cellbarcodes \
--bind $SCRATCH1:/projects/gqilab/DAESC_GPU/data/Test_run_pool1/scratch \
--env SCRATCH1="/projects/gqilab/DAESC_GPU/data/Test_run_pool1/scratch" \
/home/users/tfcui23/Stat_Gene/My_SALSA/salsa_latest.sif \
bash /projects/gqilab/DAESC_GPU/data/Test_run_pool1/Scripts/step3_sample1.sh >> sample1_SALSA.log 2>&1








