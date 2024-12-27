#!/bin/bash
#SBATCH --job-name=sample2   # 作业名称
#SBATCH --partition=gqi-32c256g,all-12c128g
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-00:00:00           # 作业运行的最大时间 (1天)
#SBATCH --output=sample2.out   # 标准输出日志
#SBATCH --error=sample2.err    # 错误日志

# 添加 Cell Ranger 到系统路径
export PATH=/home/users/tfcui23/Stat_Gene/GATK/cellranger-8.0.1:$PATH
source ~/.bashrc # 激活sra

# 切换到 Sample2 路径下
cd /projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample2

# 循环处理 SRR18028398 到 SRR18028417
for RUN_ID in $(seq -w 18028398 18028417); do
  RUN="SRR${RUN_ID}"
  echo "Processing $RUN"

  # 下载 SRA 数据
  prefetch -O /projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample2 $RUN

  # 切换到下载好的样本文件夹
  cd /projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample2/$RUN

  # 将 .sra 文件转换为 FASTQ 文件
  fastq-dump --split-files /projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample2/$RUN/$RUN.sra

  # 重命名 .fastq 文件
  mv ${RUN}_1.fastq ${RUN}_S1_L001_I1_001.fastq
  mv ${RUN}_2.fastq ${RUN}_S1_L001_R1_001.fastq
  mv ${RUN}_3.fastq ${RUN}_S1_L001_R2_001.fastq

  # 运行 Cell Ranger count 命令
  cellranger count --id=cellranger_bam_${RUN} \
                   --transcriptome=/projects/gqilab/DAESC_GPU/data/reference/refdata-gex-GRCh38-2020-A \
                   --fastqs=/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample2/$RUN \
                   --sample=$RUN \
                   --create-bam=true \
                   --localcores=3 \
                   --nosecondary

  # 返回 Sample2 目录以便处理下一个样本
  cd /projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample2

done

