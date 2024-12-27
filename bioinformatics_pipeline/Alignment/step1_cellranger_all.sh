#!/bin/bash
#SBATCH --job-name=sample_process   # 作业名称
#SBATCH --partition=gqi-32c256g,all-12c128g
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-00:00:00           # 作业运行的最大时间 (1天)
#SBATCH --output=%x.out   # 标准输出日志, 以作业名命名
#SBATCH --error=%x.err    # 错误日志, 以作业名命名

# 脚本说明：这是一个针对Step1，所有Sample的一个脚本，只需要指定Sample即可
# 比如：sbatch script.sh Sample2 18028398 18028417
# 添加 Cell Ranger 到系统路径
export PATH=/home/users/tfcui23/Stat_Gene/GATK/cellranger-8.0.1:$PATH
source ~/.bashrcsource ~/.bashrc # 激活sra

# 定义 Sample 目录变量
SAMPLE_DIR="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/$1"  # 第一个参数作为样本目录
START_RUN=$2  # 第二个参数为起始 RUN_ID
END_RUN=$3    # 第三个参数为结束 RUN_ID

# 切换到指定的 Sample 目录下
cd $SAMPLE_DIR

# 循环处理指定的 RUN_ID 范围
for RUN_ID in $(seq -w $START_RUN $END_RUN); do
  RUN="SRR${RUN_ID}"
  echo "Processing $RUN"

  # 下载 SRA 数据
  prefetch -O $SAMPLE_DIR $RUN

  # 切换到下载好的样本文件夹
  cd $SAMPLE_DIR/$RUN

  # 将 .sra 文件转换为 FASTQ 文件
  fastq-dump --split-files $SAMPLE_DIR/$RUN/$RUN.sra

  # 重命名 .fastq 文件
  mv ${RUN}_1.fastq ${RUN}_S1_L001_I1_001.fastq
  mv ${RUN}_2.fastq ${RUN}_S1_L001_R1_001.fastq
  mv ${RUN}_3.fastq ${RUN}_S1_L001_R2_001.fastq

  # 运行 Cell Ranger count 命令
  cellranger count --id=cellranger_bam_${RUN} \
                   --transcriptome=/projects/gqilab/DAESC_GPU/data/reference/refdata-gex-GRCh38-2020-A \
                   --fastqs=$SAMPLE_DIR/$RUN \
                   --sample=$RUN \
                   --create-bam=true \
                   --localcores=3 \
                   --nosecondary

  # 返回 Sample 目录以便处理下一个样本
  cd $SAMPLE_DIR

done

