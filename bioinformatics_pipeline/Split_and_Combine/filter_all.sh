#!/bin/bash

# 定义输入BAM文件路径和输出目录
input_bam="/projects/gqilab/DAESC_GPU/data/Test_run1/SRR18028384/cellranger_bam/outs/possorted_genome_bam.bam"
output_dir="/projects/gqilab/DAESC_GPU/data/Test_run1/SRR18028384/cellranger_bam/outs/filter_bam"

# 创建输出目录（如果不存在）
# mkdir -p $output_dir

# 提取BAM文件的header
samtools view -H $input_bam > header.sam

# 循环处理所有txt文件
for txt_file in ind_*.txt
do
  base_name=$(basename $txt_file .txt)
  echo "Processing $txt_file"
  samtools view $input_bam | grep -Ff $txt_file > ${base_name}_reads.sam
  cat header.sam ${base_name}_reads.sam | samtools view -bS - > $output_dir/${base_name}.bam
  rm ${base_name}_reads.sam
done

# 清理临时文件
rm header.sam
