#!/bin/bash

source ~/.bashrc
conda activate daesc_qc

# 从命令行参数中获取 SRR_ID
Sample=$1
SRR_ID=$2

# 定义相对输入BAM文件路径和输出目录
input_bam="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/${Sample}/${SRR_ID}/cellranger_bam_${SRR_ID}/outs/possorted_genome_bam.bam"
output_dir="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/${Sample}/${SRR_ID}/cellranger_bam_${SRR_ID}/outs/filter_bam"
txt_dir="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/${Sample}/${Sample}_cellbarcode"

# 创建输出目录（如果不存在）
mkdir -p $output_dir
cd /projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/${Sample}/${SRR_ID}/cellranger_bam_${SRR_ID}/outs/

# 提取BAM文件的header
echo "Extracting BAM header..."
samtools view -H $input_bam > header.sam
if [ $? -ne 0 ]; then
  echo "Failed to extract BAM header. Exiting."
  exit 1
fi

# 缓存 BAM 文件到一个临时的SAM文件，避免重复读取
echo "Caching BAM file as SAM..."
samtools view $input_bam > temp.sam
if [ $? -ne 0 ]; then
  echo "Failed to cache BAM file. Exiting."
  exit 1
fi

# 循环处理txt文件目录下的所有txt文件
for txt_file in $txt_dir/individual_*.txt
do
  base_name=$(basename $txt_file .txt)
  echo "Processing $txt_file..."

  # 使用txt文件中列出的barcodes过滤缓存的SAM文件
  grep -Ff $txt_file temp.sam > ${base_name}_reads.sam
  #while read barcode; do
  #  echo "Processing barcode: $barcode"
  #  grep "$barcode" temp.sam >> ${base_name}_reads.sam
  #done < $txt_file

  if [ $? -ne 0 ]; then
    echo "Failed to filter reads for $txt_file. Skipping."
    continue
  fi

  # 合并header和过滤后的reads，并生成新的BAM文件
  cat header.sam ${base_name}_reads.sam | samtools view -bS - > $output_dir/${base_name}.bam
  if [ $? -ne 0 ]; then
    echo "Failed to create BAM file for $txt_file. Skipping."
    rm ${base_name}_reads.sam
    continue
  fi
  
  # 删除临时SAM文件
  rm ${base_name}_reads.sam
  echo "Successfully created $output_dir/${base_name}.bam"
done

# 清理临时文件
rm header.sam temp.sam
echo "Finished processing $SRR_ID"

