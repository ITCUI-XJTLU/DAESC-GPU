#!/bin/bash
#SBATCH --job-name=s2_csv   # 作业名称
#SBATCH --partition=gqi-32c256g,all-12c128g
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-00:00:00           # 作业运行的最大时间 (1天)
#SBATCH --output=s2_csv.out   # 标准输出日志
#SBATCH --error=s2_csv.err    # 错误日志

# 定义原始CSV文件路径
input_file="GSM5899882_OneK1K_scRNA_Sample10_Individual_Barcodes.csv"

# 获取所有unique的individual IDs
individuals=$(awk -F, 'NR>1 {print $1}' "$input_file" | sort | uniq)

# 循环遍历每个individual
for individual in $individuals; do
  # 创建新文件，文件名为 individual_{individual}.csv
  output_file="individual_${individual}.csv"
  
  # 提取当前individual的数据并去掉barcode末尾的"-1"
  awk -F, -v individual="$individual" 'BEGIN {OFS=","; print "barcode","orig.ident","celltype"}
      NR>1 && $1 == individual {gsub(/-1$/, "", $2); print $2,"pbmc_1k",1}' "$input_file" > "$output_file"
  
  
  echo "Saved $output_file"
done
