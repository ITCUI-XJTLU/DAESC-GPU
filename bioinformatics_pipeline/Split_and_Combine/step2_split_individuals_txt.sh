#!/bin/bash
#SBATCH --job-name=s2_txt   # 作业名称
#SBATCH --partition=gqi-32c256g,all-12c128g
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-00:00:00           # 作业运行的最大时间 (1天)
#SBATCH --output=s2_txt.out   # 标准输出日志
#SBATCH --error=s2_txt.err    # 错误日志

# 定义原始CSV文件路径
input_file="GSM5899880_OneK1K_scRNA_Sample8_Individual_Barcodes.csv"

# 获取所有unique的individual IDs
individuals=$(awk -F, 'NR>1 {print $1}' "$input_file" | sort | uniq)

# 循环遍历每个individual
for individual in $individuals; do
  # 创建新的txt文件，文件名为 individual_{individual}.txt
  output_file="individual_${individual}.txt"
  
  # 提取当前individual的数据，并将cell barcodes格式化输出为CB:Z:{Barcode}
  awk -F, -v individual="$individual" 'NR>1 && $1 == individual {print "CB:Z:" $2}' "$input_file" > "$output_file"
  
  echo "Saved $output_file"
done
