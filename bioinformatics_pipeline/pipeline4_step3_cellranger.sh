#!/bin/bash

# ---------------- 参数解析 ----------------
root_path=""
individual_list=""
reference_path=""
line_number=""
cellranger_path=""

for ARG in "$@"; do
    case $ARG in
        root_path=*)
            root_path="${ARG#*=}"
            ;;
        individual_list=*)
            individual_list="${ARG#*=}"
            ;;
        reference_path=*)
            reference_path="${ARG#*=}"
            ;;
        line_number=*)
            line_number="${ARG#*=}"
            ;;
        cellranger_path=*)
            cellranger_path="${ARG#*=}"
            ;;
        *)
            echo "❌ Unknown argument: $ARG"
            echo "Usage: bash $0 root_path=... individual_list=... reference_path=... line_number=... cellranger_path=..."
            exit 1
            ;;
    esac
done

# ---------------- 参数检查 ----------------
if [[ -z "$root_path" || -z "$individual_list" || -z "$reference_path" || -z "$line_number" || -z "$cellranger_path" ]]; then
    echo "❌ Missing one or more required arguments."
    echo "Usage: bash $0 root_path=... individual_list=... reference_path=... line_number=... cellranger_path=..."
    exit 1
fi

# ---------------- 获取当前行的 individual ID ----------------
line=$(sed -n "${line_number}p" "$individual_list")
if [[ -z "$line" ]]; then
    echo "❌ Line $line_number is empty or does not exist in $individual_list"
    exit 1
fi
individual=$(echo "$line" | awk '{print $1}' | sed 's/[:,]//g')  # 清理非法字符

# ---------------- 设置路径 ----------------
indiv_dir="${root_path}/${individual}"
fastq_dir="$indiv_dir"
sample_name="Ind${individual}"
output_subdir="${fastq_dir}/${sample_name}"
log_dir="${root_path}/logs"
log_file="${log_dir}/cellranger_${sample_name}.log"

mkdir -p "$log_dir"
cd ${root_path}/${individual}

# ---------------- 执行 CellRanger ----------------
echo "🧬 Running CellRanger for: $individual"
echo "📁 FASTQ directory: $fastq_dir"
echo "📦 Sample name: $sample_name"
echo "📚 Reference: ${reference_path}/refdata-gex-GRCh38-2020-A"
echo "🚀 CellRanger path: $cellranger_path"
echo "📝 Log file: $log_file"

"${cellranger_path}/cellranger" count \
    --id="$sample_name" \
    --transcriptome="${reference_path}/refdata-gex-GRCh38-2020-A" \
    --fastqs="$fastq_dir" \
    --sample="$sample_name" \
    --create-bam=true \
    --localcores=3 \
    --nosecondary \
    &> "$log_file"

echo "✅ Finished: $individual"

