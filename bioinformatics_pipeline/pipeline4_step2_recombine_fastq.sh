#!/bin/bash

# ---------------- 参数解析 ----------------
root_path=""
individual_list=""
line_number=""

for ARG in "$@"; do
    case $ARG in
        root_path=*)
            root_path="${ARG#*=}"
            ;;
        individual_list=*)
            individual_list="${ARG#*=}"
            ;;
        line_number=*)
            line_number="${ARG#*=}"
            ;;
        *)
            echo "Unknown argument: $ARG"
            exit 1
            ;;
    esac
done

# 参数检查
if [[ -z "$root_path" || -z "$individual_list" || -z "$line_number" ]]; then
    echo "Usage: bash $0 root_path=... individual_list=... line_number=..."
    exit 1
fi

# 读取第 line_number 行
line=$(sed -n "${line_number}p" "$individual_list")
if [[ -z "$line" ]]; then
    echo "Line $line_number is empty or does not exist in $individual_list"
    exit 1
fi

# 提取 individual ID（忽略括号后面的东西）
individual=$(echo "$line" | awk '{print $1}')
echo "🔄 Merging FASTQ files for individual: $individual"

# 创建目标目录
indiv_dir="${root_path}/${individual}"
mkdir -p "$indiv_dir"

for read_type in R1 R2 I1; do
    output_file="${indiv_dir}/Ind${individual}_S1_${read_type}_001.fastq"
    > "$output_file"

    # 合并所有匹配的 FASTQ 文件
    find "$root_path" -type f -name "Ind${individual}_S*_${read_type}_001.fastq" | sort | while read -r file; do
        echo "  📎 Adding: $file"
        cat "$file" >> "$output_file"
    done
done

echo "✅ Done: $individual"

