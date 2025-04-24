#!/bin/bash

# ------------------ 解析命名参数 -------------------
root_path=""
txt_1=""
barcodes_dir=""
line_number=""

for ARG in "$@"; do
    case $ARG in
        root_path=*)
            root_path="${ARG#*=}"
            ;;
        txt_1=*)
            txt_1="${ARG#*=}"
            ;;
        barcodes_dir=*)
            barcodes_dir="${ARG#*=}"
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

# ------------------ 参数检查 -------------------
if [[ -z "$root_path" || -z "$txt_1" || -z "$line_number" ]]; then
    echo "Usage: bash $0 root_path=... txt_1=... line_number=... [barcodes_dir=...]"
    exit 1
fi

# ------------------ 设置 Conda 环境 PATH -------------------
export PATH="${root_path}/reference_test/miniconda3/envs/scASE_conda/bin:$PATH"
echo "🔧 PATH set to use Conda environment in ${root_path}/reference_test/miniconda3/envs/scASE_conda"

# 读取指定行
line=$(sed -n "${line_number}p" "$txt_1")
if [[ -z "$line" ]]; then
    echo "❌ Line $line_number is empty or does not exist in $txt_1"
    exit 1
fi

read -r sample srr_id individuals_raw <<< "$line"
individuals=$(echo "$individuals_raw" | sed -E 's/[\(\)]//g')

# 路径设置
data_dir="${root_path}/${sample}"
fastq_dir="${data_dir}/${srr_id}"
split_fastq_dir="${fastq_dir}/split_fastq"
result_dir="${split_fastq_dir}/splited_result"

mkdir -p "$split_fastq_dir" "$result_dir"
rm -f "$split_fastq_dir"/* "$result_dir"/*

# FASTQ 前缀
sample_id=$(echo "$sample" | sed 's/Sample/S/')
fastq_prefix="${sample}_${sample_id}"

echo "[$srr_id] 生成 reads_information.txt"
awk 'NR%4==1 {srr=$1} NR%4==2 {print srr"\t"$1}' \
    "${fastq_dir}/${fastq_prefix}_R1_001.fastq" > "$split_fastq_dir/reads_information.txt"

for ind in $individuals; do
    bc_file="${barcodes_dir}/${ind}.csv"
    if [[ ! -f "$bc_file" ]]; then
        echo "❗ Barcode file not found for individual $ind, skipping."
        continue
    fi

    echo "[$srr_id][$ind] 生成 cell_barcodes_only_${ind}.txt"
    tail -n +2 "$bc_file" | cut -d',' -f1 > "$split_fastq_dir/cell_barcodes_only_${ind}.txt"

    echo "[$srr_id][$ind] 提取 reads_information_${ind}.txt"
    awk 'NR==FNR {bc[$1]; next} {cell_bc=substr($2, 1, 16); if (cell_bc in bc) print $0}' \
        "$split_fastq_dir/cell_barcodes_only_${ind}.txt" "$split_fastq_dir/reads_information.txt" > \
        "$split_fastq_dir/reads_information_${ind}.txt"

    cut -f1 "$split_fastq_dir/reads_information_${ind}.txt" | sed 's/^@//' > "$split_fastq_dir/read_ids_${ind}.txt"

    echo "[$srr_id][$ind] split FASTQ：R1, R2, I1"
    for read_type in R1 I1 R2; do
        filterbyname.sh \
            in="${fastq_dir}/${fastq_prefix}_${read_type}_001.fastq" \
            out="${result_dir}/Ind${ind}_${sample_id}_${read_type}_001.fastq" \
            names="$split_fastq_dir/read_ids_${ind}.txt" \
            include=t
    done
done

echo "✅ Finished SRR $srr_id"

