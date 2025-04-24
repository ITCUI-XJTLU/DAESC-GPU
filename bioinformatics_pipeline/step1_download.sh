#!/bin/bash

# =============================
# 👇 参数解析
# =============================
root_path=""
txt_1=""
line_number=""

for ARG in "$@"; do
    case $ARG in
        root_path=*)
            root_path="${ARG#*=}"
            ;;
        txt_1=*)
            txt_1="${ARG#*=}"
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

# =============================
# 参数检查
# =============================
if [[ -z "$root_path" || -z "$txt_1" || -z "$line_number" ]]; then
    echo "Usage: bash $0 root_path=... txt_1=... line_number=..."
    exit 1
fi

# =============================
# 👇 设置 Conda 环境路径
# =============================
export PATH="${root_path}/reference_test/miniconda3/envs/scASE_conda/bin:$PATH"
echo "🔧 PATH set to use Conda environment in ${root_path}/reference_test/miniconda3/envs/scASE_conda"

# =============================
# 下载逻辑
# =============================
echo "📥 Starting download..."
line=$(sed -n "${line_number}p" "$txt_1")
if [[ -z "$line" ]]; then
    echo "❌ Line $line_number is empty or does not exist in $txt_1"
    exit 1
fi

read -r SAMPLE SRR_ID _ <<< "$line"

SAMPLE_DIR="${root_path}/${SAMPLE}"
mkdir -p "$SAMPLE_DIR"

echo "🚀 Downloading $SRR_ID for $SAMPLE"
prefetch "$SRR_ID" -O "$SAMPLE_DIR"

SRR_DIR="${SAMPLE_DIR}/${SRR_ID}"
cd "$SRR_DIR" || { echo "❌ Failed to cd into $SRR_DIR"; exit 1; }

fastq-dump --split-files "$SRR_ID.sra"

SAMPLE_ID=$(echo "$SAMPLE" | sed 's/Sample/S/')
[[ -f "${SRR_ID}_1.fastq" ]] && mv "${SRR_ID}_1.fastq" "${SAMPLE}_${SAMPLE_ID}_I1_001.fastq"
[[ -f "${SRR_ID}_2.fastq" ]] && mv "${SRR_ID}_2.fastq" "${SAMPLE}_${SAMPLE_ID}_R1_001.fastq"
[[ -f "${SRR_ID}_3.fastq" ]] && mv "${SRR_ID}_3.fastq" "${SAMPLE}_${SAMPLE_ID}_R2_001.fastq"

rm -f "$SRR_ID.sra"
echo "✅ Done: $SRR_ID"

