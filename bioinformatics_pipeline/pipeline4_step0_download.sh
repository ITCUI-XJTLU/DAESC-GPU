#!/bin/bash
# =============================
# 👇 参数解析 + 支持 * 表示空格
# =============================
conda_activate_cmd=""
root_path=""
txt_1=""
line_number=""

for ARG in "$@"; do
    case $ARG in
        conda_activate_cmd=*)
            conda_activate_cmd="${ARG#*=}"
            ;;
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
if [[ -z "$conda_activate_cmd" || -z "$root_path" || -z "$txt_1" || -z "$line_number" ]]; then
    echo "Usage: bash $0 conda_activate_cmd=source*/path/to/conda.sh root_path=... txt_1=... line_number=..."
    exit 1
fi

# =============================
# 👇 将 * 替换为空格后执行激活命令
# =============================
conda_activate_cmd="${conda_activate_cmd//\*/ }"
echo "activate conda by $conda_activate_cmd"
eval "$conda_activate_cmd"
conda activate scASE

# =============================
# 下载逻辑保持不变
# =============================
echo "start to download"
line=$(sed -n "${line_number}p" "$txt_1")
if [[ -z "$line" ]]; then
    echo "Line $line_number is empty or does not exist in $txt_1"
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

