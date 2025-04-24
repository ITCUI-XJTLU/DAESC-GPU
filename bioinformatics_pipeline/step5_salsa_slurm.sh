#!/bin/bash

# ====== SLURM array job header (日志输出) ======
#SBATCH --job-name=sl
#SBATCH --output=./logs/pipeline4_%A_%a.log
#SBATCH --error=./logs/pipeline4_%A_%a.log

# ====== 参数解析 ======
for arg in "$@"; do
    case $arg in
        root_path=*)
            ROOT_PATH="${arg#*=}" ;;
        input_list=*)
            INPUT_LIST="${arg#*=}" ;;
        *)
            echo "❌ 未知参数：$arg"
            exit 1 ;;
    esac
done

# ====== 参数检查 ======
if [[ -z "${ROOT_PATH:-}" || -z "${INPUT_LIST:-}" ]]; then
    echo "❌ 错误：必须传入参数 root_path=... input_list=..."
    exit 1
fi

# ====== 推导默认路径 ======
REFERENCE_PATH="${ROOT_PATH}/reference"
SIF_PATH="${ROOT_PATH}/reference_test/salsa_latest.sif"

# ====== 检查路径是否存在 ======
if [[ ! -f "$SIF_PATH" ]]; then
    echo "❌ 错误：未找到 SIF 文件：$SIF_PATH"
    exit 1
fi

if [[ ! -d "$REFERENCE_PATH" ]]; then
    echo "❌ 错误：未找到参考基因组目录：$REFERENCE_PATH"
    exit 1
fi

# ====== 读取当前任务行 ======
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$INPUT_LIST")
Individual_ID=$(echo "$LINE" | awk '{print $1}')
CHR_LIST=$(echo "$LINE" | sed -e 's/^[^ ]* //' -e 's/[()]//g')

# ====== 构建路径 ======
OUTPUT_DIR_BASE="${ROOT_PATH}/${Individual_ID}"
BARCODE_DIR="${ROOT_PATH}/barcodes"
LOG_DIR="${OUTPUT_DIR_BASE}/logs_${Individual_ID}"
SCRATCH1="${OUTPUT_DIR_BASE}/scratch_${Individual_ID}"

mkdir -p "$LOG_DIR" "$SCRATCH1"

# ====== 执行主分析流程 ======
singularity exec \
--bind $ROOT_PATH:$ROOT_PATH \
--env SCRATCH1="$SCRATCH1" \
"$SIF_PATH" \
bash "$ROOT_PATH/pipeline4_step4_salsa.sh" \
"$Individual_ID" "$REFERENCE_PATH" "$OUTPUT_DIR_BASE" "$LOG_DIR" "$BARCODE_DIR" "$SCRATCH1" "$CHR_LIST"

#singularity exec \
#--bind $ROOT_PATH:$ROOT_PATH \
#--bind $LOG_DIR:$LOG_DIR \
#--bind $SCRATCH1:$SCRATCH1 \
#--bind $BARCODE_DIR:$BARCODE_DIR \
#--env SCRATCH1="$SCRATCH1" \
#"$SIF_PATH" \
#bash "$ROOT_PATH/pipeline4_step4_salsa.sh" \
#"$Individual_ID" "$REFERENCE_PATH" "$OUTPUT_DIR_BASE" "$LOG_DIR" "$BARCODE_DIR" "$SCRATCH1" "$CHR_LIST"

