#!/bin/bash

# -----------------------------
# 参数解析
# -----------------------------
for arg in "$@"; do
    case $arg in
        list_file=*)
            LIST_FILE="${arg#*=}"
            shift
            ;;
        root_path=*)
            ROOT_PATH="${arg#*=}"
            shift
            ;;
        *)
            echo "Unknown argument: $arg"
            echo "Usage: ./split_by_csvlist.sh list_file=your_list.txt root_path=/your/output/dir"
            exit 1
            ;;
    esac
done

# -----------------------------
# 参数检查
# -----------------------------
if [ -z "$LIST_FILE" ] || [ -z "$ROOT_PATH" ]; then
    echo "Error: list_file and root_path must be specified."
    exit 1
fi

if [ ! -f "$LIST_FILE" ]; then
    echo "Error: list file '$LIST_FILE' does not exist."
    exit 1
fi

# 创建输出目录
OUTPUT_DIR="${ROOT_PATH}/cellbarcodes"
mkdir -p "$OUTPUT_DIR"

# -----------------------------
# 获取唯一 CSV 路径
# -----------------------------
csv_paths=$(awk '{print $3}' "$LIST_FILE" | sort | uniq)

# -----------------------------
# 生成新的 summary TXT 文件
# -----------------------------
SUMMARY_OUTPUT="${ROOT_PATH}/list_test_srr.txt"
> "$SUMMARY_OUTPUT"  # 清空或创建文件

while read -r line; do
    sample=$(echo "$line" | awk '{print $1}')
    srr=$(echo "$line" | awk '{print $2}')
    csv_path=$(echo "$line" | awk '{print $3}')

    if [ ! -f "$csv_path" ]; then
        echo "Warning: File '$csv_path' not found, skipping." >&2
        continue
    fi

    # 提取第一列（individual ID），去重并清理非字母数字字符
    cleaned_ids=$(awk -F, 'NR>1 {gsub(/[^a-zA-Z0-9]/, "", $1); print $1}' "$csv_path" | sort | uniq | tr '\n' ' ')

    # 去掉尾部多余的空格
    cleaned_ids=$(echo "$cleaned_ids" | sed 's/ *$//')

    # 拼接成目标格式
    echo "$sample $srr (${cleaned_ids})" >> "$SUMMARY_OUTPUT"

done < "$LIST_FILE"

echo "Saved summary to $SUMMARY_OUTPUT"


SUMMARY_INPUT="${ROOT_PATH}/list_test_srr.txt"
NEW_OUTPUT="${ROOT_PATH}/list_test_individuals.txt"

# 构造目标括号内容
CHR_LIST=$(seq 1 22 | tr '\n' ' ')
CHR_LIST=$(echo "$CHR_LIST" | sed 's/ *$//')

# 提取括号内所有 ID 并去重
all_ids=$(grep -oP '\([^)]+\)' "$SUMMARY_INPUT" | tr -d '()' | tr ' ' '\n' | sort | uniq)

# 输出结果
> "$NEW_OUTPUT"
for id in $all_ids; do
    echo "$id ($CHR_LIST)" >> "$NEW_OUTPUT"
done

echo "Saved mapping to $NEW_OUTPUT"


# -----------------------------
# 对每个 CSV 文件进行处理
# -----------------------------
for INPUT_FILE in $csv_paths; do
    if [ ! -f "$INPUT_FILE" ]; then
        echo "Warning: CSV file '$INPUT_FILE' does not exist, skipping."
        continue
    fi

    echo "Processing $INPUT_FILE..."

    # 获取 unique individual IDs
    individuals=$(awk -F, 'NR>1 {print $1}' "$INPUT_FILE" | sort | uniq)

    for individual in $individuals; do
        cleaned_id=$(echo "$individual" | tr -d '_')

        csv_output="${OUTPUT_DIR}/${cleaned_id}.csv"
        txt_output="${OUTPUT_DIR}/${cleaned_id}.txt"

        # 生成 CSV 文件
        awk -F, -v individual="$individual" '
            BEGIN {
                OFS=",";
                print "barcode","orig.ident","celltype"
            }
            NR>1 && $1 == individual {
                sub(/-1$/, "", $2);
                print $2,"pbmc_1k",1
            }' "$INPUT_FILE" > "$csv_output"

        # 生成 TXT 文件
        awk -F, -v individual="$individual" '
            NR>1 && $1 == individual {
                print "CB:Z:" $2
            }' "$INPUT_FILE" > "$txt_output"

        echo "Saved ${csv_output} and ${txt_output}"
    done
done

