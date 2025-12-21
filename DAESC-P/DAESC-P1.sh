#!/bin/bash
# =============================================================================
# DAESC-P1.sh - SRA to FASTQ Conversion and Per-Individual Splitting
# =============================================================================
# Purpose:
#   1. (Optional) Download SRA file using prefetch
#   2. Convert SRA files to compressed FASTQ.gz format using fastq-dump
#   3. Intelligently rename FASTQ files based on 10X Genomics version detection
#   4. Split FASTQ files by individual based on cell barcodes
#
# Usage: 
#   bash DAESC-P1.sh root_path=... txt_input=... line_number=... cellbarcodes_file=... [download_srr=True/False]
#
# Input:
#   - SRA file: Data_SRR/{SRR_ID}/{SRR_ID}.sra (or will be downloaded)
#   - cellbarcodes_file: CSV with (individual_id, barcode, ...) columns
#
# Output:
#   - FASTQ files: {SRR_ID}_S1_{I1,R1,R2}_001.fastq.gz
#   - Split FASTQ files: split_fastq/{SRR_ID}_Ind{individual}_{I1,R1,R2}_001.fastq.gz
#
# Dependencies:
#   - SRA-Toolkit (prefetch, fastq-dump)
#   - BBMap (filterbyname.sh)
# =============================================================================

# =============================
# Argument Parsing
# =============================
root_path=""
txt_input=""
line_number=""
cellbarcodes_file=""
download_srr="False"

for ARG in "$@"; do
    case $ARG in
        root_path=*) root_path="${ARG#*=}" ;;
        txt_input=*) txt_input="${ARG#*=}" ;;
        line_number=*) line_number="${ARG#*=}" ;;
        cellbarcodes_file=*) cellbarcodes_file="${ARG#*=}" ;;
        download_srr=*) download_srr="${ARG#*=}" ;;
        *) echo "❌ Unknown argument: $ARG"; exit 1 ;;
    esac
done

# =============================
# Argument Validation
# =============================
if [[ -z "$root_path" || -z "$txt_input" || -z "$line_number" || -z "$cellbarcodes_file" ]]; then
    echo "Usage: bash $0 root_path=... txt_input=... line_number=... cellbarcodes_file=... [download_srr=True/False]"
    exit 1
fi

# =============================
# Read SRR ID from Specified Line
# =============================
SRR_ID=$(sed -n "${line_number}p" "$txt_input")

if [[ -z "$SRR_ID" ]]; then
    echo "❌ Line $line_number is empty or does not exist"
    exit 1
fi

# Remove any whitespace characters from SRR ID
SRR_ID=$(echo "$SRR_ID" | tr -d '[:space:]')

echo "=========================================="
echo "📋 Task $line_number: Processing $SRR_ID"
echo "=========================================="
echo "📥 Download SRR: $download_srr"
echo ""

# =============================
# Path Configuration
# =============================
DATA_DIR="${root_path}/Data_SRR"
SRR_DIR="${DATA_DIR}/${SRR_ID}"
SRA_FILE="${SRR_DIR}/${SRR_ID}.sra"

# Create directory if it doesn't exist
mkdir -p "$SRR_DIR"

# =============================
# Step 0: Download SRA File (Optional)
# =============================
if [[ "$download_srr" == "True" || "$download_srr" == "true" || "$download_srr" == "TRUE" ]]; then
    echo "=========================================="
    echo "📌 Step 0: Download SRA file"
    echo "=========================================="
    
    # Check if already downloaded (skip to save time)
    if [[ -f "$SRA_FILE" ]]; then
        echo "⏭️  Skipping $SRR_ID (already exists)"
        echo "📁 SRA file: $SRA_FILE"
    else
        # Use prefetch to download SRA file with progress display
        echo "📥 Downloading $SRR_ID..."
        prefetch "$SRR_ID" --max-size unlimited --progress -O "$DATA_DIR"
        PREFETCH_EXIT_CODE=$?

        # Check download success based on exit code and file existence
        if [[ $PREFETCH_EXIT_CODE -eq 0 ]] && [[ -f "$SRA_FILE" ]]; then
            echo "✅ Success: $SRR_ID"
        elif [[ -f "$SRA_FILE" ]]; then
            echo "⚠️  Warning: File exists but prefetch returned error code $PREFETCH_EXIT_CODE"
            echo "✅ Treating as success: $SRR_ID"
        else
            echo "❌ Failed: $SRR_ID (exit code: $PREFETCH_EXIT_CODE)"
            echo "   Cannot proceed without SRA file"
            exit 1
        fi
    fi
    
    echo "✅ Step 0 Complete: SRA file ready"
    echo ""
fi

# =============================
# Step 1: Convert SRA to FASTQ.gz
# =============================
echo "=========================================="
echo "📌 Step 1: SRA to FASTQ.gz"
echo "=========================================="

# Verify SRA file exists
if [[ ! -f "$SRA_FILE" ]]; then
    echo "❌ SRA file not found: $SRA_FILE"
    echo "   Please download the data first or use download_srr=True"
    exit 1
fi

echo "📁 SRA file: $SRA_FILE"

# Check if FASTQ conversion was already completed (skip to save time)
if [[ -f "${SRR_DIR}/${SRR_ID}_S1_R1_001.fastq.gz" ]]; then
    echo "⏭️  FASTQ.gz files already exist, skipping fastq-dump..."
else
    # Change to SRR directory for fastq-dump output
    cd "$SRR_DIR" || { echo "❌ Failed to cd into $SRR_DIR"; exit 1; }

    echo "🔄 Running fastq-dump with --gzip on $SRR_ID..."

    # Run fastq-dump with --split-files and --gzip to generate compressed FASTQ files
    if ! fastq-dump --split-files --gzip "${SRR_ID}.sra"; then
        echo "❌ fastq-dump failed for $SRR_ID"
        exit 1
    fi

    # Count how many files were generated
    NUM_FILES=$(ls -1 "${SRR_ID}_"*.fastq.gz 2>/dev/null | wc -l)
    
    if [[ "$NUM_FILES" -eq 0 ]]; then
        echo "❌ Error: No FASTQ.gz files generated"
        exit 1
    fi
    
    echo "📊 Generated $NUM_FILES FASTQ.gz file(s)"
    
    # =============================
    # Intelligent Renaming Based on 10X Version
    # =============================
    echo "🔍 Detecting 10X Genomics version and renaming files..."
    
    # Function to check if a file contains sequences of length 26-28 bp (R1 indicator)
    check_read_length() {
        local file=$1
        # Extract first sequence from compressed file and check length
        local seq_length=$(zcat "$file" | head -n 2 | tail -n 1 | tr -d '\n' | wc -c)
        if [[ $seq_length -ge 26 && $seq_length -le 28 ]]; then
            return 0  # True - this is R1
        else
            return 1  # False - this is not R1
        fi
    }
    
    if [[ "$NUM_FILES" -eq 2 ]]; then
        # v3 Chemistry: 2 files
        echo "📌 Detected 10X v3 chemistry (2 files)"
        
        # Determine which file is R1 by checking sequence length
        if check_read_length "${SRR_ID}_1.fastq.gz"; then
            echo "   _1.fastq.gz is R1 (26-28 bp sequences detected)"
            R1_SOURCE="${SRR_ID}_1.fastq.gz"
            R2_SOURCE="${SRR_ID}_2.fastq.gz"
        else
            echo "   _2.fastq.gz is R1 (26-28 bp sequences detected)"
            R1_SOURCE="${SRR_ID}_2.fastq.gz"
            R2_SOURCE="${SRR_ID}_1.fastq.gz"
        fi
        
        # Rename files
        echo "📝 Renaming files with prefix: ${SRR_ID}_S1"
        [[ -f "$R1_SOURCE" ]] && mv "$R1_SOURCE" "${SRR_ID}_S1_R1_001.fastq.gz"
        [[ -f "$R2_SOURCE" ]] && mv "$R2_SOURCE" "${SRR_ID}_S1_R2_001.fastq.gz"
        
    elif [[ "$NUM_FILES" -eq 3 ]]; then
        # v2 Chemistry: 3 files
        echo "📌 Detected 10X v2 chemistry (3 files)"
        
        # For v2, typically:
        # _1.fastq.gz = I1 (Index - barcode + UMI)
        # _2.fastq.gz = R1 (Read 1 - should be 26-28 bp)
        # _3.fastq.gz = R2 (Read 2)
        
        # Verify _2 is R1 by checking sequence length
        if check_read_length "${SRR_ID}_2.fastq.gz"; then
            echo "   _2.fastq.gz confirmed as R1 (26-28 bp sequences detected)"
            echo "📝 Renaming files with prefix: ${SRR_ID}_S1"
            [[ -f "${SRR_ID}_1.fastq.gz" ]] && mv "${SRR_ID}_1.fastq.gz" "${SRR_ID}_S1_I1_001.fastq.gz"
            [[ -f "${SRR_ID}_2.fastq.gz" ]] && mv "${SRR_ID}_2.fastq.gz" "${SRR_ID}_S1_R1_001.fastq.gz"
            [[ -f "${SRR_ID}_3.fastq.gz" ]] && mv "${SRR_ID}_3.fastq.gz" "${SRR_ID}_S1_R2_001.fastq.gz"
        else
            echo "⚠️  Warning: _2.fastq.gz doesn't match expected R1 length (26-28 bp)"
            echo "   Proceeding with standard v2 naming convention"
            [[ -f "${SRR_ID}_1.fastq.gz" ]] && mv "${SRR_ID}_1.fastq.gz" "${SRR_ID}_S1_I1_001.fastq.gz"
            [[ -f "${SRR_ID}_2.fastq.gz" ]] && mv "${SRR_ID}_2.fastq.gz" "${SRR_ID}_S1_R1_001.fastq.gz"
            [[ -f "${SRR_ID}_3.fastq.gz" ]] && mv "${SRR_ID}_3.fastq.gz" "${SRR_ID}_S1_R2_001.fastq.gz"
        fi
        
    else
        echo "❌ Error: Unexpected number of files ($NUM_FILES)"
        echo "   Expected 2 (v3) or 3 (v2) files"
        exit 1
    fi
fi

# Verify final R1 file has correct sequence length
echo ""
echo "✅ Verifying R1 sequence length..."
R1_FILE="${SRR_DIR}/${SRR_ID}_S1_R1_001.fastq.gz"
if [[ -f "$R1_FILE" ]]; then
    SAMPLE_SEQ=$(zcat "$R1_FILE" | head -n 2 | tail -n 1)
    SEQ_LEN=${#SAMPLE_SEQ}
    echo "   Sample R1 sequence: $SAMPLE_SEQ"
    echo "   Length: $SEQ_LEN bp"
    
    if [[ $SEQ_LEN -ge 26 && $SEQ_LEN -le 28 ]]; then
        echo "   ✅ R1 length verification passed (26-28 bp)"
    else
        echo "   ⚠️  Warning: R1 length is $SEQ_LEN bp (expected 26-28 bp)"
    fi
fi

# Display generated file statistics
echo ""
echo "📊 FASTQ.gz files:"
ls -lh "${SRR_DIR}/${SRR_ID}_S1_"*.fastq.gz 2>/dev/null || echo "   No files found!"

echo "✅ Step 1 Complete: $SRR_ID"

# =============================
# Step 2: Split FASTQ by Individual
# =============================
echo ""
echo "=========================================="
echo "📌 Step 2: Split FASTQ by Individual"
echo "=========================================="

# Verify cellbarcodes file exists
if [[ ! -f "$cellbarcodes_file" ]]; then
    echo "❌ Cellbarcodes file not found: $cellbarcodes_file"
    exit 1
fi

echo "📁 Cellbarcodes file: $cellbarcodes_file"

# Path configuration for FASTQ splitting
FASTQ_DIR="${SRR_DIR}"
SPLIT_FASTQ_DIR="${FASTQ_DIR}/split_fastq"
FASTQ_PREFIX="${SRR_ID}_S1"

mkdir -p "$SPLIT_FASTQ_DIR"

# =============================
# Extract Unique Individual IDs
# =============================
echo ""
echo "🔍 Extracting unique individual IDs from cellbarcodes file..."

# Skip header, extract first column, remove special characters (_-*), sort and deduplicate
INDIVIDUALS=$(tail -n +2 "$cellbarcodes_file" | cut -d',' -f1 | sed 's/[_\-\*]//g' | sort -u)

echo "📋 Found individuals:"
for ind in $INDIVIDUALS; do
    echo "   - $ind"
done
echo ""

# =============================
# Generate reads_information.txt
# =============================
# This file maps read IDs to their barcode sequences for efficient filtering
echo "[$SRR_ID] Generating reads_information.txt from compressed FASTQ..."
zcat "${FASTQ_DIR}/${FASTQ_PREFIX}_R1_001.fastq.gz" | \
    awk 'NR%4==1 {read_id=$1} NR%4==2 {print read_id"\t"$1}' > \
    "$SPLIT_FASTQ_DIR/reads_information.txt"

TOTAL_READS=$(wc -l < "$SPLIT_FASTQ_DIR/reads_information.txt")
echo "   Total reads: $TOTAL_READS"

# =============================
# Process Each Individual
# =============================
for ind in $INDIVIDUALS; do
    echo ""
    echo "----------------------------------------"
    echo "[$SRR_ID][$ind] Processing individual: $ind"
    echo "----------------------------------------"
    
    # Extract barcodes for this individual from cellbarcodes file
    # Note: Original individual IDs may contain special characters, need to match original format
    # Using awk to match: rows where cleaned ID equals $ind
    echo "[$SRR_ID][$ind] Extracting cell barcodes..."
    tail -n +2 "$cellbarcodes_file" | awk -F',' -v ind="$ind" '
    {
        # Remove special characters from first column for matching
        cleaned = $1
        gsub(/[_\-\*]/, "", cleaned)
        if (cleaned == ind) {
            # Extract barcode (remove -1 suffix if present, 10X convention)
            split($2, arr, "-")
            print arr[1]
        }
    }' > "$SPLIT_FASTQ_DIR/cell_barcodes_only_${ind}.txt"
    
    BC_COUNT=$(wc -l < "$SPLIT_FASTQ_DIR/cell_barcodes_only_${ind}.txt")
    echo "   Found $BC_COUNT barcodes for individual $ind"
    
    if [[ "$BC_COUNT" -eq 0 ]]; then
        echo "   ⚠️ No barcodes found, skipping..."
        continue
    fi

    # Match reads to barcodes
    # Compare first 16bp of read sequence (cell barcode) against barcode list
    echo "[$SRR_ID][$ind] Matching reads to barcodes..."
    awk 'NR==FNR {bc[$1]; next} {cell_bc=substr($2, 1, 16); if (cell_bc in bc) print $0}' \
        "$SPLIT_FASTQ_DIR/cell_barcodes_only_${ind}.txt" \
        "$SPLIT_FASTQ_DIR/reads_information.txt" > \
        "$SPLIT_FASTQ_DIR/reads_information_${ind}.txt"

    MATCHED_READS=$(wc -l < "$SPLIT_FASTQ_DIR/reads_information_${ind}.txt")
    echo "   Matched reads: $MATCHED_READS"

    if [[ "$MATCHED_READS" -eq 0 ]]; then
        echo "   ⚠️ No reads matched, skipping..."
        continue
    fi

    # Extract read IDs for FASTQ filtering (remove @ prefix)
    cut -f1 "$SPLIT_FASTQ_DIR/reads_information_${ind}.txt" | sed 's/^@//' > "$SPLIT_FASTQ_DIR/read_ids_${ind}.txt"

    # Split FASTQ files using BBMap's filterbyname.sh
    # Process all three read types: R1 (read 1), R2 (read 2), I1 (index)
    echo "[$SRR_ID][$ind] Splitting FASTQ.gz files (R1, R2, I1)..."
    for read_type in R1 R2; do
        INPUT_FASTQ="${FASTQ_DIR}/${FASTQ_PREFIX}_${read_type}_001.fastq.gz"
        OUTPUT_FASTQ="${SPLIT_FASTQ_DIR}/${SRR_ID}_Ind${ind}_${read_type}_001.fastq.gz"
        
        if [[ ! -f "$INPUT_FASTQ" ]]; then
            echo "   ⚠️ Input file not found: $INPUT_FASTQ"
            continue
        fi
        
        # Use BBMap filterbyname.sh to extract reads by ID
        # include=t means keep reads that match the name list
        # BBMap can handle .gz files directly
        filterbyname.sh \
            in="$INPUT_FASTQ" \
            out="$OUTPUT_FASTQ" \
            names="$SPLIT_FASTQ_DIR/read_ids_${ind}.txt" \
            include=t 2>/dev/null
        
        if [[ -f "$OUTPUT_FASTQ" ]]; then
            # Count reads in compressed file (divide by 4 since each read has 4 lines)
            OUT_READS=$(($(zcat "$OUTPUT_FASTQ" | wc -l) / 4))
            echo "   ✓ ${read_type}: $OUT_READS reads"
        fi
    done
done

# =============================
# Cleanup Temporary Files (Optional)
# =============================
# Uncomment the following lines to remove intermediate files and save disk space
# rm -f "$SPLIT_FASTQ_DIR"/reads_information*.txt
# rm -f "$SPLIT_FASTQ_DIR"/cell_barcodes_only_*.txt
# rm -f "$SPLIT_FASTQ_DIR"/read_ids_*.txt

# =============================
# Final Summary
# =============================
echo ""
echo "=========================================="
echo "📊 Final Summary for $SRR_ID"
echo "=========================================="
echo "Split FASTQ files in: $SPLIT_FASTQ_DIR"
ls -lh "$SPLIT_FASTQ_DIR"/${SRR_ID}_Ind*_*.fastq.gz 2>/dev/null || echo "   No split files found!"

echo ""
echo "=========================================="
echo "✅ All steps complete for $SRR_ID"
echo "=========================================="