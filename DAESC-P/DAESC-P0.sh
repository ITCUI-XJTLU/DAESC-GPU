#!/bin/bash
# =============================================================================
# DAESC-P0.sh - Preprocessing Script (Run on login node)
# =============================================================================
# Purpose:  
#   1. Generate List_Individuals.txt with unique individual IDs
#   2. Generate cell barcode files (CSV and TXT) for each individual
#   3. Optionally download SRA files using prefetch
#
# Usage: 
#   bash DAESC-P0.sh root_path=... cellbarcodes_file=... [list_input_file=...] [download_srr=True/False]
# 
# Example:
#   bash DAESC-P0.sh \
#       root_path=/scratch/user/tfcui2025/DAESC/mypipeline_2 \
#       cellbarcodes_file=/scratch/user/tfcui2025/DAESC/mypipeline_2/cellbarcodes/Sample11_three.csv \
#       list_input_file=/path/to/my_srr_list.txt \
#       download_srr=True
#
# Input Requirements:
#   - cellbarcodes_file: CSV with columns (individual_id, barcode, ...)
#   - list_input_file: (optional) File containing SRR IDs for download
#
# Output:
#   - Lists/List_Individuals.txt: Individual IDs with chromosome list
#   - Lists/List_Input.txt: Copy of the input SRR list file
#   - cellbarcodes/{individual}.csv: Per-individual barcode CSV files
#   - cellbarcodes/{individual}.txt: Per-individual barcode TXT files
#   - Data_SRR/{SRR_ID}/: Downloaded SRA files (if download_srr=True)
#
# Package Requirements:
#  - SRA Toolkit (for prefetch command)
# =============================================================================

# =============================
# Argument Parsing
# =============================
root_path=""
cellbarcodes_file=""
list_input_file=""
download_srr="False"

for ARG in "$@"; do
    case $ARG in
        root_path=*) root_path="${ARG#*=}" ;;
        cellbarcodes_file=*) cellbarcodes_file="${ARG#*=}" ;;
        list_input_file=*) list_input_file="${ARG#*=}" ;;
        download_srr=*) download_srr="${ARG#*=}" ;;
        *) echo "❌ Unknown argument: $ARG"; exit 1 ;;
    esac
done

# =============================
# Argument Validation
# =============================
if [[ -z "$root_path" || -z "$cellbarcodes_file" ]]; then
    echo "Usage: bash $0 root_path=... cellbarcodes_file=... [list_input_file=...] [download_srr=True/False]"
    exit 1
fi

if [[ ! -f "$cellbarcodes_file" ]]; then
    echo "❌ Cellbarcodes file not found: $cellbarcodes_file"
    exit 1
fi

echo "=========================================="
echo "🚀 DAESC-P0: Preprocessing"
echo "=========================================="
echo "📁 Root path: $root_path"
echo "📁 Cellbarcodes file: $cellbarcodes_file"
if [[ -n "$list_input_file" ]]; then
    echo "📁 List input file: $list_input_file"
fi
echo "📥 Download SRR: $download_srr"
echo ""

# =============================
# Create Required Directories
# =============================
LISTS_DIR="${root_path}/Lists"
CELLBARCODES_DIR="${root_path}/cellbarcodes"
DATA_DIR="${root_path}/Data_SRR"

mkdir -p "$LISTS_DIR"
mkdir -p "$CELLBARCODES_DIR"
mkdir -p "$DATA_DIR"

# =============================
# Copy List_Input.txt if provided
# =============================
if [[ -n "$list_input_file" ]]; then
    if [[ ! -f "$list_input_file" ]]; then
        echo "❌ List input file not found: $list_input_file"
        exit 1
    fi
    
    echo "=========================================="
    echo "📋 Copying List_Input.txt"
    echo "=========================================="
    echo "📁 Source: $list_input_file"
    echo "📁 Destination: ${LISTS_DIR}/List_Input.txt"
    
    cp "$list_input_file" "${LISTS_DIR}/List_Input.txt"
    
    if [[ $? -eq 0 ]]; then
        echo "✅ List_Input.txt copied successfully"
    else
        echo "❌ Failed to copy List_Input.txt"
        exit 1
    fi
    echo ""
fi

# =============================
# Step 1: Extract Unique Individual IDs and Generate List_Individuals.txt
# =============================
echo "=========================================="
echo "📌 Step 1: Generate List_Individuals.txt"
echo "=========================================="

# Build chromosome list (1 2 3 ... 22) - used for downstream SALSA analysis
CHR_LIST=$(seq 1 22 | tr '\n' ' ' | sed 's/ *$//')

# Extract unique individual IDs (remove special characters like _-*)
echo "🔍 Extracting unique individual IDs..."
INDIVIDUALS=$(tail -n +2 "$cellbarcodes_file" | cut -d',' -f1 | sed 's/[^a-zA-Z0-9]//g' | sort -u)

LIST_INDIVIDUALS="${LISTS_DIR}/List_Individuals.txt"
> "$LIST_INDIVIDUALS"

echo "📋 Found individuals:"
for ind in $INDIVIDUALS; do
    echo "   - $ind"
    # Format: individual_id (chr1 chr2 ... chr22)
    echo "$ind ($CHR_LIST)" >> "$LIST_INDIVIDUALS"
done

echo ""
echo "✅ Saved to: $LIST_INDIVIDUALS"
echo ""

# =============================
# Step 2: Generate Per-Individual Cellbarcode Files
# =============================
echo "=========================================="
echo "📌 Step 2: Generate cellbarcodes files"
echo "=========================================="

# Get unique individual IDs from original file (preserve original format for matching)
ORIGINAL_INDIVIDUALS=$(tail -n +2 "$cellbarcodes_file" | cut -d',' -f1 | sort -u)

for original_ind in $ORIGINAL_INDIVIDUALS; do
    # Clean ID (remove special characters for output filename)
    cleaned_ind=$(echo "$original_ind" | sed 's/[^a-zA-Z0-9]//g')
    
    CSV_OUTPUT="${CELLBARCODES_DIR}/${cleaned_ind}.csv"
    TXT_OUTPUT="${CELLBARCODES_DIR}/${cleaned_ind}.txt"
    
    echo "[$cleaned_ind] Generating CSV and TXT files..."
    
    # Generate CSV file (format: barcode,orig.ident,celltype)
    # This format is required by SALSA for barcode filtering
    awk -F',' -v individual="$original_ind" '
        BEGIN {
            OFS=","
            print "barcode","orig.ident","celltype"
        }
        NR>1 && $1 == individual {
            # Remove -N suffix from barcode (10X Genomics convention, e.g., -1, -11, -123)
            barcode = $2
            sub(/-[0-9]+$/, "", barcode)
            print barcode, "pbmc_1k", 1
        }
    ' "$cellbarcodes_file" > "$CSV_OUTPUT"
    
    # Generate TXT file (format: CB:Z:barcode)
    # This format is used for read filtering in DAESC-P1
    awk -F',' -v individual="$original_ind" '
        NR>1 && $1 == individual {
            print "CB:Z:" $2
        }
    ' "$cellbarcodes_file" > "$TXT_OUTPUT"
    
    # Count lines for summary
    CSV_LINES=$(($(wc -l < "$CSV_OUTPUT") - 1))  # Subtract header line
    TXT_LINES=$(wc -l < "$TXT_OUTPUT")
    echo "   ✓ ${cleaned_ind}.csv: $CSV_LINES barcodes"
    echo "   ✓ ${cleaned_ind}.txt: $TXT_LINES barcodes"
done

echo ""
echo "✅ Cellbarcodes files saved to: $CELLBARCODES_DIR"
echo ""

# =============================
# Step 3: Download SRA Files (Optional)
# =============================
if [[ "$download_srr" == "True" || "$download_srr" == "true" || "$download_srr" == "TRUE" ]]; then
    echo "=========================================="
    echo "📌 Step 3: Download SRA files"
    echo "=========================================="
    
    TXT_INPUT="${LISTS_DIR}/List_Input.txt"
    
    if [[ ! -f "$TXT_INPUT" ]]; then
        echo "❌ List_Input.txt not found: $TXT_INPUT"
        echo "   Please create this file with SRR IDs (one per line)"
        echo "   Or provide list_input_file parameter"
        exit 1
    fi
    
    
    # Count total files to download
    TOTAL_LINES=$(wc -l < "$TXT_INPUT")
    echo "📊 Total files to download: $TOTAL_LINES"
    echo ""
    
    # Process each SRR ID
    COUNTER=0
    while read -r SRR_ID; do
        ((COUNTER++))
        
        # Skip empty lines
        [[ -z "$SRR_ID" ]] && continue
        
        # Remove whitespace characters
        SRR_ID=$(echo "$SRR_ID" | tr -d '[:space:]')
        
        SRR_DIR="${DATA_DIR}/${SRR_ID}"
        SRA_FILE="${SRR_DIR}/${SRR_ID}.sra"
        
        # Check if already downloaded (skip to save time)
        if [[ -f "$SRA_FILE" ]]; then
            echo "⏭️  [$COUNTER/$TOTAL_LINES] Skipping $SRR_ID (already exists)"
            continue
        fi
        
        # Use prefetch to download SRA file with progress display
        echo "📥 [$COUNTER/$TOTAL_LINES] Downloading $SRR_ID..."
        prefetch "$SRR_ID" --max-size unlimited --progress -O "$DATA_DIR"
        PREFETCH_EXIT_CODE=$?

        # Check download success based on exit code and file existence
        if [[ $PREFETCH_EXIT_CODE -eq 0 ]] && [[ -f "$SRA_FILE" ]]; then
            echo "✅ [$COUNTER/$TOTAL_LINES] Success: $SRR_ID"
        elif [[ -f "$SRA_FILE" ]]; then
            echo "⚠️  [$COUNTER/$TOTAL_LINES] Warning: File exists but prefetch returned error code $PREFETCH_EXIT_CODE"
            echo "✅ [$COUNTER/$TOTAL_LINES] Treating as success: $SRR_ID"
        else
            echo "❌ [$COUNTER/$TOTAL_LINES] Failed: $SRR_ID (exit code: $PREFETCH_EXIT_CODE)"
        fi
        echo ""
        
    done < "$TXT_INPUT"
    
    echo "✅ Download complete!"
    echo ""
fi

# =============================
# Final Summary
# =============================
echo "=========================================="
echo "📊 Final Summary"
echo "=========================================="
echo ""
echo "📁 Directory structure:"
echo "   ${root_path}/"
echo "   ├── Lists/"
echo "   │   ├── List_Input.txt (user provided)"
echo "   │   └── List_Individuals.txt (generated)"
echo "   ├── cellbarcodes/"

for ind in $INDIVIDUALS; do
    echo "   │   ├── ${ind}.csv"
    echo "   │   └── ${ind}.txt"
done

echo "   └── Data_SRR/"

if [[ "$download_srr" == "True" || "$download_srr" == "true" || "$download_srr" == "TRUE" ]]; then
    if [[ -f "${LISTS_DIR}/List_Input.txt" ]]; then
        while read -r SRR_ID; do
            [[ -z "$SRR_ID" ]] && continue
            SRR_ID=$(echo "$SRR_ID" | tr -d '[:space:]')
            if [[ -f "${DATA_DIR}/${SRR_ID}/${SRR_ID}.sra" ]]; then
                echo "       ├── ${SRR_ID}/ ✓"
            else
                echo "       ├── ${SRR_ID}/ ✗"
            fi
        done < "${LISTS_DIR}/List_Input.txt"
    fi
fi

echo ""
echo "=========================================="
echo "✅ DAESC-P0 Complete!"
echo "=========================================="
echo ""
echo "Next step: Run DAESC-P1 with SLURM"
echo "   sbatch DAESC-P1_slurm.sh"
echo ""