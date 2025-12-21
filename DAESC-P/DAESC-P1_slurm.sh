#!/bin/bash
#SBATCH --job-name=daescp1
#SBATCH --output=logs/daescp1_%A_%a.out
#SBATCH --error=logs/daescp1_%A_%a.out
#SBATCH --array=1-2%2          # Adjust range based on List_Input.txt line count; %5 limits to 5 concurrent jobs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G               # Increased memory for FASTQ splitting operations
#SBATCH --time=10:00:00         # Increased time for large SRA files

# =========================
# Basic Configuration (Modify these paths according to your project setup)
# =========================
ROOT_PATH="/scratch/user/tfcui2025/DAESC/DAESC-P"
TXT_INPUT="${ROOT_PATH}/Lists/List_Input_OneK1K.txt"
CELLBARCODES_FILE="${ROOT_PATH}/cellbarcodes/Sample11_three.csv"

ml purge
ml GCC/12.3.0 OpenMPI/4.1.5
module load SRA-Toolkit/3.0.10   # Required for downloading SRA files
module load BBMap/39.19        # Required for FASTQ file processing 

# =========================
# Execute Processing Script (Do not modify below this line)
# =========================
LINE_NUMBER="${SLURM_ARRAY_TASK_ID}"
mkdir -p "${ROOT_PATH}/logs"
cd "${ROOT_PATH}" || exit 1
bash ${ROOT_PATH}/scr/DAESC-P1.sh \
    root_path="${ROOT_PATH}" \
    txt_input="${TXT_INPUT}" \
    line_number="${LINE_NUMBER}" \
    download_srr=False \
    cellbarcodes_file="${CELLBARCODES_FILE}"