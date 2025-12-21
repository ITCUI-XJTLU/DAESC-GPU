#!/bin/bash
#SBATCH --job-name=daescp2
#SBATCH --output=logs/daescp2_%A_%a.out
#SBATCH --error=logs/daescp2_%A_%a.out
#SBATCH --array=1-2%2          
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=360G
#SBATCH --time=12:00:00    
 
# =========================
# Basic Configuration (Modify these paths according to your project setup)
# =========================
ROOT_PATH="/scratch/user/tfcui2025/DAESC/DAESC-P"
 
ml purge
ml GCC/12.3.0 OpenMPI/4.1.5
module load Singularity/3.10.2 

# =========================
# Execute Processing Script (Do not modify below this line)
# =========================
TXT_INDIVIDUALS="${ROOT_PATH}/Lists/List_Individuals.txt"
LINE_NUMBER="${SLURM_ARRAY_TASK_ID}"
mkdir -p "${ROOT_PATH}/logs"
cd "${ROOT_PATH}" || exit 1
line=$(sed -n "${LINE_NUMBER}p" "$TXT_INDIVIDUALS")
Individual_ID=$(echo "$line" | awk '{print $1}')
CHR_STRING=$(echo "$line" | sed -E 's/^[^ ]+ \(([^)]+)\)/\1/')
bash ${ROOT_PATH}/scr/DAESC-P2.sh \
    root_path="${ROOT_PATH}" \
    individual_id="${Individual_ID}" \
    chr_list="${CHR_STRING}"