#!/bin/bash

set -euo pipefail

#########################################
# Parse named arguments
#########################################
for arg in "$@"; do
  case $arg in
    -root_path=*)
      ROOT_DIR="${arg#*=}"
      shift
      ;;
    -conda_cmd=*)
      CONDA_ACTIVATION_CMD="${arg#*=}"
      shift
      ;;
    -google_cloud_sdk=*)
      GOOGLE_SDK_CMD="${arg#*=}"
      shift
      ;;
    *)
      echo "❌ Unknown argument: $arg"
      echo "Usage: bash setup_reference.sh -root_path=PATH -conda_cmd=true|\"source*/...\" -google_cloud_sdk=true|\"export*/...\""
      exit 1
      ;;
  esac
done

if [ -z "${ROOT_DIR:-}" ] || [ -z "${CONDA_ACTIVATION_CMD:-}" ] || [ -z "${GOOGLE_SDK_CMD:-}" ]; then
  echo "❌ Missing required arguments."
  echo "Usage: bash setup_reference.sh -root_path=PATH -conda_cmd=... -google_cloud_sdk=..."
  exit 1
fi

echo "🔧 Base directory: $ROOT_DIR"
mkdir -p "$ROOT_DIR"
cd "$ROOT_DIR" || { echo "❌ Failed to enter $ROOT_DIR"; exit 1; }


#########################################
# Safe downloading function
#########################################
download_if_missing() {
  local url=$1
  local output=$2
  if [ ! -f "$output" ]; then
    echo "📥 Downloading $output ..."
    wget -O "$output" "$url"
  else
    echo "✅ $output already exists. Skipping."
  fi
}

#########################################
# Install or activate Google Cloud SDK
#########################################
if [ "$GOOGLE_SDK_CMD" == "true" ]; then
  echo "📦 Installing Google Cloud SDK under $ROOT_DIR/google-cloud-sdk"

  SDK_VERSION="456.0.0"
  SDK_FILENAME="google-cloud-sdk-${SDK_VERSION}-linux-x86_64.tar.gz"
  SDK_URL="https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/${SDK_FILENAME}"

  download_if_missing "$SDK_URL" "$SDK_FILENAME"
  tar -xzf "$SDK_FILENAME"
  ./google-cloud-sdk/install.sh --quiet

  echo "✅ Google Cloud SDK installed."
  export PATH="$ROOT_DIR/google-cloud-sdk/bin:$PATH"
else
  CLEANED_SDK_CMD="${GOOGLE_SDK_CMD//\*/ }"
  echo "🧪 Activating Google SDK using: $CLEANED_SDK_CMD"
  eval "$CLEANED_SDK_CMD"
fi

#########################################
# Conda install or activation
#########################################
if [ "$CONDA_ACTIVATION_CMD" == "true" ]; then
  echo "📦 Installing Miniconda under $ROOT_DIR/miniconda3"

  INSTALL_DIR="$ROOT_DIR/miniconda3"
  MINICONDA_SCRIPT="$ROOT_DIR/miniconda.sh"
  mkdir -p "$INSTALL_DIR"

  echo "📥 Downloading Miniconda installer..."
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "$MINICONDA_SCRIPT"

  echo "⚙️ Installing Miniconda to $INSTALL_DIR..."
  bash "$MINICONDA_SCRIPT" -b -p "$INSTALL_DIR"

  echo "🔧 Initializing Conda..."
  source "$INSTALL_DIR/etc/profile.d/conda.sh"
  conda init bash
  source ~/.bashrc

  echo "✅ Conda installed and base environment activated."
else
  # 清除掉 * 符号（让它变成合法 shell 命令）
  CLEANED_CONDA_CMD="${CONDA_ACTIVATION_CMD/\*/ }"

  echo "🧪 Activating Conda using: $CLEANED_CONDA_CMD"
  eval "$CLEANED_CONDA_CMD"

  if ! command -v conda &> /dev/null; then
    echo "❌ Conda command not available after activation."
    exit 1
  fi
fi

conda --version && echo "✅ Conda is ready."

#########################################
# Create Conda environment
#########################################
echo "📦 Checking if conda environment 'scASE2' exists..."
if conda env list | grep -qE "^\s*scASE2\s"; then
  echo "⚠️ Conda environment 'scASE' already exists. Skipping creation."
else
  echo "📦 Creating conda environment 'scASE2' with bioinformatics tools..."
  conda create -n scASE2 -c bioconda -c conda-forge \
    samtools bamtools htslib bbmap sra-tools \
    bedtools bwa pysam picard openjdk=17 -y

  echo "✅ Conda environment 'scASE' created."
fi

#########################################
# Environment test
#########################################
echo "🔬 Testing environment activation..."
conda activate scASE2
echo "📊 Packages in scASE:"
conda list | head -n 10

echo "🎉 All setup completed successfully under $ROOT_DIR"
#########################################
# Download Cell Ranger Reference
#########################################
download_if_missing \
  "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz" \
  "refdata-gex-GRCh38-2020-A.tar.gz"
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz
picard CreateSequenceDictionary \
  R=./refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  O=./refdata-gex-GRCh38-2020-A/fasta/genome.dict
conda deactivate # in this session, we need to use conda to run 'picard'

#########################################
# Download GATK Reference Files
#########################################
VCF_DIR="$ROOT_DIR/gatk"
mkdir -p "$VCF_DIR"
cd "$VCF_DIR"

conda activate scASE2
echo "activate conda environment"

#VCF_FILES=(
# "Homo_sapiens_assembly38.known_indels.vcf.gz"
#  "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
#  "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
#)
#
#for file in "${VCF_FILES[@]}"; do
#  gsutil cp "gs://genomics-public-data/resources/broad/hg38/v0/$file" .
#  gsutil cp "gs://genomics-public-data/resources/broad/hg38/v0/${file}.tbi" .
   
   mv "$file" "$target_dir/resources_broad_hg38_v0_$file"
   mv "${file}.tbi" "$target_dir/resources_broad_hg38_v0_${file}.tbi"
#done

DBSNP_VCF="Homo_sapiens_assembly38.dbsnp138.vcf"
gsutil cp "gs://genomics-public-data/resources/broad/hg38/v0/$DBSNP_VCF" .
bgzip "$DBSNP_VCF"
tabix -p vcf "${DBSNP_VCF}.gz"
mv "${dbsnp_vcf}.gz" "$target_dir/resources_broad_hg38_v0_${dbsnp_vcf}.gz"
mv "${dbsnp_vcf}.gz.tbi" "$target_dir/resources_broad_hg38_v0_${dbsnp_vcf}.gz.tbi"

interval_file="wgs_calling_regions.hg38.interval_list"
gsutil cp "gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list" .
mv "$interval_file" "$target_dir/resources_broad_hg38_v0_$interval_file"

echo "✅ GATK reference files saved to: $VCF_DIR"
cd "$ROOT_DIR"
conda deactivate

#########################################
# Download Singularity container
#########################################
#echo "📦 Building or pulling Singularity container..."
#if ! command -v singularity &> /dev/null; then
#  echo "❌ Singularity not available."
#  exit 1
#fi

#if ! singularity build salsa_latest.sif docker://p4rkerw/salsa:latest; then
#  echo "⚠️ Build failed. Falling back to pull..."
#  singularity pull salsa_latest.sif docker://p4rkerw/salsa:latest
#fi
#echo "✅ Singularity container ready."

#########################################
# Download Funcotator data sources
#########################################
#FUNC_TAR="funcotator_dataSources.v1.6.20190124g.tar.gz"
#download_if_missing \
#  "https://storage.googleapis.com/broad-public-datasets/funcotator/$FUNC_TAR" \
#  "$FUNC_TAR"
#tar -xzf "$FUNC_TAR"
#echo "✅ Funcotator data extracted."

#########################################
# Download hg38 chromInfo
#########################################
#download_if_missing \
#  "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz" \
#  "$ROOT_DIR/hg38_chromInfo.txt.gz"

#########################################
# Download Phasing reference (chr1–22)
#########################################
#PHASING_DIR="$ROOT_DIR/phasing/biallelic_SNV"
#mkdir -p "$PHASING_DIR"
#cd "$PHASING_DIR"

#echo "📥 Downloading phasing reference files..."
#for chr in {1..22}; do
  #base="ALL.chr${chr}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
  #wget -c "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/${base}"
  #wget -c "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/${base}.tbi"
#done
#echo "✅ Phasing reference files saved to: $PHASING_DIR"

#########################################
# Download SHAPEIT5 binary
#########################################
#SHAPEIT5_DIR="$ROOT_DIR/shapeit5"
#mkdir -p "$SHAPEIT5_DIR"
#cd "$SHAPEIT5_DIR"

#SHAPEIT5_BIN="SHAPEIT5_phase_common_static_v1.0.0"
#SHAPEIT5_URL="https://github.com/odelaneau/shapeit5/releases/download/v1.0.0/$SHAPEIT5_BIN"

#if [ ! -f "$SHAPEIT5_BIN" ]; then
#  echo "📥 Downloading SHAPEIT5 binary..."
#  wget "$SHAPEIT5_URL"
#  chmod +x "$SHAPEIT5_BIN"
#  echo "✅ SHAPEIT5 downloaded and made executable."
#else
#  echo "✅ SHAPEIT5 binary already exists. Skipping download."
#fi


#echo "🎉 ALL SETUP COMPLETED SUCCESSFULLY"

