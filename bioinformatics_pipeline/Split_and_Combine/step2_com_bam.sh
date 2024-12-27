#!/bin/bash
#SBATCH --job-name=st2_s5
#SBATCH --output=sample5_com.log
#SBATCH --error=sample5_com_err.log
#SBATCH --partition=all-12c128g
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00           # 作业运行的最大时间
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=16G

# 定义输入和输出目录的父路径
output_dir="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample10"
samples=("Sample10")
srrs=("SRR18029330" "SRR18029331" "SRR18029332" "SRR18029333" "SRR18029334" "SRR18029335" "SRR18029336" "SRR18029337" "SRR18029338" "SRR18029339" "SRR18029340" "SRR18029341" "SRR18029342" "SRR18029343" "SRR18029344" "SRR18029345" "SRR18029346" "SRR18029347" "SRR18029348" "SRR18029349")
individuals=("823_824" "824_825" "825_826" "826_827" "827_828" "829_830" "830_831" "831_832" "832_833" "842_843" "843_844")

# 循环处理每个 individual，生成对应的 BAM 文件
for individual in "${individuals[@]}"; do
    # 初始化 BAM 文件的路径列表
    bam_files_to_merge=""
    
    # 遍历样本和 SRR 列表，查找每个 BAM 文件
    for sample in "${samples[@]}"; do
        for srr in "${srrs[@]}"; do
            bam_path="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/${sample}/${srr}/cellranger_bam_${srr}/outs/filter_bam/individual_${individual}.bam"
            if [ -f "${bam_path}" ]; then
                # 统一读取组信息，确保所有文件的样本名相同
                fixed_bam_path="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/${sample}/${srr}/cellranger_bam_${srr}/outs/filter_bam/fixed_individual_${individual}.bam"
                samtools addreplacerg -r ID:${individual} -r SM:${individual} -r LB:${individual} -r PL:ILLUMINA -r PU:${individual} -o "${fixed_bam_path}" "${bam_path}"
                bam_files_to_merge="${bam_files_to_merge} ${fixed_bam_path}"
            else
                echo "Warning: ${bam_path} does not exist."
            fi
        done
    done
    
    # 定义输出文件夹，使用 individual 的名字
    output_folder="${output_dir}/com_${individual}"
    
    # 创建输出文件夹
    mkdir -p "${output_folder}"
    
    # 定义输出 BAM 文件路径，文件名为 merged_xxx_xxx.bam
    output_bam="${output_folder}/merged_${individual}.bam"
    
    # 合并 BAM 文件
    if [ -n "${bam_files_to_merge}" ]; then
        samtools merge -f "${output_bam}" ${bam_files_to_merge}
        echo "Merged individual_${individual}.bam into ${output_bam}"
        
        # 创建 BAM 文件索引
        samtools index "${output_bam}"
    else
        echo "No BAM files found for individual_${individual}."
    fi
done

