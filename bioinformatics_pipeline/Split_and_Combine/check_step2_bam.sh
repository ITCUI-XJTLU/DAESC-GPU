#!/bin/bash

# 定义基目录
base_dir="/projects/gqilab/DAESC_GPU/data/OneK1K_Alldata/Sample2"

# 获取所有SRR开头的文件夹列表
srr_dirs=$(ls -d ${base_dir}/SRR*)

# 循环遍历每个SRR文件夹
for srr_dir in $srr_dirs; do
    # 提取SRR编号
    srr_id=$(basename $srr_dir)

    # 定义要计算BAM文件数量的路径
    bam_dir="${srr_dir}/cellranger_bam_${srr_id}/outs/filter_bam"

    # 检查路径是否存在
    if [ -d "$bam_dir" ]; then
        # 计算BAM文件的数量
        bam_count=$(ls -1 ${bam_dir}/*.bam 2>/dev/null | wc -l)

        # 输出结果
        echo "${srr_id} 下的 BAM 文件数量：${bam_count}"
    else
        echo "路径 ${bam_dir} 不存在"
    fi
done

