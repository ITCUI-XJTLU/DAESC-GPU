#!/bin/bash

# 定义你要处理的所有SRR_ID
Sample="Sample10"
SRR_IDS=("SRR18029348" "SRR18029349")

for SRR_ID in "${SRR_IDS[@]}"
do
    echo "Submitting job for $SRR_ID"
    sbatch --job-name=sp2_s10_${SRR_ID} \
           --partition=all-12c128g \
           --ntasks=1 \
           --cpus-per-task=1 \
           --mem-per-cpu=8G \
           --time=1-00:00:00 \
           --output=${SRR_ID}.out \
           --error=${SRR_ID}.err \
           step2_split_bam.sh $Sample $SRR_ID
done

# gqi-32c256g,all-12c128g
