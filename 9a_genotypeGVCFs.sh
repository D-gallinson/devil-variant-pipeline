#!/bin/bash
#SBATCH --job-name=9a_GenotypeGVCFs-host_1
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/9a_genotypeGVCFs/host/out/1/%a.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/9a_genotypeGVCFs/host/err/1/%a.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-49

#########################################################################
# This script must point to directories indexed from 0 and corresponding
# to the batch number (e.g., 0=batch 1). These directories are used as
# multipliers to rename the err/out files properly upon completion
# 
# File names should then be changed using the following:
# 
# let batch_len=50
# ext=err
# for file in $(find -name *.err)
# do
#     file=${file:2}
#     let mult=$(echo $file | grep -E -o '^[0-9]')
#     let add=$(echo $file | grep -P -o '(?<=/)[0-9]')
#     let num=$mult*$batch_len+$add
#     mv $file final/$num.err
# done
#########################################################################



source main.env
source refs.env
source tools.env

sample_type="host"
first_batch="false"
progress=$SCRIPTS/logs/9a_genotypeGVCFs/$sample_type/progress.txt
let start=0

if [[ $first_batch == "true" ]]
then
    batch=1
    end=$SLURM_ARRAY_TASK_MAX
else
    let start=$(cut -f3 $progress | tail -1)+1
    let end=$SLURM_ARRAY_TASK_MAX+$start
    let batch=$(cut -f1 $progress | tail -1)+1
fi

if [[ $SLURM_ARRAY_TASK_ID == $SLURM_ARRAY_TASK_MAX ]]
then
    let final=$SLURM_ARRAY_TASK_MAX+$start
    echo -e "$batch\t$start\t$SLURM_ARRAY_TASK_MAX" >> $progress
fi

let file_num=$start+$SLURM_ARRAY_TASK_ID
interval_file="$(printf "%04d" $file_num)-scattered.interval_list"
output=$DATA/joint-variants
geno=${file_num}-${sample_type}.vcf
intervals=${INTERVAL_ARRAY}/${interval_file}

$GATK --java-options "-Xmx10G -Xms2G" GenotypeGVCFs \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir $output/tmp \
    -V gendb://$output/genomicsDB-${sample_type} \
    -R $REF \
    -L $intervals \
    -O $output/parallel/$geno