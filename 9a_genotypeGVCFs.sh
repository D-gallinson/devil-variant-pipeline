#!/bin/bash
#SBATCH --job-name=9a_GenotypeGVCFs-host_2-1
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/9a_genotypeGVCFs/host/out/2-1.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/9a_genotypeGVCFs/host/err/2-1.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=187G
#SBATCH --time=21-00:00:00

source main.env
source refs.env
source tools.env

sample_type="host"
interval_num=2
restore_point="1"

interval=${INTERVAL_ARRAY}/interval_list_chrom/${interval_num}.list
output=$DATA/joint-variants
geno=${interval_num}_${sample_type}.vcf

if [[ $restore_point ]]
then
    interval=$INTERVAL_ARRAY/restore_points/${interval_num}_${sample_type}.interval_list
    geno=${interval_num}-${restore_point}_${sample_type}.vcf
fi

$GATK --java-options "-Xmx150G -Xms150G -XX:ParallelGCThreads=2" GenotypeGVCFs \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir $output/tmp \
    -V gendb://$output/genomicsDB-${sample_type} \
    -R $REF \
    -L $interval \
    -O $output/parallel/$geno