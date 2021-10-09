#!/bin/bash
#SBATCH --job-name=9a_GenotypeGVCFs-host_1
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/9a_genotypeGVCFs/host/out/3.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/9a_genotypeGVCFs/host/err/3.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=187G
#SBATCH --time=21-00:00:00

source main.env
source refs.env
source tools.env

sample_type="host"
interval_num=3

interval=$INTERVAL_ARRAY/${interval_num}.list
output=$DATA/joint-variants
geno=${interval_num}-${sample_type}.vcf

$GATK --java-options "-Xmx150G -Xms150G -XX:ParallelGCThreads=2" GenotypeGVCFs \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir $output/tmp \
    -V gendb://$output/genomicsDB-${sample_type} \
    -R $REF \
    -L $interval \
    -O $output/parallel/$geno
