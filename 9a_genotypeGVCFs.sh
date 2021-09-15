#!/bin/bash
#SBATCH --job-name=9a_GenotypeGVCFs
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/out/9a_genotypeGVCFs/%a.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/err/9a_genotypeGVCFs/%a.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7800M
#SBATCH --time=4-00:00:00
#SBATCH --array=0-71

source main.env
source refs.env
source tools.env

sample_type="tumor"
interval_file="$(printf "%04d" $SLURM_ARRAY_TASK_ID)-scattered.interval_list"
output=$DATA/joint-variants
geno=${SLURM_ARRAY_TASK_ID}-${sample_type}.vcf
intervals=${INTERVAL_ARRAY}/${interval_file}

$GATK --java-options "-Xmx7700M -Xms2G" GenotypeGVCFs \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir $output/tmp \
    -V gendb://$output/genomicsDB-${sample_type} \
    -R $REF \
    -L $intervals \
    -O $output/parallel/$geno