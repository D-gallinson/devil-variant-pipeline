#!/bin/bash
#SBATCH --job-name=parallel_GenotypeGVCFs
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture1_6-11-21/out/genotypeGVCFs/%a_genotypeGVCFs.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture1_6-11-21/err/genotypeGVCFs/%a_genotypeGVCFs.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=14-00:00:00
#SBATCH --array=1-7
#SBATCH --nodelist=mdc-1057-4-7

source main.env
source refs.env
source tools.env

sample_type=$1
output=$DATA/joint-variants
geno=chrom_${SLURM_ARRAY_TASK_ID}-${sample_type}.vcf
intervals=${INTERVAL_ARRAY}/${SLURM_ARRAY_TASK_ID}.list

$GATK --java-options "-Xmx25G" GenotypeGVCFs \
    -V gendb://$output/genomicsDB-${sample_type} \
    -L  $intervals \
    -R $REF \
    -O $output/parallel/$geno