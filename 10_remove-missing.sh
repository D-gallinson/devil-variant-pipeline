#!/bin/bash
#SBATCH --job-name=remove-missing
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/10_remove-missing.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/10_remove-missing.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=2:00:00

module purge
module add apps/vcftools

variants=${WORK_BGFS}/outputs/intermediates/8_joint-variants
input=${variants}/FINAL_SNPs-host.vcf.gz
output=${variants}/SNPs-host_FINAL_FILTER

vcftools \
    --gzvcf $input \
    --max-missing 0.75 \
    --out $output