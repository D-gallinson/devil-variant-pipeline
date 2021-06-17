#!/bin/bash
#SBATCH --job-name=3_filter
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/9_vcf/test_runs/out/3_filter.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/9_vcf/test_runs/err/3_filter.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=2:00:00

###############################################
# Extra filters (such as --maf) can be added
# to this script. Values here should also be
# changed contigent on the figures/descriptives
# generated from 1_stats.sh and based on
# testing various ATOMM run
###############################################

module purge
module add apps/vcftools

variant_mode=SNP
sample_mode=host
dir=${WORK_BGFS}/outputs/intermediates/8_joint-variants
input=${dir}/FINAL_SNPs-${sample_mode}.vcf
output=${dir}/FINAL_FILTER_${variant_mode}-${sample_mode}

vcftools \
    --vcf $input \
    --max-missing 1 \
    --out $output