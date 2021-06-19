#!/bin/bash
#SBATCH --job-name=1_stats
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/9_vcf/test_runs/out/1_stats.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/9_vcf/test_runs/err/1_stats.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=01:00:00

# TODO
# Add R script to generate plots/descriptives
# The idea of this script is that it generates the necessary metric files/figs for whatever vcf file 
# (e.g., raw_snps, filtered_snps) is used as input

module purge
module add apps/vcftools

variants=${WORK_BGFS}/outputs/intermediates/8_joint-variants
stats=${WORK_BGFS}/outputs/intermediates/8_joint-variants/vcftools_metrics
input=${variants}/FINAL_SNPs-host.vcf.gz
output=${stats}/SNPs-host

#Get INFO column metrics
vcftools \
    --gzvcf $input \
    --get-INFO DP \
    --get-INFO QD \
    --get-INFO FS \
    --get-INFO SOR \
    --get-INFO MQ \
    --get-INFO MQRankSum \
    --get-INFO ReadPosRankSum \
    --out $output

#Get sample level metrics
vcftools \
    --gzvcf $input \
    --depth \
    --out ${output}

vcftools \
    --gzvcf $input \
    --site-mean-depth \
    --out ${output}

vcftools \
    --gzvcf $input \
    --freq2 \
    --out ${output}

vcftools \
    --gzvcf $input \
    --site-quality \
    --out ${output}

vcftools \
    --gzvcf $input \
    --missing-indv \
    --out ${output}

vcftools \
    --gzvcf $input \
    --missing-site \
    --out ${output}

# R scripts to be added down here