#!/bin/bash
#SBATCH --job-name=isec
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/isec.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/isec.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=186G
#SBATCH --time=6-23:00:00

variants=${WORK_BGFS}/outputs/intermediates/8_joint-variants

start=`date +%s`

${HOME}/tools/bcftools-1.12/./bcftools \
    isec \
    --threads 24 \
    ${variants}/FINAL_SNPs-host.vcf.gz \
    ${variants}/FINAL_SNPs-tumor.vcf.gz \
    -p ${variants}/isec

end=`date +%s`
printf "Execution time: $((end-start))s"