#!/bin/bash
#SBATCH --job-name=2_isec
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/9_vcf/test_runs/out/2_isec.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/9_vcf/test_runs/err/2_isec.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=186G
#SBATCH --time=6-23:00:00

dir=${WORK_BGFS}/outputs/intermediates/8_joint-variants
mkdir ${dir}/isec

start=`date +%s`

${HOME}/tools/bcftools-1.12/./bcftools \
    isec \
    --threads 24 \
    ${dir}/FINAL_SNPs-host.vcf.gz \
    ${dir}/FINAL_SNPs-tumor.vcf.gz \
    -p ${dir}/isec

end=`date +%s`
printf "Execution time: $((end-start))s"