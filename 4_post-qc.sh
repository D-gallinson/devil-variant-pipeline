#!/bin/bash
#SBATCH --job-name=4_post-qc
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=scripts/master/logs/Capture1_6-11-21/out/4_post-qc.out
#SBATCH --error=scripts/master/logs/Capture1_6-11-21/err/4_post-qc.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=186G
#SBATCH --time=45:00

module purge
module add apps/python/3.8.5

source ${WORK_BGFS}/scripts/master/main.env

src=${DATA}/${batch}/3_trim
dest=${RESULTS}/${batch}/qc/post
out_dir=MultiQC

#mv fastqc and trimming reports from intermediates directory to results directory
mv ${src}/*fastqc* ${dest}/fastqc
mv ${src}/*trimming_report* ${dest}/trimming

#Multiqc fastqc and trimming reports
multiqc \
    -n postQC_multiqc \
    -o ${dest}/${out_dir} \
    ${dest}/fastqc \

multiqc \
    -n trimming_multiqc \
    -o ${dest}/${out_dir} \
    ${dest}/trimming \