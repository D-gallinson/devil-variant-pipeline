#!/bin/bash
#SBATCH --job-name=2_pre-qc
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=scripts/master/logs/Capture1_6-11-21/out/2_pre-qc.out
#SBATCH --error=scripts/master/logs/Capture1_6-11-21/err/2_pre-qc.err
#SBATCH --ntasks=72
#SBATCH --nodes=3
#SBATCH --mem=186G
#SBATCH --time=04:30:00

module purge
module add apps/fastqc/0.11.5

source ${WORK_BGFS}/scripts/master/main.env

input=${DATA}/${batch}/reads_1/*.fastq.gz
output=${RESULTS}/${batch_name}/qc

#mv 1_cp.sh log files to proper logs/batch
mv ${LOGS}/1_cp.out ${LOGS}/${batch}/out
mv ${LOGS}/1_cp.err ${LOGS}/${batch}/err

#FastQC
fastqc \
    -t 48 \
    $input \
    -o ${output}/pre

#Necessary to prevent conflict between loading FastQC (which loads python2.7) and multiqc
module purge
module add apps/python/3.8.5

multiqc \
    -n preQC_multi \
    -o $output \
    ${output}/pre