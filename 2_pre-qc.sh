#!/bin/bash
#SBATCH --job-name=2_pre-qc
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/out/2_pre-qc.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/err/2_pre-qc.err
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --mem=186G
#SBATCH --time=1-00:00:00

module purge
module add apps/fastqc/0.11.5

source main.env

input=${DATA}/${batch}/1_reads/*.fastq.gz
output=${RESULTS}/${batch}/qc

# mv 1_cp.sh log files to proper logs/batch
mv ${LOGS}/1_setup.out ${LOGS}/${batch}/out
mv ${LOGS}/1_setup.err ${LOGS}/${batch}/err

# FastQC
fastqc \
    -t 72 \
    $input \
    -o ${output}/pre

# Necessary to prevent conflict between loading FastQC (which loads python2.7) and multiqc
module purge
module add apps/python/3.8.5

multiqc \
    -n preQC_multi \
    -o $output \
    ${output}/pre