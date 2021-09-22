#!/bin/bash
#SBATCH --job-name=4_post-qc
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture5/out/4_post-qc.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture5/err/4_post-qc.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=186G
#SBATCH --time=30:00

module purge
module add apps/python/3.8.5

source main.env

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