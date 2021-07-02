#!/bin/bash
#SBATCH --job-name=6_combine
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=scripts/master/logs/Capture1_6-11-21/out/6_combine.out
#SBATCH --error=scripts/master/logs/Capture1_6-11-21/err/6_combine.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=186G
#SBATCH --time=20:00

module purge
module add apps/python/3.8.5

source ${WORK_BGFS}/scripts/master/main.env

align=${RESULTS}/${batch}/align
out_dir=MultiQC

#Generate HSMetrics summary file (my custom python script)
python3 ${SCRIPTS}/utility/combine-hs.py ${align}/HS-metrics ${align}/HS_summary.csv

#MultiQC flagstat, duplicates stats, and HsMetrics
multiqc \
    -n flagstat_multiqc \
    -o ${align}/${out_dir} \
    ${align}/flagstat \

multiqc \
    -n duplicates_multiqc \
    -o ${align}/${out_dir} \
    ${align}/duplicates \

multiqc \
    -n HS_multiqc \
    -o ${align}/${out_dir} \
    ${align}/HS-metrics \