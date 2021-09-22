#!/bin/bash
#SBATCH --job-name=6_combine
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture5/out/6_combine.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture5/err/6_combine.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=186G
#SBATCH --time=20:00

module purge
module add apps/python/3.8.5

source main.env

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