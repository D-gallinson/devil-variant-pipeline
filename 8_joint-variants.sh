#!/bin/bash
#SBATCH --job-name=8_genotypeGVCFs
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture1_6-11-21/out/8_genotypeGVCFs.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture1_6-11-21/err/8_genotypeGVCFs.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=14-00:00:00

module purge
module add apps/python/3.8.5

source ${WORK_BGFS}/scripts/master/main.env

start=`date +%s`

#Run joint variant calling (GenomicsDBImport, GenotypeGVCFs, and hard filters)
mode=T #change to "H" for host run and "T" for tumor run
python3 ${SCRIPTS}/utility/gatk_second_loop.py $mode

end=`date +%s`
printf "Execution time: $((end-start))s"
