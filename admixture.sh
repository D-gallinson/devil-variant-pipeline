#!/bin/bash
#SBATCH --job-name=admixture_%a
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/admixture/out/prelim/%a.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/admixture/err/prelim/%a.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --array=1-10

source main.env
source tools.env

input=$RESULTS/admixture/input/FINAL_SNPs-host.minDP_10.maxDP_100.alleles.missing_100.mac_2.chr.vcf
output=$RESULTS/admixture

$ADMIXTURE \
    --cv \
    ${input}.bed \
    $SLURM_ARRAY_TASK_ID \
    -j4

mv ${input}.${SLURM_ARRAY_TASK_ID}.* $output/$SLURM_ARRAY_TASK_ID