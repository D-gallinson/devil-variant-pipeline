#!/bin/bash
#SBATCH --job-name=gvcf_bgzip
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/tmp/out/bgzip.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/tmp/err/bgzip.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00


source main.env
source tools.env

gvcf_dir=$DATA/gvcf

for file in $gvcf_dir/*.g.vcf
do
    $BGZIP $file
done