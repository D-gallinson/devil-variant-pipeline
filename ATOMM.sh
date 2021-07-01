#!/bin/bash
#SBATCH --job-name=ATOMM-test
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/ATOMM-test.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/ATOMM-test.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=6-23:00:00

module purge
module add apps/matlab

# This must be run from the ATOMM/ directory
matlab -nodisplay -nosplash -r "devil_model"