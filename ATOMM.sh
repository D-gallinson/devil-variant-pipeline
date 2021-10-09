#!/bin/bash
#SBATCH --job-name=ATOMM-marginal
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/ATOMM/prelim/devil_survival/rand_seq_estimate.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/ATOMM/prelim/devil_survival/rand_seq_estimate.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00

module purge
module add apps/matlab

# This must be run from the ATOMM/ directory
matlab -nodisplay -nosplash -r "devil_marginal"

r=''
r+='host <- read.table("output/marginal_host.txt"); patho <- read.table("output/marginal_pathogen.txt"); '
r+='sink("output/MARGINAL_STATS.txt"); '
r+='cat(paste("Host NAs:", sum(is.na(host$V4)), "\n")); cat(paste("Pathogen NAs:", sum(is.na(patho$V4)), "\n\n")); '
r+='host <- host[complete.cases(host),]; patho <- patho[complete.cases(patho),]; '
r+='host_sig <- host[host$V4 < 0.05,]; patho_sig <- patho[patho$V4 < 0.05,]; '
r+='cat(paste("Host SNPs (p < 0.05):", dim(host_sig)[1], "\n")); '
r+='summary(host_sig$V4); '
r+='cat("\n"); '
r+='cat(paste("Pathogen SNPs (p < 0.05):", dim(patho_sig)[1], "\n")); '
r+='summary(patho_sig$V4); '
Rscript -e "$r"
