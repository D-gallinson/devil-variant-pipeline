#!/bin/bash
#SBATCH --job-name=ATOMM-convergence
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/scripts/master/logs/ATOMM/prelim/devil_survival/subset_filter_marginal.out
#SBATCH --error=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/scripts/master/logs/ATOMM/prelim/devil_survival/subset_filter_marginal.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00

module purge
module add apps/matlab

# This must be run from the ATOMM/ directory
# matlab -nodisplay -nosplash -r "devil_marginal"

r=''
r+='host <- read.table("output/marginal_host.txt"); patho <- read.table("output/marginal_pathogen.txt"); '
r+='sink("output/MARGINAL_STATS.txt"); '
r+='cat("Host NAs:", sum(is.na(host$V4)), "\n"); cat("Pathogen NAs:", sum(is.na(patho$V4)), "\n\n"); '
r+='host <- host[complete.cases(host),]; patho <- patho[complete.cases(patho),]; '
r+='cat("===== Host p-values =====\n"); cat("Total host SNPs:", nrow(host), "\n"); '
r+='summary(host$V4); cat("\n"); '
r+='host_sig <- host[host$V4 < 0.05,]; patho_sig <- patho[patho$V4 < 0.05,]; '
r+='cat("Host SNPs (p < 0.05): ", nrow(host_sig), " (", round((nrow(host_sig) / nrow(host)) * 100, 2), "%)", "\n", sep=""); '
r+='summary(host_sig$V4); '
r+='cat("\n"); '
r+='cat("===== Pathogen p-values =====\n"); cat("Total pathogen SNPs:", nrow(patho), "\n"); '
r+='summary(patho$V4); cat("\n"); '
r+='cat("Pathogen SNPs (p < 0.05):", nrow(patho_sig), " (", round((nrow(patho_sig) / nrow(patho)) * 100, 2), "%)", "\n", sep=""); '
r+='summary(patho_sig$V4); '
Rscript -e "$r"
