#!/bin/bash
#SBATCH --job-name=admixture_demo
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/admixture/out/demo/%a.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/admixture/err/demo/%a.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --array=1-10

source main.env
source tools.env

module add compilers/gcc

input=$RESULTS/admixture/demo
file=FINAL_SNPs-host.minDP_10.maxDP_100.alleles.missing_100.mac_2.chr.vcf
file=hapmap3
output=$RESULTS/admixture/host_cv
log_dir="prelim"
plotting="true"
K=10

# PART 1: Run admixture in --cv mode (THIS SHOULD BE AN SBATCH)
if [[ $plotting != "true" ]]
then
    $ADMIXTURE \
        --cv \
        $input/${file}.bed \
        $SLURM_ARRAY_TASK_ID \
        -j4

    mv ${file}.${SLURM_ARRAY_TASK_ID}.* $output
# PART 2: plot the results of admixture (THIS SHOULD BE SOURCED)
else
    printf "CV\tK\n" > $output/CV.tsv
    for i in $(seq 1 $k)
    do
        grep "CV" $LOGS/admixture/out/$log_dir/$i.out | cut -d " " -f4 | xargs printf "%f\t$i\n" >> $output/CV.tsv
    done

    r=''
    r+='library(ggplot2); '
    r+='args <- commandArgs(trailingOnly = T); '
    r+='df <- read.table(paste(args[1], "/CV.tsv", sep = ""), header = T); '
    r+='plot <- ggplot(data = df, aes(x = K, y = CV)) + geom_line() + geom_point() + scale_x_continuous(breaks = seq(length(df$CV))) + ylab("Cross-validation error"); '
    r+='ggsave(paste(args[1], "/CV_plot.png", sep = ""), plot, height=1200, width=1200, units="px"); '
    Rscript -e "$r" $output
fi