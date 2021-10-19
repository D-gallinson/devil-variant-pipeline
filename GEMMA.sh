#!/bin/bash
#SBATCH --job-name=GEMMA_survival
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/scripts/master/logs/GEMMA/out/prelim/xval_classic/%a_survival.out
#SBATCH --error=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/scripts/master/logs/GEMMA/err/prelim/xval_classic/%a_survival.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-4

module purge
module load apps/vcftools
module load apps/python/3.8.5
module load compilers/gcc

source main.env
source refs.env
source tools.env

set -euo pipefail

##########################################################################
# --- GENERAL DESCRIPTION ---
# Automatically run GEMMA, generating the pheno file and ordering
# the VCF samples according to the pheno file. If running in xval mode,
# then this should be run as an --array job where the number of jobs
# equals the number of sites to xval or iterations to xval
# 
# --- INTENDED USAGE ---
# This script is meant to achieve two goals:
# 1) Fit a GEMMA BSLMM to data in order to determine genomic architecture
# 2) Generate EBVs (i.e., phenotype predictions) using xval
# The workflow should thus be:
# utility/gemma_prep.R -> GEMMA.sh [mode=default] -> GEMMA.sh [mode=xval]
# Where "mode=xval" is either xval_site or xval_classic, dependent on the
# distribution of samples throughout each site. Running this script in
# mode=default will attain goal 1, generate necessary inputs for *both*
# xval modes, and generate GEMMA output files used for plotting. The
# subsequent runs of this script in either xval mode will attain goal 2,
# although it will be further necessary to determine the prediction
# accuracy (see utility/gemma_acc.py)
# 
# --- PHENOTYPE FILE FORMAT ---
# Microchip Pheno   Site
# chip_id1  val_1   site1
# chip_id2  NA      site2
# chip_id3  val_3   site1
# NOTE: val_N is numeric, and NA should be placed where values
#       are unknown (or where a prediction is intended to be made).
#       Optionally, site can be supplied if using the script in mode=xval
# 
# --- MODES ---
# default:      Fit a BSLMM to a phenotype using genotype data without 
#               making phenotype predictions. Best for determining genome
#               architecture
# xval_site:    Perform a leave-one-out blocked xval where each site is
#               left out on a given run. Not advised if some sites are 
#               highly sampled and others are sparse. Should be parallelized
#               by site (e.g., 5 sites means --array=0-4)
# xval_classic: A typical xval whereby the data are randomly split into 
#               a test and training set. The test percent must be 
#               divisible into 100 and this should be parallelized by
#               --array=0-(100/test_percent)-1. This ensures that all samples
#               are at one point used in the test set. The script 
#               "gemma_prep.R" should also be run beforehand to generate
#               the test_training_file
##########################################################################

# GEMMA settings
MIN_MAF=0.01            # GEMMA drops sites with MAF < MIN_MAF [default = 0.01]
MAX_MISS=0.05           # GEMMA drops sites with missingness > MAX_MISS and imputes sites with missingness between MAX and 1 using the mean genotype [default = 0.05]
MODEL=1                 # MCMC BSLMM model [default = 1]
BURNIN=1500000          # Burnin iterations (typically 10% of the chain) [default = 100,000]
CHAIN_LENGTH=15000000   # Chain length (i.e., number of iterations) [default = 1,000,000]
REL_MAT=1               # Relatedness matrix type (1 = centered, 2 = standardized) [default = 1]

# ============= PARAMETERS =============
# These should be changed to meet the needs of a specific run.
# DEFINITIONS
#    output_dir:         Directory to output GEMMA files to. If run in xval mode, site subdirs will be automatically made
#    out_prefix:         The prefix used for all output files. Only necessary for default mode runs
#    vcf_input:          The input VCF file. GEMMA does some basic filtering/imputation but this should be done beforehand
#    phenotype:          The phenotype input file, which should contain, in this order, cols: Microchip, my_pheno, Site (can accept NA values)
#    mode:               Can be "xval_site", "xval_classic" or "default"
#    sites:              Sites to use when running in xval_site mode (sites with no samples should not be included)
#    test_percent:       The percent of samples to be used in the test set in xval_classic mode (should be an integer evenly divisible into 100)
#    test_training_file: Path to the file containing each run's test split when in xval_classic mode (run utility/gemma_prep.R to generate)
output_dir=$DATA/results/GEMMA/survival
out_prefix=survival
vcf_input=$DATA/test/toy_host.vcf.gz
phenotype=$RESULTS/GEMMA/survival/phenotype_survival.txt
mode="xval_classic"
sites=("Arthur River" "Black River" "Freycinet" "Takone" "WPP")
test_percent=20
test_training_file=$RESULTS/GEMMA/survival/test_training.txt

# Intermediate files
template=$output_dir/tmp_template.vcf
bimbam_input=$output_dir/tmp_input.vcf
mean_geno=$output_dir/input.mg
geno_input=$output_dir/geno.mg


# The mean genotype file is only generated in default mode, meaning default mode must either be run
# before any xval or this file must be generated separately
if [[ $mode == "default" ]]
then
    # ============= ORDER VCF SAMPLES BY PHENOTYPE FILE SAMPLES =============
    # Generate a template VCF header with the sample order being determined by $phenotype.
    # The template VCF header is then used to shuffle the input VCF, ensuring the order
    # of samples in the VCF matches the order in the phenotype file. Example:
    # VCF:                  sample1    sample2    sample3
    # Phenotype:            sample2    sample1    sample3
    # VCF after shuffle:    sample2    sample1    sample3
    # NOTE: The $phenotype file is in list order (s1, s2, etc).T but is converted to a single
    #       tabular row by the script
    zgrep -m 1 '^##fileformat' $vcf_input > $template
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" $(cut -f 1 $phenotype | tail -n +2) | sed -r "s/\s+/\t/g" >> $template
    vcf-shuffle-cols -t $template $vcf_input > $bimbam_input

    # ===== QCTOOL BIMBAM GENOTYPE FILE GENERATION =====
    $HOME/tools/qctool/./qctool \
        -g $bimbam_input \
        -ofiletype bimbam_dosage \
        -og $mean_geno
fi

# If xval, make dirs for specific sites/test-sets to store gemma outputs
if [[ $mode == "xval_site" ]]
then
    site=${sites[$SLURM_ARRAY_TASK_ID]}
    out_prefix="${site/ /_}"
    output_dir=$output_dir/$out_prefix
    mkdir -p $output_dir
elif [[ $mode == "xval_classic" ]]
then
    let xval_iter=$SLURM_ARRAY_TASK_ID+1
    out_prefix="xval_${xval_iter}"
    output_dir=$output_dir/$out_prefix
    mkdir -p $output_dir   
else
    pheno_na=$phenotype
fi

# Convenience variables used during or after BSLMM fitting
full_out=$output_dir/$out_prefix
pheno_na=$output_dir/tmp_pheno_na.txt
pheno_input=$output_dir/tmp_pheno.txt

if [[ $mode == "xval_site" ]]
then
    # awk with delim=\t and replace col 2 with NA where col 3 matches the current site
    awk -v site="$site" 'BEGIN{FS=OFS="\t"} $3~site {$2="NA"}1' $phenotype > $pheno_na
elif [[ $mode == "xval_classic" ]]
then
    # This uses test_training_file to set phenotype values to NA for a given xval_classic array job
    r=''
    r+='args <- commandArgs(trailingOnly = T); '
    r+='na_df <- read.table(args[1], header = F, sep = "\t"); pheno_df <- read.table(args[2], header = F, sep = "\t", skip = 1); '
    r+='na_col <- na_df[, as.numeric(args[3])]; na_col <- na_col[!is.na(na_col)]; '
    r+='pheno_df[na_col, "V2"] <- NA; write.table(pheno_df, args[4], sep = "\t", row.names = F, quote = F); '
    Rscript -e "$r" $test_training_file $phenotype $xval_iter $pheno_na
fi

# Extract just the phenotype col (col 2), omitting the header
cut -f 2 $pheno_na | tail -n +2 > $pheno_input

# Delete this if I use allele_swap.py
geno_input=$mean_geno

# ===== SWAP REF/ALT TO MAJOR/MINOR =====
# The BIMBAM mean genotype format specifies that minor alleles=1 and
# major alleles=0, but VCF codes ref=0 and alt=1. I am uncertain if this
# discrepancy will affect the model output. This script takes the mean
# genotype file and a REF freq file and swaps 1s and 0s where REF<0.5.
# For exploratory purposes, the -a flag swaps all 0s and 1s (DO NOT USE
# FOR AN ACTUAL ANALYSIS).
# Generate freq file: VCFtools --freq2 | cut -f 5 | tail -n +2\

# python3 allele_swap.py -f input_ref.frq -o $geno_input $mean_geno
# python3 allele_swap.py -a -o $geno_input $mean_geno

# Generate a relatedness matrix. This is unnecessary for BSLMMM model fitting but
# is needed for predictions. As such, it is necessary to generate this in "xval"
# mode but not in "default" mode, although it is generated irrespective of the mode
$GEMMA \
    -gk 1 \
    -g $geno_input \
    -p $pheno_input \
    -outdir $output_dir \
    -o $out_prefix

# ===== RUN GEMMA =====
$GEMMA \
    -bslmm 1 \
    -g $geno_input \
    -p $pheno_input \
    -w $BURNIN \
    -s $CHAIN_LENGTH \
    -outdir $output_dir \
    -o $out_prefix

# The prefix.hyp.txt file produced by GEMMA has a trailing tab which messes with my R stats script
sed -i 's/\t$//g' $full_out.hyp.txt

# In either xval mode, make a phenotype prediction for each blocked xval
# (note that the relatedness matrix is either a cXX centered or sXX standardized file)
if [[ $mode != "default" ]]
then
    $GEMMA \
        -predict 1 \
        -g $geno_input \
        -p $pheno_input \
        -epm $full_out.param.txt \
        -emu $full_out.log.txt \
        -ebv $full_out.bv.txt \
        -k $full_out.*XX.txt \
        -outdir $output_dir \
        -o $out_prefix
fi

# Generate convergence and hyperparam stats when running in default mode
if [[ $mode == "default" ]]
then
    Rscript $SCRIPTS/utility/gemma_stats.R $full_out.hyp.txt $output_dir
fi

# rm all tmp files
rm $output_dir/tmp_*