#!/bin/bash
#SBATCH --job-name=GEMMA_YOB_blocked
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture2_7-29-21/out/GEMMA_YOB.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture2_7-29-21/err/GEMMA_YOB.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00

module purge
module load apps/vcftools
module load apps/python/3.8.5

source main.env
source refs.env
source tools.env

set -euo pipefail

##########################################################################
# Automatically run GEMMA, generating the pheno file and ordering
# the VCF samples according to the pheno file.
# 
# output_dir=output directory of choice
# input=filtered VCF generated at the end of the variant pipeline
# phenotype=phenotype file containing the chip IDs, numeric phenotype, and optional site (a header is expected)
# 
# PHENOTYPE FILE FORMAT
# chip_id1  val_1   site1
# chip_id2  NA      site2
# chip_id3  val_3   site1
# NOTE: val_N is numeric, and NA should be placed where values
#       are unknown (or where a prediction is intended to be made).
#       Optionally, site can be supplied if using the script in mode=xval
##########################################################################

# PARAMETERS: change these as necessary for a run
output_dir=$RESULTS/GEMMA/YOB
out_prefix=YOB
vcf_input=$DATA/joint-variants/preliminary/filtered/FINAL_SNPs-host.minDP_10.maxDP_100.alleles.missing_50.mac_2.vcf.gz
phenotype=../../phenotype_YOB.txt
mode="xval" #change to "xval" for a leave one out blocked xval or "default" for a standard run

# Intermediate files
template=$output_dir/tmp_template.vcf
bimbam_input=$output_dir/tmp_input.vcf
pheno_na=$output_dir/tmp_pheno_na.txt
pheno_input=$output_dir/tmp_pheno.txt
mean_geno=$output_dir/input.mg
geno_input=$output_dir/geno.mg

# ============= ORDER VCF SAMPLES BY PHENOTYPE FILE SAMPLES =============
# Generate a template VCF header with the sample order being determined by $phenotype.
# The template VCF header is then used to shuffle the input VCF, ensuring the order
# of samples in the VCF matches the order in the phenotype file. Example:
# VCF:                  sample1    sample2    sample3
# Phenotype:            sample2    sample1    sample3
# VCF after shuffle:    sample2    sample1    sample3
# NOTE: The $phenotype file is in list order (vertical s1, s2, etc) but is converted to a single
#       tabular row by the script
zgrep -m 1 '^##fileformat' $vcf_input > $template
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" $(cut -f 1 $phenotype | tail -n +2) | sed -r "s/\s+/\t/g" >> $template
vcf-shuffle-cols -t $template $vcf_input > $bimbam_input

# ===== QCTOOL BIMBAM GENOTYPE FILE GENERATION =====
$HOME/tools/qctool/./qctool \
    -g $bimbam_input \
    -ofiletype bimbam_dosage \
    -og $mean_geno

# In xval mode, GEMMA will be run 5 times, with leave on out xval by iterating through each site.
# Default mode will iterate once, with this hacky workaround ensuring a single iteration
# whereby the "NA" will either match nothing or be replaced with itself, thus using the unmodified
# phenotype file as input to GEMMA
if [[ $mode == "xval" ]]
then
    exclude_sites=("Arthur River" "Black River" "Freycinet" "Takone" "WPP")
else
    exclude_sites=("NA")
fi

# Processing loop (5 times for xval, once for default)
for ((i=0; i < ${#exclude_sites[@]}; i++))
do
    site="${exclude_sites[$i]}"
    gemma_out=$output_dir

    # If xval, make dirs for specific sites to store gemma outputs
    if [[ $mode == "xval" ]]
    then
        site_dir="${site/ /_}"
        mkdir -p $output_dir/$site_dir
        gemma_out=$output_dir/$site_dir
        out_prefix=$site_dir
    fi

    # awk with delim=\t and replace col 2 with NA where col 3 matches the current site
    awk -v site="$site" 'BEGIN{FS=OFS="\t"} $3~site {$2="NA"}1' $phenotype > $pheno_na
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
        -outdir $gemma_out \
        -o $out_prefix

    # ===== RUN GEMMA =====
    $GEMMA \
        -bslmm 1 \
        -g $geno_input \
        -p $pheno_input \
        -outdir $gemma_out \
        -o $out_prefix

    # In "xval" mode, make a phenotype prediction for each blocked xval
    # (note that the relatedness matrix is either a cXX centered or sXX standardized file)
    if [[ $mode == "xval" ]]
    then
        $GEMMA \
            -predict 1 \
            -g $geno_input \
            -p $pheno_input \
            -epm $out_prefix.param.txt \
            -emu $out_prefix.log.txt \
            -ebv $out_prefix.bv.txt \
            -k $out_prefix.*XX.txt \
            -outdir $gemma_out \
            -o $out_prefix
    fi

# rm all tmp files
rm $output_dir/tmp_*
done