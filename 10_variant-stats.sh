#!/bin/bash
#SBATCH --job-name=10_variant-stats
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/scripts/master/logs/ATOMM/prelim/devil_survival/vcf_stats_FINAL_SNPs-tumor.subset.out
#SBATCH --error=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/scripts/master/logs/ATOMM/prelim/devil_survival/vcf_stats_FINAL_SNPs-tumor.subset.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=01:00:00

module purge
module add apps/vcftools
module add compilers/gcc

set -euo pipefail

source main.env
source refs.env
source tools.env

sample_type="tumor" #Change to "tumor" or "host"
input_dir="${RESULTS}/filter/ATOMM"    #Change to relevant input folder
input_file="FINAL_SNPs-tumor.isec_raw.subset.missing_100.vcf.gz" #Change to relevant input file
output_folder="ATOMM/subset"  #Change to relevant output folder (assumed to be in $$RESULTS/vcf_stats)

fname=$(echo $output_folder | sed 's/.*\///g')
id=${sample_type}_${fname}
input="${input_dir}/${input_file}"
output="${RESULTS}/vcf_stats/${sample_type}/${output_folder}"

# Check to ensure that the VCF has a non-empty INFO column
info_check=$(zgrep -P -m1 '^(?!#)' $input | cut -f8)
if [[ $info_check == "." ]]
then
    echo "!!!!! INFO COLUMN IS EMPTY !!!!!"
    echo "Please re-run VCFTools with the \"--recode-INFO-all flag\""
    echo "Exiting"
    exit 1
fi

# Generate output folders
mkdir -p $output/data
mkdir -p $output/tables
mkdir -p $output/plots

# Generate INFO column metrics
vcftools \
    --gzvcf $input \
    --get-INFO DP \
    --get-INFO QD \
    --get-INFO FS \
    --get-INFO SOR \
    --get-INFO MQ \
    --get-INFO MQRankSum \
    --get-INFO ReadPosRankSum \
    --out $output/data/$id

# Generate tmp sample level metrics
vcftools \
    --gzvcf $input \
    --depth \
    --out $output/tmp

vcftools \
    --gzvcf $input \
    --missing-indv \
    --out $output/tmp

# Generate tmp site level metrics
vcftools \
    --gzvcf $input \
    --site-mean-depth \
    --out $output/tmp

vcftools \
    --gzvcf $input \
    --freq2 \
    --out $output/tmp

vcftools \
    --gzvcf $input \
    --counts2 \
    --out $output/tmp

vcftools \
    --gzvcf $input \
    --site-quality \
    --out $output/tmp

vcftools \
    --gzvcf $input \
    --missing-site \
    --out $output/tmp

# Generate SNP density
vcftools \
    --gzvcf $input \
    --SNPdensity 125000 \
    --out $output/data/$id

# Generate tmp pos file
$BGZIP -c -d $input | grep -P '^(?!##)' | cut -f1-2 | sed 's/#//' > $output/tmp.pos

# Generate the sample-level metrics file
cut -f 2,4,5 $output/tmp.imiss > $output/tmp.cut.imiss
paste $output/tmp.idepth $output/tmp.cut.imiss > $output/data/$id.sample

# Generate the site-level metrics file
echo -e "REF_FREQ\tALT_FREQ" > $output/tmp.cut.frq
echo -e "REF_COUNT\tALT_COUNT" > $output/tmp.cut.count
cut -f 5,6 $output/tmp.frq | tail -n +2 >> $output/tmp.cut.frq
cut -f 5,6 $output/tmp.frq.count | tail -n +2 >> $output/tmp.cut.count
cut -f 3,4 $output/tmp.ldepth.mean > $output/tmp.cut.ldepth.mean
cut -f 6 $output/tmp.lmiss > $output/tmp.cut.lmiss
paste $output/tmp.lqual $output/tmp.cut.frq $output/tmp.cut.count $output/tmp.cut.ldepth.mean $output/tmp.cut.lmiss > $output/data/$id.site

# Obtain SNPs chr and pos for karyotype density
grep -E '^[1-6X]' $CHROMS > $output/tmp.chroms
echo -e "chr\tstart\tend" > $output/tmp.grange.density
$BGZIP -c -d $input | cut -f1-2 | grep -P '^[1-6X]' > $output/tmp.del
paste <(cut -f1-2 $output/tmp.del) <(cut -f2 $output/tmp.del) >> $output/tmp.grange.density

# Generate basic metric files
Rscript utility/vcf_metrics.R $id $output

# Generate plots
Rscript utility/vcf_plot.R $id $output

# Generate the README file describing each file made from this script
echo "===== SUMMARY OF FILES GENERATED =====
data/
    $fname.INFO     Values from the INFO column of a VCF
    $fname.sample   Per-sample VCF statistics
    $fname.site     Per-site VCF statistics
    $fname.snpden   Number of SNPs within 125kb intervals spanning the genome
plots/
    info_plots.$fname.pdf           Plots generated from the INFO column of a VCF    
    karyotype_density.$fname.pdf    Density plots over the devil autosomes
    sample_plots.$fname.pdf         Plots generated from statistics generated per-sample
    site_plots.$fname.pdf           Plots generated from statistics generated per-site
    snp_density.$fname.pdf          A histogram of SNPs contained within 125kb intervals
tables/
    all_chroms.bp.$fname.snpden      Statistics regarding SNP distances and the number of SNPs per chromosome (SNP_count). Generated over all chromosomes
    big_chroms.bp.$fname.snpden      Statistics regarding SNP distances and the number of SNPs per chromosome (SNP_count). Generated over large chromosomes (1-6 and X)
    bottom_mac.$fname.csv            Table containing the number of sites with a minor allele count (MAC) of 0-10
    info_summary.$fname.csv          Statistics generated from the VCF INFO column
    little_chroms.bp.$fname.snpden   Statistics regarding SNP distances and the number of SNPs per chromosome (SNP_count). Generated over small chromosomes (Y and unplaced scaffolds)
    per_chrom.bp.$fname.csv          Mean SNP distances and the number of SNPs per chromosome (SNP_count). Generated for all chromosomes on a per-chromosome basis
    sample_summary.$fname.csv        Statistics generated per-sample
    site_summary.$fname.csv          Statistics generated per-site
    SNP_density.$fname.csv           Counts of the 125kb intervals containing n SNPs
" > $output/README.txt

# Remove tmp files
rm $output/tmp*