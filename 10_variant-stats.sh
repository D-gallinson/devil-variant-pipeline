#!/bin/bash
#SBATCH --job-name=10_variant-stats_isec
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/results/joint-variants/final/GEMMA/tumor/survival/stats/log.out
#SBATCH --error=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/results/joint-variants/final/GEMMA/tumor/survival/stats/log.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module purge
module add apps/vcftools
module add compilers/gcc

set -euo pipefail

source main.env
source refs.env
source tools.env

# Variables to be modified between runs
# sample_type = "host" or "tumor"
# input       = location of input .vcf.gz
# outdir      = location of the output directory
# id_num      = optional unique ID to add to file names to disambiguate names between runs
sample_type="tumor"
input="${RESULTS}/joint-variants/final/GEMMA/tumor/survival/FINAL_SNPs-tumor.isec.survival.maf_05.missing_95.alleles.vcf.gz "
outdir="${RESULTS}/joint-variants/final/GEMMA/tumor/survival/stats"
id_num=""

# Constants
dir=$(basename $outdir)
id=${sample_type}_${dir}${id_num}

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
mkdir -p $outdir/data
mkdir -p $outdir/tables
mkdir -p $outdir/plots

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
    --out $outdir/data/$id

# Generate tmp sample level metrics
vcftools \
    --gzvcf $input \
    --depth \
    --out $outdir/tmp

vcftools \
    --gzvcf $input \
    --missing-indv \
    --out $outdir/tmp

# Generate tmp site level metrics
vcftools \
    --gzvcf $input \
    --site-mean-depth \
    --out $outdir/tmp

vcftools \
    --gzvcf $input \
    --freq2 \
    --out $outdir/tmp

vcftools \
    --gzvcf $input \
    --counts2 \
    --out $outdir/tmp

vcftools \
    --gzvcf $input \
    --site-quality \
    --out $outdir/tmp

vcftools \
    --gzvcf $input \
    --missing-site \
    --out $outdir/tmp

# Generate SNP density
vcftools \
    --gzvcf $input \
    --SNPdensity 125000 \
    --out $outdir/data/$id

# Generate tmp pos file
$BGZIP -c -d $input | grep -P '^(?!##)' | cut -f1-2 | sed 's/#//' > $outdir/tmp.pos

# Generate the sample-level metrics file
cut -f 2,4,5 $outdir/tmp.imiss > $outdir/tmp.cut.imiss
paste $outdir/tmp.idepth $outdir/tmp.cut.imiss > $outdir/data/$id.sample

# Generate the site-level metrics file
echo -e "REF_FREQ\tALT_FREQ" > $outdir/tmp.cut.frq
echo -e "REF_COUNT\tALT_COUNT" > $outdir/tmp.cut.count
cut -f 5,6 $outdir/tmp.frq | tail -n +2 >> $outdir/tmp.cut.frq
cut -f 5,6 $outdir/tmp.frq.count | tail -n +2 >> $outdir/tmp.cut.count
cut -f 3,4 $outdir/tmp.ldepth.mean > $outdir/tmp.cut.ldepth.mean
cut -f 6 $outdir/tmp.lmiss > $outdir/tmp.cut.lmiss
paste $outdir/tmp.lqual $outdir/tmp.cut.frq $outdir/tmp.cut.count $outdir/tmp.cut.ldepth.mean $outdir/tmp.cut.lmiss > $outdir/data/$id.site

# Obtain SNPs chr and pos for karyotype density
grep -E '^[1-6X]' $CHROMS > $outdir/tmp.chroms
echo -e "chr\tstart\tend" > $outdir/tmp.grange.density
$BGZIP -c -d $input | cut -f1-2 | grep -P '^[1-6X]' > $outdir/tmp.del
paste <(cut -f1-2 $outdir/tmp.del) <(cut -f2 $outdir/tmp.del) >> $outdir/tmp.grange.density

# Generate basic metric files
Rscript utility/vcf_metrics.R $id $outdir

# Generate plots
Rscript utility/vcf_plot.R $id $outdir

# Generate the README file describing each file made from this script
echo "===== SUMMARY OF FILES GENERATED =====
data/
    $id.INFO     Values from the INFO column of a VCF
    $id.sample   Per-sample VCF statistics
    $id.site     Per-site VCF statistics
    $id.snpden   Number of SNPs within 125kb intervals spanning the genome
plots/
    info_plots.$id.pdf           Plots generated from the INFO column of a VCF    
    karyotype_density.$id.pdf    Density plots over the devil autosomes
    sample_plots.$id.pdf         Plots generated from statistics generated per-sample
    site_plots.$id.pdf           Plots generated from statistics generated per-site
    snp_density.$id.pdf          A histogram of SNPs contained within 125kb intervals
tables/
    all_chroms.bp.$id.snpden      Statistics regarding SNP distances and the number of SNPs per chromosome (SNP_count). Generated over all chromosomes
    big_chroms.bp.$id.snpden      Statistics regarding SNP distances and the number of SNPs per chromosome (SNP_count). Generated over large chromosomes (1-6 and X)
    bottom_mac.$id.csv            Table containing the number of sites with a minor allele count (MAC) of 0-10
    info_summary.$id.csv          Statistics generated from the VCF INFO column
    little_chroms.bp.$id.snpden   Statistics regarding SNP distances and the number of SNPs per chromosome (SNP_count). Generated over small chromosomes (Y and unplaced scaffolds)
    per_chrom.bp.$id.csv          Mean SNP distances and the number of SNPs per chromosome (SNP_count). Generated for all chromosomes on a per-chromosome basis
    sample_summary.$id.csv        Statistics generated per-sample
    site_summary.$id.csv          Statistics generated per-site
    SNP_density.$id.csv           Counts of the 125kb intervals containing n SNPs
" > $outdir/README.txt

# Remove tmp files
rm $outdir/tmp*