#!/bin/bash
#SBATCH --job-name=10_variant-stats
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/outputs/results/vcf_stats/host/missing/run_log.out
#SBATCH --error=/work_bgfs/d/dgallinson/outputs/results/vcf_stats/host/missing/run_log.err
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

sample_type="host" #Change to "tumor" or "host"
input_dir="${RESULTS}/filter"    #Change to relevant input folder
input_file="FINAL_SNPs-host.minDP_10.maxDP_100.alleles.missing_100.mac_2.vcf.gz" #Change to relevant input file
output_folder="missing"  #Change to relevant output folder (assumed to be in $WORK_BGFS/outputs/results/vcf_stats)

id=${sample_type}_${output_folder}
input="${input_dir}/${input_file}"
output="${RESULTS}/vcf_stats/${sample_type}/${output_folder}"

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

# Remove tmp files
rm $output/tmp*
