#!/bin/bash
#SBATCH --job-name=9_joint-genotype_tumor
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/out/9b_joint-genotype.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/err/9b_joint_genotype.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00

set -euo pipefail

source main.env
source refs.env
source tools.env

sample_type="tumor"  # Set to "host" or "tumor"

output=$DATA/joint-variants
parallel=$output/parallel
geno_out=$output/output-${sample_type}.vcf

snp=$output/snps-${sample_type}.vcf
snp_filter=$output/snps-${sample_type}.filter.vcf
snp_final=$output/FINAL_SNPs-${sample_type}.vcf

indel=$output/indels-${sample_type}.vcf
indel_filter=$output/indels-${sample_type}.filter.vcf
indel_final=$output/FINAL_INDELs-${sample_type}.vcf

# Combine GenotypeGVCFs chrom outputs into a single VCF
head -1000 ${parallel}/0-${sample_type}.vcf | grep -E '^#' > $geno_out
n=($(ls $parallel/*-${sample_type}.vcf | wc -l))
for (( i=0; i<$n; i++ ))
do
    file=${parallel}/${i}-${sample_type}.vcf
    grep -P '^(?!#)' $file >> $geno_out
done

# Extract/filter SNPs
$GATK SelectVariants \
    -R $REF \
    -V $geno_out \
    --select-type-to-include SNP \
    -O $snp

$GATK VariantFiltration \
    -R $REF \
    -V $snp \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "FILTER" \
    -O $snp_filter

grep -E '^#|PASS' $snp_filter > $snp_final

# Extract/filter indels
$GATK SelectVariants \
    -R $REF \
    -V $geno_out \
    --select-type-to-include INDEL \
    -O $indel

$GATK VariantFiltration \
    -R $REF \
    -V $indel \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filter-name "FILTER" \
    -O $indel_filter

grep -E '^#|PASS' $indel_filter > $indel_final

# bgzip all vcf files
files=($geno_out $snp $snp_filter $snp_final $indel $indel_filter $indel_final)
for file in ${files[@]}
do
    $BGZIP -@ 8 $file
    $BCFTOOLS index --threads 8 $file.gz
done
