#!/bin/bash
#SBATCH --job-name=ATOMM_prep
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/results/joint-variants/final/ATOMM/infection_age/filter.out
#SBATCH --error=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/results/joint-variants/final/ATOMM/infection_age/filter.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=10:00:00

module purge
module load apps/vcftools
module load apps/python/3.8.5

set -euo pipefail

source main.env
source refs.env
source tools.env

host=$RESULTS/joint-variants/final/FINAL_SNPs-host.vcf.gz
pathogen=$RESULTS/joint-variants/final/FINAL_SNPs-tumor.isec.vcf.gz
pheno_name=infection_age
rm_vcf_subset="true"

vcf_dir=$RESULTS/joint-variants/final/ATOMM/infection_age
ATOMM_dir=$RESULTS/ATOMM/input
log=$ATOMM_dir/atomm_prep.log
phenotype=$ATOMM_dir/master_phenotype.txt

host_tmp=$ATOMM_dir/tmp_host
host_subset=$vcf_dir/host.$pheno_name.subset.vcf.gz
host_final=$vcf_dir/host.$pheno_name.subset.missing_100.vcf.gz
pathogen_tmp=$ATOMM_dir/tmp_pathogen
pathogen_subset=$vcf_dir/pathogen.$pheno_name.subset.vcf.gz
pathogen_final=$vcf_dir/pathogen.$pheno_name.subset.missing_100.vcf.gz

# Generate VCF --keep files for filtering samples without phenotype data
cut -f1 $phenotype | tail -n +2 > $vcf_dir/host_list.keep
cut -f2 $phenotype | tail -n +2 > $vcf_dir/pathogen_list.keep

# Filter out missing phenotype samples
vcftools \
    --gzvcf $host \
    --keep $vcf_dir/host_list.keep \
    --recode \
    --recode-INFO-all \
    --stdout | 
    $BGZIP -@ 8 -c > $host_subset

vcftools \
    --gzvcf $pathogen \
    --keep $vcf_dir/pathogen_list.keep \
    --recode \
    --recode-INFO-all \
    --stdout | 
    $BGZIP -@ 8 -c > $pathogen_subset

# Filter SNPs with *any* missing genotypes or multiallelic sites
# (ATOMM can't handle missingness genotypes or anything beyond biallelic sites)
vcftools \
    --gzvcf $host_subset \
    --max-missing 1 \
    --max-alleles 2 \
    --min-alleles 2 \
    --recode \
    --recode-INFO-all \
    --stdout | 
    $BGZIP -@ 8 -c > $host_final

vcftools \
    --gzvcf $pathogen_subset \
    --max-missing 1 \
    --max-alleles 2 \
    --min-alleles 2 \
    --recode \
    --recode-INFO-all \
    --stdout | 
    $BGZIP -@ 8 -c > $pathogen_final

# Optionally delete the subsetted VCF files
if [[ $rm_vcf_subset == "true" ]]
then
    rm $host_subset $pathogen_subset
fi

# Sort the filtered VCFs based on the phenotype file order
zgrep -m 1 '^##fileformat' $host > ${ATOMM_dir}/host_template.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" $(cut -f 1 $vcf_dir/host_list.keep) | sed -r "s/\s+/\t/g" >> ${ATOMM_dir}/host_template.vcf
vcf-shuffle-cols -t ${ATOMM_dir}/host_template.vcf $host_final | $BGZIP -@ 8 -c > $host_tmp.vcf.gz

zgrep -m 1 '^##fileformat' $pathogen > ${ATOMM_dir}/pathogen_template.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" $(cut -f 1 $vcf_dir/pathogen_list.keep) | sed -r "s/\s+/\t/g" >> ${ATOMM_dir}/pathogen_template.vcf
vcf-shuffle-cols -t ${ATOMM_dir}/pathogen_template.vcf $pathogen_final | $BGZIP -@ 8 -c > $pathogen_tmp.vcf.gz

# Generate the 012 genotype matrix
vcftools \
    --gzvcf $host_tmp.vcf.gz \
    --012 \
    --out $host_tmp

vcftools \
    --gzvcf $pathogen_tmp.vcf.gz \
    --012 \
    --out $pathogen_tmp

# Cut out the sample IDs (first column) and tranpose the matrix so that samples are columns and rows are variant sites
cut -f2- $host_tmp.012 | datamash transpose > $host_tmp.T.012
cut -f2- $pathogen_tmp.012 | datamash transpose > $pathogen_tmp.T.012

# Replace 2's (homozygous alt) with 1's
sed 's/2/1/g' $host_tmp.T.012 > $host_tmp.T.01
sed 's/2/1/g' $pathogen_tmp.T.012 > $pathogen_tmp.T.01

# Prepend the SNP ID column (consecutive increasing numbers from 1..N)
wc -l $host_tmp.T.01 | 
cut -d ' ' -f1 | 
xargs seq | 
paste - $host_tmp.T.01 > $host_tmp.ID.T.01

wc -l $pathogen_tmp.T.01 | 
cut -d ' ' -f1 | 
xargs seq | 
paste - $pathogen_tmp.T.01 > $pathogen_tmp.ID.T.01

# Add the chr ID from the .pos file. X, Y, and unplaced scaffolds are converted as follows:
# X = 7
# Y = 8
# unplaced = 9
sed 's/[A-Z]*[0-9]*\.[0-9]/9/g; s/X/7/g; s/Y/8/g' $host_tmp.012.pos | 
cut -f1 | 
paste - $host_tmp.ID.T.01 > $ATOMM_dir/sequence_host.txt

sed 's/[A-Z]*[0-9]*\.[0-9]/9/g; s/X/7/g; s/Y/8/g' $pathogen_tmp.012.pos | 
cut -f1 | 
paste - $pathogen_tmp.ID.T.01 > $ATOMM_dir/sequence_pathogen.txt

# Replace tab-delimited to space-delimited
sed -i -e 's/\t/ /g' ${ATOMM_dir}/sequence_host.txt
sed -i -e 's/\t/ /g' ${ATOMM_dir}/sequence_pathogen.txt

# The following displays genotype information from the input/output file to be compared.
# This is done as a sanity check to ensure that the ATOMM genotype file aligns with
# the input VCF file. Only the first and last lines of the files are ouput, and only
# the first 20 genotypes from those lines.

# Configurable width to allocate each column (min should be 10)
width=15

# Format the first twenty 0/1 genotypes of the VCF with NA samples removed and the remaining samples shuffled
format_vcf (){
    if [[ $2 == "head" ]]
    then
        line=$(zgrep -m 1 -P '^(?!#)' $1)
    else
        line=$(zcat $1 | tail -1)
    fi
    echo $line |
    grep -o -P '[01]/[01]' |
    awk '{print}' ORS=' ' |
    cut -d' ' -f1-20 |
    xargs printf "%-${width}s" >> $log
    echo "" >> $log
}

# Format the first twenty 0s and 1s of the ATOMM genotype input
format_seq (){
    if [[ $2 == "head" ]]
    then
        line=$(head -1 ${ATOMM_dir}/sequence_${1}.txt)
    else
        line=$(tail -n 1 ${ATOMM_dir}/sequence_${1}.txt)
    fi
    echo $line |
    cut -d' ' -f3-22 |
    xargs printf "%-${width}s" >> $log
    echo "" >> $log
}

echo "===== DISPLAYING INPUT/OUTPUT GENOTYPES =====" > $log
echo -e "These are of the order:
SAMPLE_HEADER
HEAD -1 VCF
HEAD -1 ATOMM
TAIL -1 VCF
TAIL -1 ATOMM

Thus, the first and last two contiguous rows of the host and pathogen should align.
The two numbers in #/# are (where anything greater than 0 = 1):
0/0 = 0
0/1 = 1
1/1 = 1
" >> $log

# HOST
echo "HOST" >> $log

# Display Sample Names
zgrep -m 1 -P '^(?!##)' $host_tmp.vcf.gz |
cut -f10-29 |
xargs printf "%-${width}s" >> $log
echo "" >> $log

# Display row 1, genotypes 1-20 from the input VCF file
format_vcf $host_tmp.vcf.gz "head"

# Display row 1, genotypes 1-20 from the final genotype file
format_seq "host" "head"

# Display the final row, genotypes 1-20 from the input VCF file
format_vcf $host_tmp.vcf.gz "tail"

# Display the final row, genotypes 1-20 from the final genotype file
format_seq "host" "tail"

echo -e "\n" >> $log

# PATHOGEN
echo "PATHOGEN" >> $log

# Display Sample Names
zgrep -m 1 -P '^(?!##)' $pathogen_tmp.vcf.gz |
cut -f10-29 |
xargs printf "%-${width}s" >> $log
echo "" >> $log

# Display row 1, genotypes 1-20 from the input VCF file
format_vcf $pathogen_tmp.vcf.gz "head"

# Display row 1, genotypes 1-20 from the final genotype file
format_seq "pathogen" "head"

# Display the final row, genotypes 1-20 from the input VCF file
format_vcf $pathogen_tmp.vcf.gz "tail"

# Display the final row, genotypes 1-20 from the final genotype file
format_seq "pathogen" "tail"


# ===================== GENERATE PHENOTYPE FILE =====================
python - << EOF
import pandas as pd

pheno_file = "$phenotype"
outdir = "$ATOMM_dir"

df = pd.read_csv(pheno_file, sep="\t")
seq = [i for i in range(1, df.shape[0]+1)]
rep = [1 for i in range(1, df.shape[0]+1)]
prepend_df = pd.DataFrame({
    "host": seq,
    "patho": seq,
    "intercept": rep
})
atomm_pheno = pd.concat([prepend_df, df[df.columns[2:]]], axis=1)
atomm_pheno.to_csv(f"{outdir}/phenotype.txt", sep=" ", index=False, header=False)
EOF

# Clean up the intermediate geno_ files
rm ${ATOMM_dir}/tmp_*
rm $ATOMM_dir/*template*