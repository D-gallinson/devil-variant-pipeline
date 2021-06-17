#!/bin/bash
#SBATCH --job-name=4_genotype
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/9_vcf/test_runs/out/4_genotype.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/9_vcf/test_runs/err/4_genotype.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module purge
module load apps/vcftools
module load apps/python/3.8.5

sample_mode=host    # Use host/pathogen (tumor=pathogen)
variant_mode=SNP    # Use SNP/indel
type=${variant_mode}-${sample_mode}
dir=${WORK_BGFS}/outputs/intermediates/8_joint-variants
input=${dir}/FINAL-FILTER_${type}.recode.vcf

echo "Sample mode: ${sample_mode}\n"
echo "Variant mode: ${variant_mode}\n"

# Generate the 012 genotype matrix
vcftools \
    --vcf $input \
    --012 \
    --out ${dir}/geno_${type}

# Tranpose the matrix so that samples are columns and rows are variant sites (generate mode flag using first char of ${var_mode}${sample_mode})
python3 ${WORK_BGFS}/scripts/master/utility/mat_T.py ${var_mode:0:1}${sample_mode:0:1}

# Remove the first line (sample ID line)
tail -n +2 ${dir}/geno_${type}.T.012 > ${dir}/tmp.txt

# Remove the first line (sample ID line)
mv ${dir}/tmp.txt ${dir}/geno_${type}.T.012

# Replace 2's (homozygous alt) with 1's
sed 's/2/1/g' ${dir}/geno_${type}.T.012 > geno_${type}.T.01

# Prepend the SNP ID column (consecutive increasing numbers from 1..N)
wc -l geno_${type}.T.01 | 
cut -d ' ' -f1 | 
xargs seq | 
paste - geno_${type}.T.01 > geno_${type}.ID.T.01

# Add the chr ID from the .pos file. X, Y, and unplaced scaffolds are converted as follows:
# X > 7
# Y > 8
# unplaced > 9
sed 's/[A-Z]*[0-9]*\.[0-9]/9/g; s/X/7/g; s/Y/8/g' geno_${type}.012.pos | 
cut -f1 | 
paste - geno_${type}.ID.T.01 > sequence_${sample_mode}.txt

# Replace tab-delimited to space-delimited
sed -i -e 's/\t/ /g' ${dir}/sequence_${sample_mode}.txt

# Clean up the intermediate geno_ files
rm ${dir}/geno_*

# The following displays genotype information from the input/output file to be compared.
# This is done as a sanity check to ensure that the ATOMM genotype file aligns with
# the input VCF file. Only the first and last lines of the files are ouput, and only
# the first 20 genotypes from those lines.

# Display row 1, genotypes 1-20 from the input VCF file
grep -m 1 -P '^(?!#)' $input | 
grep -o -P '[01]/[01]' | 
column -x -c 4 | 
head -1 | 
cut -f1-20

# Display row 1, genotypes 1-20 from the final genotype file
head -1 ${dir}/sequence_${sample_mode}.txt | 
column -x -c 4 | 
cut -f3-22

# Display the final row, genotypes 1-20 from the input VCF file
tail -n 1 $input | 
grep -o -P '[01]/[01]' | 
column -x -c 4 | 
head -1 | 
cut -f1-20

# Display the final row, genotypes 1-20 from the final genotype file
tail -n 1 ${dir}/sequence_${sample_mode}.txt | 
column -x -c 4 | 
cut -f3-22