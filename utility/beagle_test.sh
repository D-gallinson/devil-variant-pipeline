#!/bin/bash
#SBATCH --job-name=BEAGLE_acc
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=log.out
#SBATCH --error=log.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G

module purge
module load apps/vcftools
module load apps/python/3.8.5

source ${WORK_BGFS}/scripts/master/main.env

# Configurable
input_vcf="input.vcf"
py_script="${SCRIPTS}/utility/beagle_acc.py"
missing=0.5
iters=10

# Constants
TRUTH="tmp_truth"
MODEL="tmp_model"

echo "Settings"
echo "Input VCF: $input_vcf"
echo "--max-missing: $missing"
echo "Iterations: $iters"
echo ""

# Generate ${TRUTH}.vcf (ground truth file)
vcftools \
    --vcf $input_vcf \
    --max-missing 1 \
    --recode \
    --out $TRUTH

vcftools \
    --vcf ${TRUTH}.recode.vcf \
    --extract-FORMAT-info GT \
    --out ${TRUTH}

grep -P '^(?!##)' ${TRUTH}.recode.vcf | cut -f1-9 > start.vcf
cut -f3- ${TRUTH}.GT.FORMAT > end.vcf
paste start.vcf end.vcf > ${TRUTH}.vcf
sed -i 's/0\/1/1\/0/g' ${TRUTH}.vcf
rm start.vcf end.vcf ${TRUTH}.recode.vcf ${TRUTH}.GT.FORMAT

# Generate ${MODEL}.counts
vcftools \
    --vcf $input_vcf \
    --max-missing $missing \
    --recode \
    --out $MODEL

vcftools \
    --vcf ${MODEL}.recode.vcf \
    --counts \
    --out $MODEL

rm ${MODEL}.recode.vcf
cut -f4 ${MODEL}.frq.count > ${MODEL}.counts
rm ${MODEL}.frq.count

# Run python comparison script
python3 $py_script ${TRUTH}.vcf ${MODEL}.counts --iters $iters
rm tmp_*