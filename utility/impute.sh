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

source main.env

set -euo pipefail

# Processing mode ["acc", "impute", "both"]
mode="both"

# Configurable
input_vcf="../../input.vcf"
out_dir="../../del"
output=$out_dir/"input"  #Used for actual BEAGLE imputation, change the fname as needed (omit an extension)
py_script="${SCRIPTS}/utility/beagle_acc.py"
missing=0 #Set this to 0 if the VCF should not be modified
iters=2

# Constants
TRUTH=$out_dir/"tmp_truth"
MODEL=$out_dir/"tmp_model"

# Ensures that the BEAGLE imputation accuracy represents the
# accuracy obtained when imputing on the unmodified input VCF
if [[ $mode == "both" ]]
then
    missing=0
fi

echo "Settings"
echo "Mode: $mode"
echo "Input VCF: $input_vcf"
echo "--max-missing: $missing"
echo "Iterations: $iters"
echo ""

if [[ $mode != "impute" ]]
then
    # Generate ${TRUTH}.vcf (ground truth file)
    vcftools \
        --vcf $input_vcf \
        --max-missing 1 \
        --recode \
        --out $TRUTH

    vcftools \
        --vcf ${TRUTH}.recode.vcf \
        --extract-FORMAT-info GT \
        --out $TRUTH

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
    python3 $py_script ${TRUTH}.vcf ${MODEL}.counts --iters $iters --output $out_dir/impute_accuracy

    # Clean up tmp files
    rm $out_dir/tmp_*
fi

# In "both" or "impute" mode, imputation of the unmodified input is conducted
if [[ $mode != "acc" ]]
then
    # Run BEAGLE
    java -jar $BEAGLE \
    gt=$input_vcf \
    out=$output \

    # Unzip the imputed VCF (BEAGLE auto bgzips)
    bgzip -d -f $output.vcf.gz

    # Replace phased pipes with unphased forward slashes
    sed -i 's/|/\//g' $output.vcf
fi