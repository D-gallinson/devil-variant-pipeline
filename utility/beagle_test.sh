module purge
module load apps/vcftools

# Configurable
input_vcf="../FINAL_SNPs-host.vcf"
missing=0.5

# Constants
TRUTH="tmp_truth"
MODEL="tmp_model"

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
rm start.vcf end.vcf ${TRUTH}.recode.vcf ${TRUTH}.GT.FORMAT

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