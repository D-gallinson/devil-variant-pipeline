chroms="${WORK_BGFS}/data/devil_chromosomes.txt"
vcf="${WORK_BGFS}/outputs/intermediates/8_joint-variants/impute/truth_clean_final.vcf"
out="truth_clean.chroms"

printf "" > $out
while read l
do
    echo "$l = `grep -c -P "^$l" $vcf`" >> $out
done <$chroms
