#!/bin/bash

module load apps/python/3.8.5
source main.env

input=$RESULTS/filter/FINAL_SNPs-tumor.isec_raw.alleles.missing_100.vcf
tmp_dir=../../tmp
rename_key=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/Capture1_6-11-21/rename_key.csv
log=${input/.vcf/.HEADER.log}

line=$(grep -n -m1 '#CHROM' $input)
line_num=$(echo $line | cut -f1 -d:)
samples=$(echo $line | cut -f2 -d: | cut -f10- -d' ')

echo "===== ORIGINAL VCF HEADER =====" > $log
echo $line | cut -f2 -d: | sed 's/\s/\t/g' >> $log
echo -e "\n" >> $log

echo $samples | sed 's/\s/\n/g' > $tmp_dir/tmp_vcf_header.txt
python $SCRIPTS/utility/vcf_col_fix.py $tmp_dir/tmp_vcf_header.txt $rename_key $tmp_dir
replacement=$(echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" $(cut -f 1 $tmp_dir/tmp_vcf_fixed.txt) | sed -r "s/\s/\t/g")

echo "===== RENAMED VCF HEADER =====" >> $log
echo $replacement | sed 's/\s/\t/g' >> $log
echo -e "\n" >> $log
echo -e "original\trename" >> $log
paste $tmp_dir/tmp_vcf_header.txt $tmp_dir/tmp_vcf_fixed.txt >> $log

sed -i "${line_num}s/.*/$replacement/" $input

rm $tmp_dir/tmp_*