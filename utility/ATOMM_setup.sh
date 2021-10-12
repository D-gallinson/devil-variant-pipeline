#!/bin/bash
#SBATCH --job-name=ATOMM_prep
#SBATCH --partition=muma_2021
#SBATCH --qos=preempt_short
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/ATOMM/prelim/devil_survival/prep.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/ATOMM/prelim/devil_survival/prep.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module purge
module load apps/vcftools
module load apps/python/3.8.5

source main.env
source refs.env
source tools.env

host=$RESULTS/filter/FINAL_SNPs-host.minDP_10.maxDP_100.alleles.missing_100.mac_2.vcf.gz
pathogen=$RESULTS/filter/FINAL_SNPs-tumor.isec_raw.alleles.missing_100.vcf.gz
phenotype=$RESULTS/ATOMM/input/atomm_pheno.txt
subset_mode="true"

dir=$RESULTS/ATOMM/input
dir=../../tmp
host_out=$dir/tmp_host
pathogen_out=$dir/tmp_pathogen
log=$dir/atomm_prep.log

# Generate VCF --keep files for filtering NA samples
cut -f1 $phenotype | tail -n +2 > $dir/tmp_host_list.txt
cut -f2 $phenotype | tail -n +2 > $dir/tmp_pathogen_list.txt

# Filter out missing samples
vcftools \
    --gzvcf $host \
    --keep ${host_out}_list.txt \
    --recode \
    --out ${host_out}_filt

vcftools \
    --gzvcf $pathogen \
    --keep ${pathogen_out}_list.txt \
    --recode \
    --out ${pathogen_out}_filt

# Sort the filtered VCFs based on the phenotype file order
zgrep -m 1 '^##fileformat' $host > ${host_out}_template.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" $(cut -f 1 ${host_out}_list.txt) | sed -r "s/\s+/\t/g" >> ${host_out}_template.vcf
vcf-shuffle-cols -t ${host_out}_template.vcf ${host_out}_filt.recode.vcf > ${host_out}_final.vcf

zgrep -m 1 '^##fileformat' $pathogen > ${pathogen_out}_template.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" $(cut -f 1 ${pathogen_out}_list.txt) | sed -r "s/\s+/\t/g" >> ${pathogen_out}_template.vcf
vcf-shuffle-cols -t ${pathogen_out}_template.vcf ${pathogen_out}_filt.recode.vcf > ${pathogen_out}_final.vcf

if [[ $subset_mode != "true" ]]
then
    # Generate the 012 genotype matrix
    vcftools \
        --vcf ${host_out}_final.vcf \
        --012 \
        --out $host_out

    vcftools \
        --vcf ${pathogen_out}_final.vcf \
        --012 \
        --out $pathogen_out

    # Cut out the sample IDs (first column) and tranpose the matrix so that samples are columns and rows are variant sites
    cut -f2- $host_out.012 | datamash transpose > $host_out.T.012
    cut -f2- $pathogen_out.012 | datamash transpose > $pathogen_out.T.012

    # Replace 2's (homozygous alt) with 1's
    sed 's/2/1/g' $host_out.T.012 > $host_out.T.01
    sed 's/2/1/g' $pathogen_out.T.012 > $pathogen_out.T.01

    # Prepend the SNP ID column (consecutive increasing numbers from 1..N)
    wc -l $host_out.T.01 | 
    cut -d ' ' -f1 | 
    xargs seq | 
    paste - $host_out.T.01 > $host_out.ID.T.01

    wc -l $pathogen_out.T.01 | 
    cut -d ' ' -f1 | 
    xargs seq | 
    paste - $pathogen_out.T.01 > $pathogen_out.ID.T.01

    # Add the chr ID from the .pos file. X, Y, and unplaced scaffolds are converted as follows:
    # X = 7
    # Y = 8
    # unplaced = 9
    sed 's/[A-Z]*[0-9]*\.[0-9]/9/g; s/X/7/g; s/Y/8/g' $host_out.012.pos | 
    cut -f1 | 
    paste - $host_out.ID.T.01 > $dir/sequence_host.txt

    sed 's/[A-Z]*[0-9]*\.[0-9]/9/g; s/X/7/g; s/Y/8/g' $pathogen_out.012.pos | 
    cut -f1 | 
    paste - $pathogen_out.ID.T.01 > $dir/sequence_pathogen.txt

    # Replace tab-delimited to space-delimited
    sed -i -e 's/\t/ /g' ${dir}/sequence_host.txt
    sed -i -e 's/\t/ /g' ${dir}/sequence_pathogen.txt

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
            line=$(grep -m 1 -P '^(?!#)' ${1}_final.vcf)
        else
            line=$(tail -1 ${1}_final.vcf)
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
            line=$(head -1 ${dir}/sequence_${1}.txt)
        else
            line=$(tail -n 1 ${dir}/sequence_${1}.txt)
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
    grep -m 1 -P '^(?!##)' ${host_out}_final.vcf |
    cut -f10-29 |
    xargs printf "%-${width}s" >> $log
    echo "" >> $log

    # Display row 1, genotypes 1-20 from the input VCF file
    format_vcf $host_out "head"

    # Display row 1, genotypes 1-20 from the final genotype file
    format_seq "host" "head"

    # Display the final row, genotypes 1-20 from the input VCF file
    format_vcf $host_out "tail"

    # Display the final row, genotypes 1-20 from the final genotype file
    format_seq "host" "tail"

    echo -e "\n" >> $log

    # PATHOGEN
    echo "PATHOGEN" >> $log

    # Display Sample Names
    grep -m 1 -P '^(?!##)' ${pathogen_out}_final.vcf |
    cut -f10-29 |
    xargs printf "%-${width}s" >> $log
    echo "" >> $log

    # Display row 1, genotypes 1-20 from the input VCF file
    format_vcf $pathogen_out "head"

    # Display row 1, genotypes 1-20 from the final genotype file
    format_seq "pathogen" "head"

    # Display the final row, genotypes 1-20 from the input VCF file
    format_vcf $pathogen_out "tail"

    # Display the final row, genotypes 1-20 from the final genotype file
    format_seq "pathogen" "tail"


# ===================== GENERATE PHENOTYPE FILE =====================
python - << EOF
import pandas as pd

pheno_file = "$phenotype"
outdir = "$dir"

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
else
    # Change the final subsetted host/pathogen file names to avoid deletion if in subset_mode
    mv ${host_out}_final.vcf $dir/host_subset_final.vcf
    mv ${pathogen_out}_final.vcf $dir/pathogen_subset_final.vcf
fi

# Clean up the intermediate geno_ files
rm ${dir}/tmp_*