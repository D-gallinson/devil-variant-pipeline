##################################################################
# Used on vcf or vcf.gz files to count the number of SNPs per
# chrom (and scaffold). This must be run from the master
# folder. Due to gross computational inefficiencies, it is
# best not to run this on large VCFs.
# 
# USAGE:
# To console: source utility/check_chroms.sh path_to_input
# To file: source utility/check_chroms.sh path_to_input anything!
# NOTE: output mode is represented by any string as a second arg
##################################################################

source refs.env

# Remove the file extension if it exists (for output mode only)
if echo "$1" | grep -q "\."
then
    prefix=$(echo $1 | awk 'BEGIN{FS=OFS="."} NF--')
else
    prefix="$1"
fi

# Output path (for output mode only)
out="${prefix}.check_chroms"

# If in output mode, create an empty file to append (prevents appending to an existing file)
if [[ $2 ]]; then printf "" > $out; fi

# Quick and dirty grep -c of each chrom/scaffold
# This means a file is grepped ~115 times, obviously very inefficient
while read l
do
    # Count for a given chrom/scaffold
    chrom=$(echo $l | cut -f1 -d " ")
    counts="$chrom = `zgrep -c -P \"^$chrom\" $1`"
    
    # Check if printing to stdout or saving to file
    if [[ $2 ]]
    then
        echo $counts >> $out
    else
        echo $counts
    fi
done < $CHROMS