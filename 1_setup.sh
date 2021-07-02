#!/bin/bash
#SBATCH --job-name=1_setup
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=scripts/master/logs/Capture1_6-11-21/out/1_setup.out
#SBATCH --error=scripts/master/logs/Capture1_6-11-21/err/1_setup.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=30:00

module purge
module load apps/python/3.8.5

source ${WORK_BGFS}/scripts/master/main.env

data_path=${DATA}/${batch}

# Generate file containing the original read directory structure
read_dirs=($(ls $data_path | grep -v '.[csvtxt]$'))
out=${data_path}/reads_structure.txt
echo "==========Directories==========" > $out
ls ${data_path} | grep -v '.[csvtxt]$' >> $out
printf "\n==========Files==========\n" >> $out
for dir in ${read_dirs[@]}
do
    printf "$dir/\n" >> $out
    ls ${data_path}/${dir} >> $out
    echo >> $out
done

# Move R1 and R2 reads into a single directory
mkdir ${data_path}/1_reads
mv ${data_path}/*/* ${data_path}/1_reads

#Remove empty R1/R2 dirs
rm -d ${data_path}/*/

#Rename files
name_csv=${data_path}/${csv}

python3 ${SCRIPTS}/utility/rename.py \
    --id-col "Library number" \
    --rename-col "ID" \
    -r \
    $data_path \
    $name_csv

#Generate dir structure in outputs/results
subdirs=(qc qc/pre qc/post qc/post/fastqc qc/post/trimming qc/post/MultiQC align align/flagstat align/duplicates align/HS-metrics align/MultiQC)

mkdir ${results}/${batch}
for subdir in ${subdirs[@]}
do
    mkdir ${results}/${batch}/${subdir}
done

#Generate dir structure in SHARE_BGFS/batch
subdirs=(3_trim 5_align)
for subdir in ${subdirs[@]}
do
    mkdir ${data_path}/${subdir}
done

#Generate dir structure in logs
subdirs=(3_trim 5_align 7_gvcf)

mkdir ${LOGS}/${batch}
mkdir ${LOGS}/${batch}/out
mkdir ${LOGS}/${batch}/err
for subdir in ${subdirs[@]}
do
    mkdir ${LOGS}/${batch}/out/${subdir}
    mkdir ${LOGS}/${batch}/err/${subdir}
done