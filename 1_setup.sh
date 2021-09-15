#!/bin/bash
#SBATCH --job-name=1_setup
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/1_setup.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/1_setup.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=30:00

module purge
module load apps/python/3.8.5

set -euo pipefail

source main.env
source refs.env

data_path=${DATA}/${batch}

# Generate file containing the original read directory structure
# read_dirs=($(ls $data_path | grep -v '.[csvtxt]$'))
# out=${data_path}/reads_structure.txt
# echo "==========Directories==========" > $out
# ls ${data_path} | grep -v '.[csvtxt]$' >> $out
# printf "\n==========Files==========\n" >> $out
# for dir in ${read_dirs[@]}
# do
#     printf "$dir/\n" >> $out
#     ls ${data_path}/${dir} >> $out
#     echo >> $out
# done

# # Move R1 and R2 reads into a single directory
# mkdir ${data_path}/1_reads
# mv ${data_path}/*/* ${data_path}/1_reads

#Remove empty R1/R2 dirs
# rm $(find ${data_path}/* -type d -empty)

#Rename files
python3 ${SCRIPTS}/utility/rename.py \
    --id-col "Library number" \
    --rename-col "Microchip" \
    -r \
    $data_path \
    $SAMPLE_KEY

#Generate dir structure in outputs/results
subdirs=(qc qc/pre qc/post qc/post/fastqc qc/post/trimming qc/post/MultiQC align align/flagstat align/duplicates align/HS-metrics align/MultiQC)

mkdir ${RESULTS}/${batch}
for subdir in ${subdirs[@]}
do
    mkdir ${RESULTS}/${batch}/${subdir}
done

#Generate dir structure in SHARE_BGFS/batch
subdirs=(3_trim 5_align)
for subdir in ${subdirs[@]}
do
    mkdir ${data_path}/${subdir}
done

#Generate dir structure in logs
subdirs=(3_trim 5_align 7_gvcf 9a_GenotypeGVCFs)

mkdir ${LOGS}/${batch}
mkdir ${LOGS}/${batch}/out
mkdir ${LOGS}/${batch}/err
for subdir in ${subdirs[@]}
do
    mkdir ${LOGS}/${batch}/out/${subdir}
    mkdir ${LOGS}/${batch}/err/${subdir}
done