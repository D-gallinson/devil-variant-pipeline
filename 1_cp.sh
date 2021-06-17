#!/bin/bash
#SBATCH --job-name=1_cp
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=scripts/master/logs/1_cp.out
#SBATCH --error=scripts/master/logs/1_cp.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=186G
#SBATCH --time=02:00:00

module purge
module load apps/python/3.8.5

batch=Capture1_6-11-21   #change between runs
csv=Cap_libraries_pools_1-2.csv   #change between runs

share=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data
read_dir=$(ls ${share}/${batch} | grep -v ".csv$" | head -1)
read_dir=${read_dir/_R[01]/}
dest=data

mkdir ${dest}/${batch}
mkdir ${dest}/${batch}/${read_dir}

#Copy files
cp ${share}/${batch}/*/* ${dest}/${batch}/${read_dir}
cp ${share}/${batch}/*.csv ${dest}/${batch}

#Rename files
name_csv=${WORK_BGFS}/data/${batch}/${csv}
utility=${WORK_BGFS}/scripts/master/utility

python3 scripts/master/utility/rename.py \
    --id-col "Library number" \
    --rename-col "ID" \
    -r \
    ${dest}/${batch} \
    $name_csv


#Generate dir structure in outputs/results
results=${WORK_BGFS}/outputs/results
subdirs=(qc qc/pre qc/post qc/post/fastqc qc/post/trimming qc/post/MultiQC align align/flagstat align/duplicates align/HS-metrics align/MultiQC)

mkdir ${results}/${batch}
for subdir in ${subdirs[@]}
do
    mkdir ${results}/${batch}/${subdir}
done

#Generate dir structure in outputs/intermediates
intermediates=${WORK_BGFS}/outputs/intermediates
subdirs=(3_trim 5_align)

mkdir ${intermediates}/${batch}
for subdir in ${subdirs[@]}
do
    mkdir ${intermediates}/${batch}/${subdir}
done

#Generate dir structure in logs
logs=${WORK_BGFS}/scripts/master/logs
subdirs=(3_trim 5_align 7_gvcf)

mkdir ${logs}/${batch}
mkdir ${logs}/${batch}/out
mkdir ${logs}/${batch}/err
for subdir in ${subdirs[@]}
do
    mkdir ${logs}/${batch}/out/${subdir}
    mkdir ${logs}/${batch}/err/${subdir}
done