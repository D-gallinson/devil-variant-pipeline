#!/bin/bash
#SBATCH --job-name=3_trim
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --output=scripts/master/logs/Capture1_6-11-21/out/3_trim/%a_trim.out
#SBATCH --error=scripts/master/logs/Capture1_6-11-21/err/3_trim/%a_trim.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7900M
#SBATCH --time=6-23:00:00
#SBATCH --array=0-191

######################## NOTE ###############################
# $batch_name must be the name of the specific batch and
# $batch must point to the directory containing L#ID folders
#
# --array should start at 0 and go to num_reads/2 - 1.
# If there are 128 PE reads, --array=0-63. If there are 192
# samples, then (192*8)/2 - 1 = 767, thus --array=0-767
# Original run had --mem=4G, so 72 parallel jobs at 7.75G for
# a 192 batch (--array=0-767%72, --mem=7.75G)
#############################################################

module purge
module add apps/trimgalore/0.4.4

batch=Capture1_6-11-21  #Change between runs
SHARES_BGFS=/shares_bgfs/margres_lab/Devils/BEE_Probe_Data
input_dir=${SHARES_BGFS}/${batch}
output=${SHARES_BGFS}/${batch}/3_trim

#Generate an array of all R1 reads in $input_dir. The PE R2 read is generated through string replacement
forward_array=(${input_dir}/reads_1/*R1*.fastq.gz)
forward=${forward_array[$SLURM_ARRAY_TASK_ID]}
reverse=${forward//R1/R2}

start=`date +%s`

#Trimming and FastQC
trim_galore \
    -o $output \
    --fastqc \
    --paired \
    $forward $reverse

end=`date +%s`
printf "Execution time: $((end-start))s"