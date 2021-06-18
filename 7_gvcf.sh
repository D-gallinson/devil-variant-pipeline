#!/bin/bash
#SBATCH --job-name=7_gvcf
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --output=scripts/master/logs/Capture1_6-11-21/out/7_gvcf/%a_gvcf.out
#SBATCH --error=scripts/master/logs/Capture1_6-11-21/err/7_gvcf/%a_gvcf.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=23808M
#SBATCH --time=6-23:00:00
#SBATCH --array=0-191

batch=Capture1_6-11-21   #change between runs
SHARES_BGFS=/shares_bgfs/margres_lab
index=${WORK_BGFS}/data/Sarcophilus_harrisii.mSarHar1.11.dna_sm.toplevel.fa

input_array=(${SHARES_BGFS}/${batch}/5_align/*dups.bam)
input=${input_array[$SLURM_ARRAY_TASK_ID]}

output_dir=${SHARES_BGFS}/Devils/BEE_Probe_Data/gvcf
output_name=$(echo $input | awk '{printf $NF}' FS=/ | grep -o '^[^\.]*') #Grab T#ID_microchip (or just microchip)
output=${output_dir}/${output_name}.g.vcf

start=`date +%s`

${HOME}/tools/gatk-4.2.0.0/./gatk \
    HaplotypeCaller \
    -R ${index} \
    -I ${input} \
    -O ${output} \
    -ERC GVCF \
    --do-not-run-physical-phasing

end=`date +%s`
printf "Execution time: $((end-start))s"