#!/bin/bash
#SBATCH --job-name=7_gvcf
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/out/7_gvcf/%a_gvcf.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/err/7_gvcf/%a_gvcf.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=23808M
#SBATCH --time=7-00:00:00
#SBATCH --array=0-191

source main.env
source refs.env
source tools.env

input_array=(${DATA}/${batch}/5_align/*dups.bam)
input=${input_array[$SLURM_ARRAY_TASK_ID]}

output_dir=${DATA}/gvcf
output_name=$(echo $input | awk '{printf $NF}' FS=/ | grep -o '^[^\.]*') #Grab T#ID_microchip (or just microchip)
output=${output_dir}/${output_name}.g.vcf

start=`date +%s`

${GATK} HaplotypeCaller \
    -R ${REF} \
    -I ${input} \
    -O ${output} \
    -ERC GVCF \
    --do-not-run-physical-phasing

end=`date +%s`
printf "Execution time: $((end-start))s"