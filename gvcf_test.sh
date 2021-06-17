#!/bin/bash
#SBATCH --job-name=gvcf_uncompressed
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=scripts/master/logs/Capture1_6-11-21/out/7_gvcf/test_gvcf_uncompressed.out
#SBATCH --error=scripts/master/logs/Capture1_6-11-21/err/7_gvcf/test_gvcf_uncompressed.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=6-23:00:00
#SBATCH --nodelist=mdc-1057-4-6

batch=Capture1_6-11-21   #change between runs
index=${WORK_BGFS}/data/Sarcophilus_harrisii.mSarHar1.11.dna_sm.toplevel.fa

input_array=(${WORK_BGFS}/outputs/intermediates/${batch}/5_align/*dups.bam)
input=${input_array[0]}

output_dir=${WORK_BGFS}/outputs/intermediates/7_gvcf
output_name=$(echo $input | awk '{printf $NF}' FS=/ | grep -o '^[^\.]*') #Grab T#ID_microchip (or just microchip)
output=${output_dir}/TEST_${output_name}.g.vcf

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