#!/bin/bash
#SBATCH --job-name=5_align
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --output=scripts/master/logs/Capture1_6-11-21/out/5_align/%a_align.out
#SBATCH --error=scripts/master/logs/Capture1_6-11-21/err/5_align/%a_align.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=23800M
#SBATCH --array=1-192
#SBATCH --time=04:00:00

module purge
module add apps/bwa

########################## Note #############################
# This script processes each sample in parallel,
# but iterates through the PE lanes of a sample.
# It also utilizes the _S1_, _S2_, _S3_, etc convention
# within each sample's name, where _S#_ refers to the sample 
# number. Thus, --array should begin at the first _S#_ and 
# end at the last _S#_. If this naming convention changes
# then the script must be modified to reflect this.
#############################################################

PICARD=${HOME}/tools/picard/build/libs/picard.jar
target=${WORK_BGFS}/data/intervals/All_targets_combined_WashU_Tasmanian_Devil_TE-91244716_SNP_Indel_Exons_backbone_SarHar1_1_new6_Picard.interval_list
probe=${WORK_BGFS}/data/intervals/Probe_Placement_WashU_Tasmanian_Devil_TE-91244716_SNP_Indel_Exons_backbone_SarHar1_1_197175_sorted_new6_Picard.intervals
ref=${WORK_BGFS}/outputs/intermediates/bwa-ref/S_harrisii

batch=Capture1_6-11-21 #change between runs
input=${WORK_BGFS}/outputs/intermediates/${batch}/3_trim
output_intermediate=${WORK_BGFS}/outputs/intermediates/${batch}/5_align
output_result=${WORK_BGFS}/outputs/results/${batch}/align

forward_array=(${input}/*_S${SLURM_ARRAY_TASK_ID}_*R1*)
microchip_id=$(echo ${forward_array[0]} | awk '{print $NF}' FS=/ | grep -o -E '[TH]{0,1}[0-9]*-{0,1}[0-9]+' | head -n 1) #Get microchip ID (including tumorID if it exists)

start=`date +%s`

# Loop setup for NextSeq runs, allows this script to be run if a single sample has multiple lanes (e.g., multiple fastq files)
for i in {0..0}
do
    let lane=$i+1
    outname=${microchip_id}_L00${lane}
    forward=${forward_array[$i]}
    reverse=${forward/R1/R2}
    reverse=${reverse/val_1/val_2}

    #Read group headers to be added in BWA
    header=$(zcat $forward | head -n 1)
    rg=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')

    sam_out=${output_intermediate}/${outname}.sam
    sortsam_out=${output_intermediate}/${outname}.sorted.bam

    #Generate alignments
    bwa mem \
        -t 3 \
        -M \
        -v 3 \
        -R "@RG\tID:${rg}\tSM:${microchip_id}\tLB:${microchip_id}_1\tPL:ILLUMINA" \
        $ref \
        $forward $reverse \
        > $sam_out

    #Generate alignment stats
    ${HOME}/tools/samtools-1.12/./samtools flagstat \
        -@ 3 \
        $sam_out \
        > ${output_result}/flagstat/${outname}-stats.txt

    #Generate sorted BAMs
    java -jar $PICARD SortSam \
        INPUT=${sam_out} \
        OUTPUT=${sortsam_out} \
        SORT_ORDER=coordinate

    rm $sam_out
done

#Mark duplicates and merge BAMs
java -jar $PICARD MarkDuplicates \
    $(printf "INPUT=%s " ${output_intermediate}/${microchip_id}_L00*.sorted.bam) \
    OUTPUT=${output_intermediate}/${microchip_id}.dups.bam \
    METRICS_FILE=${output_result}/duplicates/${microchip_id}.dups.stats.txt \
    VALIDATION_STRINGENCY=SILENT

rm ${output_intermediate}/${microchip_id}_L00*.sorted.bam

#Generate BAM index (.csi due to large contig size)
${HOME}/tools/samtools-1.12/./samtools \
    index \
    -@ 3 \
    -c \
    ${output_intermediate}/${microchip_id}.dups.bam

#Collect probe/target alignment metrics
java -jar $PICARD CollectHsMetrics \
    INPUT=${output_intermediate}/${microchip_id}.dups.bam \
    OUTPUT=${output_result}/HS-metrics/${microchip_id}.txt \
    BAIT_INTERVALS=$probe \
    TARGET_INTERVALS=$target

end=`date +%s`
printf "Execution time: $((end-start))s"