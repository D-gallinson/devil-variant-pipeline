#!/bin/bash
#SBATCH --job-name=5_align
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture5/out/5_align/%a_align.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture5/err/5_align/%a_align.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=26G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-192%14

# mem was previously: --mem=23800M

module purge
module add apps/bwa

source main.env
source refs.env
source tools.env

########################## Note #############################
# This script processes each sample in parallel,
# but iterates through the PE lanes of a sample.
# It also utilizes the _S1_, _S2_, _S3_, etc convention
# within each sample's name, where _S#_ refers to the sample 
# number. Thus, --array should begin at the first _S#_ and 
# end at the last _S#_. If this naming convention changes
# then the script must be modified to reflect this.
#############################################################

input=${DATA}/${batch}/3_trim
output_intermediate=${DATA}/${batch}/5_align
output_result=${RESULTS}/${batch}/align

forward_array=(${input}/*_S${SLURM_ARRAY_TASK_ID}_*R1*)
microchip_id=$(echo ${forward_array[0]} | awk '{print $NF}' FS=/ | grep -o -E '[TH]{0,1}[0-9A-Za-z]*-{0,1}[0-9]+' | head -n 1) #Get microchip ID (including tumorID if it exists)

start=`date +%s`

# Loop setup for NextSeq runs, allows this script to be run if a single sample has multiple lanes (e.g., multiple fastq files)
for i in {0..0}
do
    forward=${forward_array[$i]}
    reverse=${forward/R1/R2}
    reverse=${reverse/val_1/val_2}
    lane=$(echo $forward | grep -P -o 'L\d+')
    outname=${microchip_id}_${lane}

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
        $BWA_REF \
        $forward $reverse \
        > $sam_out

    #Generate alignment stats
    $SAMTOOLS flagstat \
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
    $(printf "INPUT=%s " ${output_intermediate}/${microchip_id}_L*.sorted.bam) \
    OUTPUT=${output_intermediate}/${microchip_id}.dups.bam \
    METRICS_FILE=${output_result}/duplicates/${microchip_id}.dups.stats.txt \
    VALIDATION_STRINGENCY=SILENT

rm ${output_intermediate}/${microchip_id}_L*.sorted.bam

#Generate BAM index (.csi due to large contig size)
$SAMTOOLS index \
    -@ 3 \
    -c \
    ${output_intermediate}/${microchip_id}.dups.bam


#Collect probe/target alignment metrics
java -jar $PICARD CollectHsMetrics \
    INPUT=${output_intermediate}/${microchip_id}.dups.bam \
    OUTPUT=${output_result}/HS-metrics/${microchip_id}.txt \
    BAIT_INTERVALS=$PROBE \
    TARGET_INTERVALS=$TARGET

end=`date +%s`
printf "Execution time: $((end-start))s"