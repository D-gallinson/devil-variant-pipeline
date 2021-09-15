#!/bin/bash
#SBATCH --job-name=8_genomicsDB_tumor
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/out/8_genomicsDB_tumor.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture3/err/8_genomicsDB_tumor.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=186G
#SBATCH --time=21-00:00:00

source main.env
source refs.env
source tools.env

sample_type="tumor"  #Set to "host" or "tumor"
db_cmd="update"     #Set to "update" or "new"

genomicsDB=$DATA/joint-variants/genomicsDB-${sample_type}
samples=$DATA/gvcf/T*.g.vcf
db_workspace="--genomicsdb-update-workspace-path"

if [[ $sample_type == "host" ]]
then
    samples=$DATA/gvcf/H*.g.vcf
    samples=$WORK_BGFS/outputs/intermediates/7_gvcf/H*.g.vcf
fi

if [[ $db_cmd == "new" ]]
then
    db_workspace="--genomicsdb-workspace-path"
    interval_cmd="--intervals"
    interval_file=$CHROMS
fi

$GATK --java-options "-Xmx50G" GenomicsDBImport \
    $(for sample in $samples; do echo -n "-V $sample "; done) \
    $db_workspace $genomicsDB \
    --genomicsdb-shared-posixfs-optimizations true \
    $interval_cmd $interval_file

if [[ $db_cmd == "new" ]]
then
    chmod -R g+rwx $genomicsDB
fi