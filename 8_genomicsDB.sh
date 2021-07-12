#!/bin/bash
#SBATCH --job-name=8_genomicsDV
#SBATCH --partition=margres_2020
#SBATCH --qos=margres20
#SBATCH --mail-user=dgallinson@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/work_bgfs/d/dgallinson/scripts/master/logs/Capture1_6-11-21/out/8_genotypeGVCFs.out
#SBATCH --error=/work_bgfs/d/dgallinson/scripts/master/logs/Capture1_6-11-21/err/8_genotypeGVCFs.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=14-00:00:00

source main.env
source refs.env
source tools.env

sample_type="tumor" # Set to "host" or "tumor"
db_cmd="update"     # Set to "update" or "new"

genomicsDB=$DATA/joint-variants/genomicsDB-${sample_type}
samples=$DATA/gvcf/T*.g.vcf
db_workspace="--genomicsdb-update-workspace-path"

if [[ $sample_type == "host" ]]
then
    samples=$DATA/gvcf/[0-9]*.g.vcf
fi

if [[ db_cmd == "new" ]]
then
    db_workspace="--genomicsdb-workspace-path"
fi

$GATK --java-options "-Xmx 30G" GenomicsDBImport \
    $(printf "-V %s " $samples) \
    $db_workspace $genomicsDB \
    --genomicsdb-shared-posixfs-optimizations true