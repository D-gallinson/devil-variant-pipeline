# devil-variant-pipeline
Pipeline to process NGS probe-capture data derived from Tasmanian Devils and DFTD samples.<br/>
The pipeline begins with raw FastQ files and generates joint variant calls for the Devil<br/>
hosts and the tumors.

<h1>===Primary Pipeline===</h1><br/>
1_setup.sh [Move reads into a single folder, rename read files, and generate batch directory structure]<br/>
2_pre-qc.sh [Perform initial QC on the raw reads using FastQC and MultiQC]<br/>
3_trim.sh [Trim reads using TrimGalore! at settings determined based on 2_pre-qc results]<br/>
4_post-qc.sh [Collate FastQC files from TrimGalore! and run MultiQC on these]<br/>
5_align.sh [Align reads and generate sorted BAMs with marked duplicates as follows: BWA MEM > Picard SortSam > Picard MarkDuplicates]<br/>
6_combine.sh [Combine reports from 5_align using MultiQC and a custom script for summarizing Picard CollectHsMetrics]<br/>
7_gvcf.sh [Generate variant calls using HaplotypeCaller]<br/>
8_genomicsDB.sh [Create/update a tumor or host genomicsDB]<br/>
9a_genotypeGVCFs.sh [Run GenotypeGVCFs, processing equal bp regions of the genome in parallel]<br/>
9b_genotypeGVCFs.sh [Combine the jointly genotyped regions, separate indels and SNPs, and apply basic hard filters]<br/>
10_variant_stats.sh [Generate tables and figures for a VCF (this is essential for proper filtering)]<br/>
ATOMM.sh [TODO: automate ATOMM run, including subsetting SNPs for interaction test via marginal tests]<br/>
gemma.sh [Run GEMMA BSLMM either on all samples or using leave one out blocked xval. TODO: generate pheno predictions, calculate accuracy]<br/>