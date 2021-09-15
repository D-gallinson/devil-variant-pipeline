# devil-variant-pipeline
Pipeline to process NGS probe-capture data derived from Tasmanian Devils and DFTD samples.<br/>
The pipeline begins with raw FastQ files and generates joint variant calls for the Devil<br/>
hosts and the tumors.

<h1>===Primary Pipeline===</h1><br/>
<b>1_setup.sh</b> [Move reads into a single folder, rename read files, and generate batch directory structure]<br/>
<b>2_pre-qc.sh</b> [Perform initial QC on the raw reads using FastQC and MultiQC]<br/>
<b>3_trim.sh</b> [Trim reads using TrimGalore! at settings determined based on 2_pre-qc results]<br/>
<b>4_post-qc.sh</b> [Collate FastQC files from TrimGalore! and run MultiQC on these]<br/>
<b>5_align.sh</b> [Align reads and generate sorted BAMs with marked duplicates as follows: BWA MEM > Picard SortSam > Picard MarkDuplicates]<br/>
<b>6_combine.sh</b> [Combine reports from 5_align using MultiQC and a custom script for summarizing Picard CollectHsMetrics]<br/>
<b>7_gvcf.sh</b> [Generate variant calls using HaplotypeCaller]<br/>
<b>8_genomicsDB.sh</b> [Create/update a tumor or host genomicsDB]<br/>
<b>9a_genotypeGVCFs.sh</b> [Run GenotypeGVCFs, processing equal bp regions of the genome in parallel]<br/>
<b>9b_genotypeGVCFs.sh</b> [Combine the jointly genotyped regions, separate indels and SNPs, and apply basic hard filters]<br/>
<b>10_variant_stats.sh</b> [Generate tables and figures for a VCF (this is essential for proper filtering)]<br/>
<b>ATOMM.sh</b> [TODO: automate ATOMM run, including subsetting SNPs for interaction test via marginal tests]<br/>
<b>gemma.sh</b> [Run GEMMA BSLMM either on all samples or using leave one out blocked xval. TODO: generate pheno predictions, calculate accuracy]<br/>