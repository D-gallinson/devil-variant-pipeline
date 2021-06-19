# devil-variant-pipeline
Pipeline to process NGS probe-capture data derived from Tasmanian Devils and DFTD samples.
The pipeline begins with raw FastQ files and generates joint variant calls for the Devil
hosts and the tumors.

===Primary Pipeline===
1_setup.sh [Move reads into a single folder, rename read files, and generate batch directory structure]
2_pre-qc.sh [Perform initial QC on the raw reads using FastQC and MultiQC]
3_trim.sh [Trim reads using TrimGalore! at settings determined based on 2_pre-qc results]
4_post-qc.sh [Collate FastQC files from TrimGalore! and run MultiQC on these]
5_align.sh [Align reads and generate sorted BAMs with marked duplicates as follows: BWA MEM > Picard SortSam > Picard MarkDuplicates]
6_combine.sh [Combine reports from 5_align using MultiQC and a custom script for summarizing Picard CollectHsMetrics]
7_gvcf.sh [Generate variant calls using HaplotypeCaller]
8_joint-variants.sh [Joint variant calling using GenomicsDBImport and GenotypeGVCFs]