# This will combine variants for each sample

import sys, os
import fileinput
import sys
import subprocess
import time
import argparse as ag
import re

HOME = "/home/d/dgallinson"
WORK_BGFS = "/work_bgfs/d/dgallinson"
SHARES = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data"
BASE = f"{SHARES}/joint-variants"

gatk = f"{HOME}/tools/gatk-4.2.0.0/./gatk"
input = os.getcwd()
samples = os.listdir(input)
index = f"{WORK_BGFS}/data/Sarcophilus_harrisii.mSarHar1.11.dna_sm.toplevel.fa"
intervals = f"{WORK_BGFS}/data/intervals/intervals_chr.list"

parser = ag.ArgumentParser(description="Rename files specified by a CSV")
parser.add_argument("mode", help="Use \"T\" for processing tumor samples and \"H\" for processing host samples")
args = parser.parse_args()

if args.mode.lower() == "h":
	pattern = "^(?!T).*"
	sample_type = "host"
elif args.mode.lower() == "t":
	pattern = "T"
	sample_type = "tumor"
else:
	print("Error, please specify a \"T\" for tumor samples or an \"H\" host samples. Exiting.")
	exit()

genomicsDB_path = f"{BASE}/genomicsDB-{sample_type}"
samples = [file for file in samples if file.endswith(".g.vcf") and re.match(pattern, file)]
print("Found " + str(len(samples)) + " samples")

#GVCF consolidation
start = time.perf_counter()
sampleList = " "
for sample in samples[8:]:
    sampleList = sampleList + "-V " + sample + " "
# for sample in samples:
#     sampleList = sampleList + "-V " + sample + " "

# Check below for RAM
# command = f"{gatk} --java-options \"-Xmx30g\" GenomicsDBImport {sampleList} --genomicsdb-workspace-path {genomicsDB_path} --intervals {intervals} --genomicsdb-shared-posixfs-optimizations true"
command = f"{gatk} --java-options \"-Xmx30g\" GenomicsDBImport {sampleList} --genomicsdb-update-workspace-path {genomicsDB_path} --genomicsdb-shared-posixfs-optimizations true"
print(command)
subprocess.call(command,shell=True)
subprocess.call(f"chmod -R g+rwx {genomicsDB_path}", shell=True)

delta = time.perf_counter() - start
print(f"GenomicsDBImport execution time: {delta:.0f}s")

#Joint genotyping
start = time.perf_counter()
geno_out = f"{BASE}/output-{sample_type}.vcf"

command = f"{gatk} GenotypeGVCFs -R {index} -V gendb://{genomicsDB_path} -O {geno_out}"
print(command)
subprocess.call(command,shell=True)

delta = time.perf_counter() - start
print(f"GenotypeGVCFs execution time: {delta:.0f}s")


# Filter SNPs
start = time.perf_counter()

raw_snp_out = f"{BASE}/raw_snps-{sample_type}.vcf"
command = f"{gatk} SelectVariants -R {index} -V {geno_out} --select-type-to-include SNP -O {raw_snp_out}"
subprocess.call(command,shell=True)

filtered_snps_out = f"{BASE}/filtered_snps-{sample_type}.vcf"
command = f"{gatk} VariantFiltration -R {index} -V {raw_snp_out} --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filter-name \"FILTER\" -O {filtered_snps_out}"
subprocess.call(command,shell=True)

# Filter indels
raw_indel_out = f"{BASE}/raw_indels-{sample_type}.vcf"
command = f"{gatk} SelectVariants -R {index} -V {geno_out} --select-type-to-include INDEL -O {raw_indel_out}"
subprocess.call(command,shell=True)

filtered_indels_out = f"{BASE}/filtered_indels-{sample_type}.vcf"
command = f"{gatk} VariantFiltration -R {index} -V {raw_indel_out} --filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filter-name \"FILTER\" -O {filtered_indels_out}"
subprocess.call(command,shell=True)

# Generate final variant files without the filtered variants
command = f"grep -E \'^#|PASS\' {filtered_snps_out} > {BASE}/FINAL_SNPs-{sample_type}.vcf"
subprocess.call(command,shell=True)

command = f"grep -E \'^#|PASS\' {filtered_indels_out} > {BASE}/FINAL_INDELs-{sample_type}.vcf"
subprocess.call(command,shell=True)

delta = time.perf_counter() - start
print(f"SelectVariants/VariantFiltration execution time: {delta:.0f}s")
