import pandas as pd
import vcf
import time
import subprocess
from io import StringIO
from collections import OrderedDict



def get_lines(targets, vcf_path, host=True):
	start = time.perf_counter()
	tissue = "host" if host else "pathogen"
	print(f"1) Searching for lines in {tissue} VCF...", end="", flush=True)
	cat = "cat"
	if ".gz" in vcf_path:
		cat = "zcat"
	vcf_version = subprocess.check_output(f"{cat} {vcf_path} | head -1", shell=True).decode('utf-8').strip("\n")
	colnames = subprocess.check_output(f"{cat} {vcf_path} | grep -v -m1 '^##'", shell=True).decode('utf-8').strip("\n")
	sample_names = colnames.split("\t")[9:]
	header = [vcf_version, colnames]
	header_lines = subprocess.check_output(f"{cat} {vcf_path} | head -1000 | grep -c '^#'", shell=True)
	header_lines = int(header_lines.decode('utf-8').strip("\n"))
	target_str = [f"{target + header_lines}p;" for target in targets]
	target_str = "".join(target_str)
	target_str = target_str[:len(target_str)-1]
	found_lines = subprocess.check_output(f"{cat} {vcf_path} | sed -ne '{target_str}'", shell=True)
	found_lines = found_lines.decode('utf-8').split("\n")
	found_lines = found_lines[:len(found_lines) - 1]
	found_lines = [line for _, line in sorted(zip(targets, found_lines))]
	found_lines = header + found_lines
	found_str = "\n".join(found_lines)
	found_file_obj = StringIO(found_str)
	delta = time.perf_counter() - start
	print(f"done ({delta:.3f}s)")
	print(f"   Found {len(found_lines)-2} lines")
	return found_file_obj, sample_names


def vcf_df(file_like, sample_names):
	start = time.perf_counter()
	print(f"2) Converting lines to DF...", flush=True, end="")
	vcf_reader = vcf.Reader(file_like)
	vcf_dict = OrderedDict([(sample, []) for sample in sample_names])
	vcf_dict["CHROM"] = []
	vcf_dict["POS"] = []
	vcf_dict["REF"] = []
	vcf_dict["ALT"] = []
	vcf_dict.move_to_end("ALT", last=False)
	vcf_dict.move_to_end("REF", last=False)
	vcf_dict.move_to_end("POS", last=False)
	vcf_dict.move_to_end("CHROM", last=False)

	for record in vcf_reader:
		vcf_dict["CHROM"].append(record.CHROM)
		vcf_dict["POS"].append(record.POS)
		vcf_dict["REF"].append(record.REF)
		vcf_dict["ALT"].append(record.ALT[0])
		haploid = [int(sample["GT"][0]) + int(sample["GT"][-1]) for sample in record.samples]
		haploid = [geno if geno < 2 else 1 for geno in haploid]
		samples = [sample.sample for sample in record.samples]
		for i in range(len(haploid)):
			sample = samples[i]
			haploid_geno = haploid[i]
			vcf_dict[sample].append(haploid_geno)

	delta = time.perf_counter() - start
	df = pd.DataFrame(vcf_dict)
	print(f"done ({delta:.3f}s)")
	return df


def geno_check(seq_df, vcf_df):
	start = time.perf_counter()
	print("4) Checking if VCF and ATOMM sequence haploid genotypes match...", flush=True, end="")
	seq_genos = seq_df.iloc[:, 2:(seq_df.shape[1]-1)].reset_index(drop=True)
	vcf_genos = vcf_df.iloc[:, 4:(vcf_df.shape[1]-1)]
	vcf_genos.columns = seq_genos.columns
	full_df_check = vcf_genos.equals(seq_genos)
	check_list = [True for i in range(seq_genos.shape[0])]
	if not full_df_check:
		check_list = []
		for i in range(seq_genos.shape[0]):
			check_list.append((seq_genos.iloc[1] == vcf_genos.iloc[1]).all())
	delta = time.perf_counter() - start
	print(f"done ({delta:.3f}s)")
	if not any(check_list):
		print("   All lines failed to match")
		return []
	else:
		sums = seq_genos.sum(axis=1)
		all_ref = sums[sums == 0]
		all_alt = sums[sums == seq_genos.shape[1]]
		successful_lines = [i+1 for i in range(len(check_list)) if check_list[i]]
		not_list = [not val for val in check_list]
		failed_lines = [i+1 for i in range(len(not_list)) if not_list[i]]
		if not all_ref.empty:
			all_ref_i = all_ref.index
			all_ref_i = [i for i in all_ref_i if i in successful_lines]
			if all_ref_i:
				print(f"   NOTE: lines {','.join(all_ref_i)} are all 0s (match(s) could be spurious)")
		if not all_alt.empty:
			all_alt_i = all_alt.index
			all_alt_i = [i for i in all_alt_i if i in successful_lines]
			if all_alt_i:
				print(f"   NOTE: lines {','.join(all_alt_i)} are all 1s (match(s) could be spurious)")
		if not all(check_list):
			print(f"   Failed matches in the following lines: {','.join(failed_lines)} ({len(failed_lines)} failures)")
		else:
			print("   All lines matched successfully")
		return [line-1 for line in successful_lines]

# TODO
# 1) CLI-friendly
# 2) Hap_maf ID key compatible. Also fix hap_maf.R so it produces the key
# 3) Ensure geno_check works, no idea if all those damned if/elses do what I want
# 4) Allow single-file (i.e., marginal tests) inputs
# 5) A lot of the code feels redundant
# 6) Any error checking
# 7) Comments
# 8) Also this currentlu just prints the results, do something better
# 9) Better reporting of failed matches
# 10) Perhaps line counting to see if a failure will immediately happen or not I don't know how I feel about this

# SETTINGS
interaction_file = True
target_path = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/results/ATOMM/output/top_5_interaction.txt"
host_vcf = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/results/ATOMM/infection_age/attempt3/input_unfiltered/VCFs/host_subset_final.missing_100.vcf.gz"
pathogen_vcf = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/results/ATOMM/infection_age/attempt3/input_unfiltered/VCFs/pathogen_subset_final.missing_100.alleles.vcf.gz"
host_ATOMM = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/results/ATOMM/infection_age/attempt3/input_unfiltered/sequence_host_ORIGINAL.txt"
pathogen_ATOMM = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data/results/ATOMM/infection_age/attempt3/input_unfiltered/sequence_pathogen_ORIGINAL.txt"


target_df = pd.read_csv(target_path, sep="\t")
if interaction_file:
	host_targets = target_df.iloc[:, 1].unique()
	pathogen_targets = target_df.iloc[:, 3].unique()

	print(f"Found {len(host_targets)} unique host and {len(pathogen_targets)} unique pathogen targets")
	print("===== HOST =====")
	host_hits, host_names = get_lines(host_targets, host_vcf)
	host_df = vcf_df(host_hits, host_names)
	seq_start = time.perf_counter()
	print("3) Reading sequence_host and extracting target lines...", flush=True, end="")
	seq_host_df = pd.read_csv(host_ATOMM, header=None, sep=" ")
	seq_host_df = seq_host_df.loc[host_targets-1]
	seq_delta = time.perf_counter() - seq_start
	print(f"done ({seq_delta:.3f}s)")
	matching_lines = geno_check(seq_host_df, host_df)
	host_final = host_df.iloc[matching_lines][["CHROM", "POS", "REF", "ALT"]]
	host_final["ATOMM_ID"] = host_targets[matching_lines]
	print(host_final.to_string(index=False))

	print("\n===== PATHOGEN =====")
	pathogen_hits, pathogen_names = get_lines(pathogen_targets, pathogen_vcf, host=False)
	pathogen_df = vcf_df(pathogen_hits, pathogen_names)
	seq_start = time.perf_counter()
	print("3) Reading sequence_pathogen and extracting target lines...", flush=True, end="")
	seq_pathogen_df = pd.read_csv(pathogen_ATOMM, header=None, sep=" ")
	seq_pathogen_df = seq_pathogen_df.loc[pathogen_targets-1]
	seq_delta = time.perf_counter() - seq_start
	print(f"done ({seq_delta:.3f}s)")
	matching_lines = geno_check(seq_pathogen_df, pathogen_df)
	pathogen_final = pathogen_df.iloc[matching_lines][["CHROM", "POS", "REF", "ALT"]]
	pathogen_final["ATOMM_ID"] = pathogen_targets[matching_lines]
	print(pathogen_final.to_string(index=False))
else:
	targets = target_df.iloc[:, 1]