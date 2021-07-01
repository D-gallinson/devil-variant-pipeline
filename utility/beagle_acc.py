import numpy as np
import pandas as pd
import argparse as ag
import math
import subprocess
import time


pd.options.mode.chained_assignment = None

# =========================== FUNCTIONS ===========================
def test_even(s):
	evens = s % 2
	evens = evens.to_numpy()
	return (evens[0] == evens).all()


def unique_counts(series):
	if not isinstance(series, np.ndarray):
		series = series.to_numpy()
	unique, counts = np.unique(series, return_counts=True)
	count_dict = dict(zip(unique, counts))
	return count_dict


def get_freq(missing, samples, prec=3, display=False):
	count_dict = unique_counts(missing)
	total = sum(count_dict.values())
	if display:
		print(f"Total: {total}")
		for k, v in count_dict.items():
			print(f"{k}: {v} [{round(v/total, 3)}]")
	return {(samples - k): v/total for k, v in count_dict.items()}


def missing_vcf(genos, missing_list):
	missing_genos = genos.copy()
	for i in range(len(missing_list)):
		missing = missing_list[i]
		if missing.size == 0:
			continue
		missing_genos.iloc[i, missing] = "./."
	return missing_genos


# =========================== MAIN ===========================
def main(input, freq_model_file, output="testset.vcf", display_mode=False, show_stats=False, save_test_set=False, nrows=None):
	# Constants
	HOME = "/home/d/dgallinson"
	WORK_BGFS = "/work_bgfs/d/dgallinson"
	BEAGLE = f"{HOME}/tools/beagle.jar"

	if display_mode:
		start = time.perf_counter()
		print(f"Loading VCF file \"{input}\" as dataframe...", end="", flush=True)

	# Read in VCF to be filled with missing samples (this file should have no missing values to begin with)
	subprocess.call("sed -i '0,/#CHROM/{s/#CHROM/CHROM/}' " + input, shell=True)
	subprocess.call("sed -i 's/0\/1/1\/0/g' " + input, shell=True)
	df = pd.read_csv(input, sep="\t", comment="#", dtype=str, nrows=nrows)
	subprocess.call("sed -i '0,/CHROM/{s/CHROM/#CHROM/}' " + input, shell=True)
	df = df.rename(columns={"CHROM": "#CHROM"})
	samples = len(df.columns[9:])
	truth_geno_df = df[df.columns[9:]]
	total_genos = truth_geno_df.shape[0] * truth_geno_df.shape[1]

	if display_mode:
		delta = time.perf_counter() - start
		print(f"done in {delta:.2f}s")

	if display_mode:
		print(f"Obtaining missing frequencies (\"{freq_model_file}\" used as model)")
		print("---genotypes_missing: count [freq]---")

	# Use a file to obtain frequencies of missing samples (this should be the same file used to generate "input" 
	# but with --max-missing set to the desired level of allowed missing genotypes)
	missing_df = pd.read_csv(freq_model_file, nrows=nrows) // 2
	missing_freqs = get_freq(missing_df, samples, display=display_mode)

	if display_mode:
		print(f"Generating missing values...", end="", flush=True)

	# Artificially create missing samples (genotypes) in "input"
	missing_vector = np.random.choice(list(missing_freqs.keys()), p=list(missing_freqs.values()), size=len(truth_geno_df))
	total_missing = sum(missing_vector)
	missing_list = [np.random.choice(np.arange(samples), size=missing, replace=False) if missing != 0 else np.array([]) for missing in missing_vector]
	missing_genos = missing_vcf(truth_geno_df, missing_list)
	final_df = pd.concat([df[df.columns[:9]], missing_genos], axis=1)

	if display_mode:
		delta = time.perf_counter() - start
		print(f"done in {delta:.2f}s")
		print(f"Saved to {output}")

	if save_test_set:
		final_df.to_csv(output, index=False, sep="\t")

	# Impute missing genotypes and compare the imputed test VCF to the truth VCF
	unique_id = f"{time.time()}".replace(".", "")
	tmp_missing = f"missing_{unique_id}.vcf"
	tmp_impute = f"impute_{unique_id}"
	final_df.to_csv(tmp_missing, index=False, sep="\t")
	command = f"java -jar {BEAGLE} gt={tmp_missing} out={tmp_impute} > /dev/null"
	if display_mode:
		print(command)
	subprocess.call(command, shell=True)
	subprocess.call("gunzip " + f"{tmp_impute}.vcf.gz", shell=True)
	subprocess.call("sed -i 's/|/\//g' " + f"{tmp_impute}.vcf", shell=True)
	subprocess.call("sed -i '0,/#CHROM/{s/#CHROM/CHROM/}' " + f"{tmp_impute}.vcf", shell=True)
	subprocess.call("sed -i 's/0\/1/1\/0/g' " + f"{tmp_impute}.vcf", shell=True)
	impute_df = pd.read_csv(f"{tmp_impute}.vcf", dtype=str, sep="\t", comment="#", nrows=nrows)
	impute_df = impute_df.rename(columns={"CHROM": "#CHROM"})
	subprocess.call(f"rm {tmp_missing} {tmp_impute}.*", shell=True)

	impute_geno_df = impute_df[impute_df.columns[9:]]
	non_missing = total_genos - total_missing

	incorrect = (np.sum(truth_geno_df.compare(impute_geno_df).count())) // 2
	correct = total_missing - incorrect
	accuracy = (correct / total_missing) * 100
	if show_stats:
		print(f"Correctly imputed: {correct}")
		print(f"Total missing: {total_missing}")
		print(f"Accuracy: {accuracy:.2f}%")
	return accuracy


start = time.perf_counter()

parser = ag.ArgumentParser(description="Test imputation accuracy of BEAGLE")
parser.add_argument("truth", help="Path to the ground truth VCF (must have no missing samples)")
parser.add_argument("alleles", help="Path to the allele count file (generatef via: vcftools --counts | cut -f4)")
parser.add_argument("-o", "--output", dest="output", default="", help="Output name of the generated testing VCF, only specify if you wish to save this file [default=None]")
parser.add_argument("-d", "--display-mode", dest="display_flag", action="store_true", help="Flag to display processing steps")

parser.add_argument("-rc", "--rename-col", dest="rename_col_name", default="Microchip", help="Name of the rename column in the CSV file [default=Microchip]")
parser.add_argument("-tc", "--tissue-col", dest="tissue_col_name", default="Tissue", help="Name of the tissue column in the CSV file [default=Tissue]")
parser.add_argument("-tid", "--tid-col", dest="tumor_id_col_name", default="Tumour number", help="Name of the tumor ID column in the CSV file [default=Tumour number]")
parser.add_argument("-i", "--interactive", dest="interactive_flag", action="store_true", help="Use the tool in interactive mode")
parser.add_argument("-o", "--old-csv", dest="old_csv_flag", action="store_true", help="Use the tool on the old CSV format (the rename_col was spread between multiple columns)")
parser.add_argument("-r", "--recursive", dest="recur_flag", action="store_true", help="Flag to search subdirectories recursively")
parser.add_argument("-u", "--undo", dest="undo_flag", action="store_true", help="Undo renaming. Supply the parent directory containing all renamed files and the path to the rename_key.csv")
args = parser.parse_args()

iters = 10
accuracy = []
for i in range(iters):
	accuracy.append(main("truth.vcf", "SNPs.50.freq_model", "SNPs.50.test_set"))

delta = time.perf_counter() - start

mean = sum(accuracy) / len(accuracy)
std = math.sqrt(sum([(acc - mean)**2 for acc in accuracy]) / len(accuracy))
print(f"Accuracy mean: {mean:.2f}%")
print(f"Accuracy std: {std:.2f}%")
print(f"Execution time: {delta:.2f}s")