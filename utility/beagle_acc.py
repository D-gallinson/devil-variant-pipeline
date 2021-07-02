# Test BEAGLE imputation accuracy using a custom VCF input file
import numpy as np
import pandas as pd
import argparse as ag
import math
import subprocess
import time
from utility import Env


pd.options.mode.chained_assignment = None

# =========================== CLASSES/FUNCTIONS ===========================
# Class to generate the probability model
class MissingModel:
	def __init__(self):
		pass


	def load_model(self, allele_file, samples, nrows=None):
		allele_df = pd.read_csv(allele_file, nrows=nrows)
		self.__check_model(allele_df)
		allele_df = allele_df // 2
		return self.__get_freq(allele_df, samples)


	def __check_model(self, model):
		if model.shape[1] != 1:
			raise ModelError(f"allele file should have only a single column (not {model.shape[1]} cols.")
		try:
			col1 = model.columns[0]
			even_test = model % 2
			odds = even_test[even_test[col1] != 0]
			if odds.size > 0:
				print(f"Found odd allele counts in the following {odds.size} rows:")
				print(", ".join(odds.index.astype(str)))
				raise ModelError("detected odd allele counts (all should be even). Check these rows and consider removing them.", False)
		except TypeError:
			raise ModelError("some or all values are not numbers.")


	def __get_freq(self, missing, samples):
		count_dict = self.__unique_counts(missing)
		total = sum(count_dict.values())
		print("==========Counts/frequencies of samples for all sites==========")
		print(f"#Total samples: {samples}")
		print(f"#Total sites: {total}")
		print("SAMPLES\tCOUNT\tFREQ")
		for k, v in count_dict.items():
			print(f"{k}\t{v}\t{round(v/total, 3)}")
		print("===============================================================")
		return {(samples - k): v/total for k, v in count_dict.items()}


	def __unique_counts(self, series):
		if not isinstance(series, np.ndarray):
			series = series.to_numpy()
		unique, counts = np.unique(series, return_counts=True)
		count_dict = dict(zip(unique, counts))
		return count_dict



class ModelError(Exception):
	def __init__(self, error, default=True):
		self.error = error
		self.default = default

	def __str__(self):
		str = f"Problem generating frequency model: {self.error}"
		if self.default:
			str += "Try \"vcftools --counts | cut -f4\""
		return str



# Load the ground truth VCF
def truth_vcf(path, nrows=None):
	subprocess.call("sed -i '0,/#CHROM/{s/#CHROM/CHROM/}' " + path, shell=True)
	df = pd.read_csv(path, sep="\t", comment="#", dtype=str, nrows=nrows)
	subprocess.call("sed -i '0,/CHROM/{s/CHROM/#CHROM/}' " + path, shell=True)
	df = df.rename(columns={"CHROM": "#CHROM"})
	return df


# Convert genotypes to missing to generate the test VCF
def missing_vcf(genos, missing_list):
	missing_genos = genos.copy()
	for i in range(len(missing_list)):
		missing = missing_list[i]
		if missing.size == 0:
			continue
		missing_genos.iloc[i, missing] = "./."
	return missing_genos


# =========================== MAIN ===========================
def main(truth_df, freq_model, tmp_id="", output="", phased=False, nrows=None):
	# Constants
	tools_env = Env("tools.env")
	BEAGLE = tools_env.get_var("BEAGLE")

	truth_geno_df = truth_df[truth_df.columns[9:]]
	samples = truth_geno_df.shape[1]
	total_genos = truth_geno_df.shape[0] * truth_geno_df.shape[1]

	# Artificially create missing samples (genotypes) in "input"
	missing_vector = np.random.choice(list(freq_model.keys()), p=list(freq_model.values()), size=len(truth_geno_df))
	total_missing = sum(missing_vector)
	missing_list = [np.random.choice(np.arange(samples), size=missing, replace=False) if missing != 0 else np.array([]) for missing in missing_vector]
	missing_genos = missing_vcf(truth_geno_df, missing_list)
	test_df = pd.concat([truth_df[truth_df.columns[:9]], missing_genos], axis=1)

	# This is for finding errors in input files/debugging only (should be run with --iters 1)
	if output:
		test_df.to_csv(output, index=False, sep="\t")

	# Impute missing genotypes
	tmp_missing = f"tmp_missing{tmp_id}.vcf"
	tmp_impute = f"tmp_impute{tmp_id}"
	test_df.to_csv(tmp_missing, index=False, sep="\t")
	subprocess.call(f"java -jar {BEAGLE} gt={tmp_missing} out={tmp_impute} > /dev/null", shell=True)

	# Compare imputed genotypes to true genotypes, return total missing, correctly imputed, and imputation accuracy
	subprocess.call("gunzip " + f"{tmp_impute}.vcf.gz", shell=True)
	if not phased:
		subprocess.call("sed -i 's/|/\//g' " + f"{tmp_impute}.vcf", shell=True)
		subprocess.call("sed -i 's/0\/1/1\/0/g' " + f"{tmp_impute}.vcf", shell=True)
	subprocess.call("sed -i '0,/#CHROM/{s/#CHROM/CHROM/}' " + f"{tmp_impute}.vcf", shell=True)
	impute_df = pd.read_csv(f"{tmp_impute}.vcf", dtype=str, sep="\t", comment="#", nrows=nrows)
	impute_df = impute_df.rename(columns={"CHROM": "#CHROM"})
	subprocess.call(f"rm {tmp_missing} {tmp_impute}.*", shell=True)

	impute_geno_df = impute_df[impute_df.columns[9:]]
	non_missing = total_genos - total_missing
	incorrect = (np.sum(truth_geno_df.compare(impute_geno_df).count())) // 2
	correct = total_missing - incorrect
	accuracy = (correct / total_missing) * 100
	stats = [correct, total_missing, accuracy]
	return stats


start = time.perf_counter()

parser = ag.ArgumentParser(description="Test imputation accuracy of BEAGLE")
parser.add_argument("truth", help="Path to the ground truth VCF (must have no missing samples)")
parser.add_argument("alleles", help="Path to the allele count file (generated via: vcftools --counts | cut -f4)")
parser.add_argument("--iters", dest="iters", metavar="<int>", type=int, default=1, help="Number of iterations to run for [default=1]")
parser.add_argument("--output", dest="output", metavar="<string>", default="out", help="Output prefix of the test result file, formatted as <PREFIX>.beagle_test.tabular [default=out]")
parser.add_argument("--phased", dest="phased_flag", action="store_true", help="Flag to use if data are phased")
parser.add_argument("--rows", dest="nrows", type=int, metavar="<int>", default=None, help="Subset all files to N rows. This is primarily for debugging [default=None]")
parser.add_argument("--test-output", dest="test_out", metavar="<string>", default="", help="Only specify if you wish to save the generated testing VCF (--iters 1 should be set). Primarily for debugging [default=None]")
parser.add_argument("--tmp-id", dest="tmp_id", metavar="<string>", default="", help="ID to append to tmp files. Use if running in parallel [default=None]")
args = parser.parse_args()


truth_df = truth_vcf(args.truth, args.nrows)
num_samples = len(truth_df.columns[9:])
model_obj = MissingModel()
freq_model = model_obj.load_model(args.alleles, num_samples, args.nrows)
print()

stats_dict = {
	"correct": [],
	"total_missing": [],
	"accuracy": []
}

total_start = time.perf_counter()
for i in range(args.iters):
	start = time.perf_counter()
	print(f"Start iter {i+1}...", end="", flush=True)
	
	stats = main(truth_df, freq_model, args.tmp_id, args.test_out, args.phased_flag, args.nrows)
	stats_dict["correct"].append(stats[0])
	stats_dict["total_missing"].append(stats[1])
	stats_dict["accuracy"].append(stats[2])

	delta = time.perf_counter() - start
	print(f"done in {delta:.2f}s")

total_delta = time.perf_counter() - total_start
print(f"Total run time: {total_delta:.2f}s")

accuracy = stats_dict["accuracy"]
mean = sum(accuracy) / len(accuracy)
std = math.sqrt(sum([(acc - mean)**2 for acc in accuracy]) / len(accuracy))
print(f"\nAccuracy mean: {mean:.2f}%")
print(f"Accuracy std: {std:.2f}%")

output = f"{args.output}.beagle_test.tabular"
print(f"\nWriting accuracy report to {output}", end="")
pd.DataFrame(stats_dict).to_csv(output, index=False, sep="\t")