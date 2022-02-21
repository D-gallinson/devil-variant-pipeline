import argparse as ag
import pandas as pd
import vcf
import time
import subprocess
from io import StringIO
from collections import OrderedDict



# Sort the param file on descending gamma and grab the first N rows
def gamma(path, n, exclude):
	df = pd.read_csv(path, sep="\t")
	df = df.sort_values(by=["gamma"], ascending=False)
	df = df.iloc[:n]
	if exclude:
		df = df[~(df["rs"] == exclude)]
	return df


# Get genotypes from the mean genotype file
def geno_mg(gamma, geno_path):
	start = time.perf_counter()
	print(f"Searching for lines in geno MG...", end="", flush=True)
	targets = gamma["rs"]
	found = []
	print()
	for target in targets:
		snp = subprocess.check_output(f"zgrep -m1 \"{target}\" {geno_path}", shell=True).decode('utf-8').split("\n")
		found.append(snp[0].split())
	df = pd.DataFrame(found)
	df = df.rename(columns={0: "snpID", 1: "REF", 2: "ALT"})
	delta = time.perf_counter() - start
	print(f"done ({delta:.3f}s)")
	return df


def geno_pheno(geno_df, pheno_df):
	pheno_df = pheno_df.set_index("Microchip")
	pheno_df = pheno_df.iloc[:, 0]
	gamma = geno_df["gamma"].astype(float).round(3).astype(str)
	colnames = geno_df["snpID"] + " = " + gamma
	geno_df = geno_df.drop(columns=["snpID", "REF", "ALT", "gamma"])
	geno_df.columns = pheno_df.index
	geno_df = geno_df.transpose()
	geno_df.columns = colnames
	geno_df = geno_df.join(pheno_df)
	return geno_df


# Find VCF lines by looping through and zgrepping (inefficient, shouldn't be used with many samples)
def get_lines(targets, vcf_path):
	start = time.perf_counter()
	print(f"1) Searching for lines in VCF...", end="", flush=True)
	header = subprocess.check_output(f"zcat {vcf_path} | head -1000 | grep '^#'", shell=True).decode('utf-8')
	colnames = subprocess.check_output(f"zcat {vcf_path} | grep -v -m1 '^##'", shell=True).decode('utf-8').strip("\n")
	sample_names = colnames.split("\t")[9:]
	found = ""
	for i in range(targets.shape[0]):
		chrom = targets.iloc[i, 0]
		pos = targets.iloc[i, 1]
		search_str = f"{chrom}\t{pos}"
		found += subprocess.check_output(f"zgrep -m1 \"{search_str}\" {vcf_path}", shell=True).decode('utf-8')
	found = header + found
	delta = time.perf_counter() - start
	print(f"done ({delta:.3f}s)")
	return StringIO(found), sample_names


# Convert the found lines to a VCF with cols: CHROM POS REF ALT sample1..N
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
		for sample in record.samples:
			vcf_dict[sample.sample].append(sample["GT"])
	delta = time.perf_counter() - start
	df = pd.DataFrame(vcf_dict)
	print(f"done ({delta:.3f}s)")
	return df


# Write to a VEP-friendly format
def to_vep(vcf, vep_path):
	str_vcf = vcf.astype(str)
	vep_str = str_vcf["CHROM"] + " " + str_vcf["POS"] + " " + str_vcf["POS"] + " " + str_vcf["REF"] + "/" + str_vcf["ALT"] + " 1"
	vep_str = vep_str.to_string(header=False, index=False)
	vep_str = vep_str.split("\n")
	vep_str = [line.strip(" ") for line in vep_str]
	vep_str = "\n".join(vep_str)
	with open(vep_path, 'w') as f:
	    f.write(vep_str)
	print(f"Writing VEP file to {vep_path}")


def vcf_pheno(vcf_df, pheno_df):
	pheno_df = pheno_df.set_index("Microchip")
	pheno_df = pheno_df.iloc[:, 0]
	out_df = vcf_df.transpose()
	chrom = out_df.loc["CHROM"]
	pos = out_df.loc["POS"].astype(str)
	gamma = out_df.loc["gamma"].astype(float).round(3).astype(str)
	colnames = [chrom[i] + ":" + pos[i] + " = " + gamma[i] for i in range(len(chrom))]
	out_df.columns = colnames
	out_df = out_df.drop(index=["CHROM", "POS", "REF", "ALT", "gamma"])
	out_df = out_df.join(pheno_df)
	return out_df


# CLI arguments
parser = ag.ArgumentParser(description="Obtain gamma values from a GEMMA output.param file. Optionally, the genotypes of the top N gamma SNPs can be obtained")
parser.add_argument("param_path", help="Path to param file (can be gzipped)")
parser.add_argument("--exclude", dest="exclude", metavar="<string>", default="", help="SNP to exclude based on rs ID (currently only supports excluding a single SNP) [default=None]")
parser.add_argument("--geno", dest="geno_path", metavar="<string>", default="", help="Path to a BIMBAM mean genotype file (this can be used instead of --vcf) [default=None]")
parser.add_argument("-n", dest="n_gamma", metavar="<int>", type=int, default=5, help="Number of top SNPs to use from the param file [default=5]")
parser.add_argument("--output", dest="output", metavar="<string>", default="gamma.vcf", help="Output path to write the found VCF SNPs to [default=gamma.vcf]")
parser.add_argument("--pheno", dest="pheno_path", metavar="<string>", default="", help="Path to phenotype file. Specifying this transposes the matrix to facilitate plotting [default=None]")
parser.add_argument("--print", dest="print_flag", action="store_true", help="Flag to print the top N gamma SNPs")
parser.add_argument("--vcf", dest="vcf_path", metavar="<string>", default="", help="Path to the VCF file (this can be used instead of --geno) [default=None]")
parser.add_argument("--vep", dest="vep_path", metavar="<string>", default="", help="Path to output in a VEP-friendly format [default=None]")
args = parser.parse_args()


gamma_df = gamma(args.param_path, args.n_gamma, args.exclude)

if args.print_flag:
	print(gamma_df.to_string(index=False))
	exit()

if args.vep_path and not args.vcf_path:
	print("Error! A VCF path (--vcf) must be specified if outputting VEP_path format. Please re-run and specify a VCF path. Exiting.")
	exit()

if args.vcf_path:
	search = gamma_df["rs"].str.split(":", expand=True)
	file_like, names = get_lines(search, args.vcf_path)
	vcf = vcf_df(file_like, names)
	if args.vep_path:
		to_vep(vcf, args.vep_path)
	vcf.insert(loc=4, column="gamma", value=gamma_df["gamma"].reset_index(drop=True))
	index = False
	if args.pheno_path:
		pheno_df = pd.read_csv(args.pheno_path, sep="\t")
		vcf = vcf_pheno(vcf, pheno_df)
		index = True

	print(f"Writing found VCF lines to {args.output}")
	vcf.to_csv(args.output, index=index)
else:
	geno_df = geno_mg(gamma_df, args.geno_path)
	index = False
	if args.pheno_path:
		geno_df["gamma"] = gamma_df["gamma"].reset_index(drop=True)
		pheno_df = pd.read_csv(args.pheno_path, sep="\t")
		geno_df = geno_pheno(geno_df, pheno_df)
		index = True
	print(f"Writing found geno MG lines to {args.output}")
	geno_df.to_csv(args.output, index=index)