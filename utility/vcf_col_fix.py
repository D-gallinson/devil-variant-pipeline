from samples import Samples
import argparse as ag
import pandas as pd
import re



parser = ag.ArgumentParser(description="Helper to vcf_col_fix.sh")
parser.add_argument("vcf_header", help="Path to the names of each sample from the VCF")
parser.add_argument("outpath", help="Output path for the tmp file")
args = parser.parse_args()


sample = Samples()
sample.subset("Tissue", "Tumour")
haystack = sample.to_vcf_chip(inplace=False).dropna().to_list()

target_list = pd.read_csv(args.vcf_header, header=None, squeeze=True)
target_list[target_list == "T2-426081"] = "T4-426081" # This fixes a problem with the batch 1 tumors where the L0200 tumor was named as T2 instead of T4
target_list = [target if re.match("T\d", target) else target[-6:] for target in target_list]

corrected_names = []
for target in target_list:
	if re.match("T\d", target):
		corrected_names.append(target)
		continue
	for mchip in haystack:
		if re.match(f"T.*-{target}", mchip):
			corrected_names.append(mchip)
			break

output = pd.Series(corrected_names)
output.to_csv(args.outpath, header=False, index=False)