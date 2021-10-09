from samples import Samples
import argparse as ag



parser = ag.ArgumentParser(description="Helper to vcf_col_fix.sh")
parser.add_argument("vcf_header", help="Path to the names of each sample from the VCF")
parser.add_argument("rename_key", help="Path to the relevant rename_key.csv")
parser.add_argument("tmp_dir", help="Path to the tmp dir")
args = parser.parse_args()

sample = Samples()
sample.vcf_header_fix(args.vcf_header, args.rename_key, f"{args.tmp_dir}/vcf_fixed.txt")