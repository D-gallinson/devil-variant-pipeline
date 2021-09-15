import numpy as np
import subprocess
import argparse as ag


parser = ag.ArgumentParser(description="Swap minor and major allele dosages in a mean genotype file")
parser.add_argument("input", help="Path to the mean genotype file")
parser.add_argument("-a", "--all", dest="all_flag", action="store_true", help="Flag to use if all alleles should be swapped (this overrides the freq file set with -f)")
parser.add_argument("-d", "--delim", dest="delim", metavar="<string>", default=" ", help="Mean genotype file delimiter [default=\" \"]")
parser.add_argument("-f", "--freq", dest="freq_path", metavar="<file>", default="input.frq", help="Path to file containing REF allele frequencies. Generate with \"VCFtools --freq2 | cut -f 5 | tail -n +2\".[default=None]")
parser.add_argument("-o", "--output", dest="output", metavar="<file>", default="output.mg", help="Name of the output file [default=output.mg]")
args = parser.parse_args()


with open(args.input) as h:
	ncols = len(h.readline().split(args.delim))
data = np.genfromtxt(args.input, delimiter=args.delim, usecols=range(3, ncols))

if args.all_flag:
	data = (data * -1) + 2
	print(f"Swapping all {len(data)} sites")
else:
	freqs = np.genfromtxt(args.freq_path)
	swap_indices = freqs < 0.5
	print(f"Sites swapped: {np.count_nonzero(swap_indices)}")
	data[swap_indices] = (data[swap_indices] * -1) + 2

np.savetxt("allele_tmp_mg", data, fmt="%1.f", delimiter=" ")
if args.delim == " ":
	args.delim = "\" \""
subprocess.call("sed -i \"s/nan/NA/g\" allele_tmp_mg", shell=True)
subprocess.call(f"cut -d {args.delim} -f 1-3 {args.input} > allele_tmp_mg_info", shell=True)
subprocess.call(f"paste -d {args.delim} allele_tmp_mg_info allele_tmp_mg > {args.output}", shell=True)
subprocess.call("rm allele_tmp*", shell=True)