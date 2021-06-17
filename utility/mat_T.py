import numpy as np
import argparse as ag


parser = ag.ArgumentParser(description="Transpose a matrix")
parser.add_argument("mode", help="Directory path for matrix file [SH=SNP-host, SP=SNP-pathogen, IH=indel-host, IP=indel-pathogen]")
args = parser.parse_args()

if args.mode.lower() == "sh":
	type = "SNP-host"
elif args.mode.lower() == "sp":
	type = "SNP-pathogen"
elif args.mode.lower() == "ih":
	type = "indel-host"
elif args.mode.lower() == "ip":
	type = "indel-pathogen"
else:
	print(f"Erroneous argument \"{args.mode}\", exiting. Please use: [SH=SNP-host, SP=SNP-pathogen, IH=indel-host, IP=indel-pathogen]")
	exit()

file = f"geno_{type}"
root = "/work_bgfs/d/dgallinson/outputs/intermediates/8_joint-variants"

mat = np.loadtxt(f"{root}/{file}.012", dtype=int)
print(f"Starting shape: {mat.shape}")
T = mat.T
print(f"Ending shape: {T.shape}")
np.savetxt(f"{root}/{file}.T.012", T, fmt="%i", delimiter="\t")