import numpy as np
import argparse as ag


parser = ag.ArgumentParser(description="Transpose a matrix")
parser.add_argument("file", help="Path to the file. The file must be a matrix to be tranposed")
args = parser.parse_args()

mat = np.loadtxt(args.file, dtype=int)
print(f"Starting shape: {mat.shape}")
T = mat.T
print(f"Ending shape: {T.shape}")
print(f"Saved to {args.file}.T")
np.savetxt(f"{args.file}.T", T, fmt="%i", delimiter="\t")