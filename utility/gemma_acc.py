import argparse as ag
import numpy as np
import os
import pandas as pd
from pathlib import Path
from scipy import stats


class Accuracy:
	def __init__(self, truth_path):
		self.truth = pd.read_csv(truth_path, sep="\t")
		self.clean_indices = self.truth.dropna().index
		self.truth = self.truth.loc[self.clean_indices][self.truth.columns[1]]


	# Definition: 
	# 			  Find values for a specific site, automatically ignoring phenotypes with an NA value
	# 			  in the truth DF. This is done by removing NA values in the truth DF (leaving only
	# 			  values that possessed phenotype information) and removing NAs in the pred DF (leaving
	# 			  only values which were either NA in the orignal dataset or values for the specific
	# 			  site). The intersect of these thus leaves only values which were both present in the
	# 			  original dataset and specific to the site.
	# Arguments:
	# 			  pred: a DF with the predicted values and with NAs already dropped
	# Returns:
	# 			  A numpy array of the site-specific, non-NA indices
	def get_indices(self, pred):
		site_indices = np.in1d(pred.index, self.clean_indices, assume_unique=True)
		return site_indices


	def regression(self):
		stats_dict = {}
		_, _, r, p, std_err = stats.linregress(self.truth, self.pred)
		stats_dict["r"] = r
		stats_dict["r^2"] = r**2
		stats_dict["p"] = p
		stats_dict["std_err"] = std_err
		return stats_dict



class Single(Accuracy):
	def __init__(self, truth_path, pred_path):
		super().__init__(truth_path)
		pred = pd.read_csv(pred_path, header=None)
		pred = pred.dropna()
		pred = pred.loc[self.get_indices(pred)]
		self.truth = self.truth[pred.index].reset_index(drop=True)
		self.pred = pred[0].reset_index(drop=True)



class All(Accuracy):
	def __init__(self, truth_path, pred_path):
		super().__init__(truth_path)
		dfs = []
		for path in Path(pred_path).rglob("*.prdt.txt"):
			pred = pd.read_csv(path, header=None)
			pred = pred.dropna()
			site_indices = self.get_indices(pred)
			dfs.append(pred.loc[site_indices])
		pred = pd.concat(dfs).sort_index()
		self.pred = pred[0].reset_index(drop=True)
		self.truth = self.truth.reset_index(drop=True)



parser = ag.ArgumentParser(description="Obtain the prediction accuracy of GEMMA")
parser.add_argument("truth", help="Path to the phenotype file used as input for the GEMMA.sh script (this should contain headers and a Microchip, phenotype, and Site column)")
parser.add_argument("prediction", help="Path to a file containing phenotype predicted by GEMMA or the top level directory containing multiple prediction files (must be .prdt.txt files)")
parser.add_argument("--site", "-s", dest="site", metavar="<string>", help="The site for which phenotype predictions were made. This must be used if specifying a single prediction file")
args = parser.parse_args()


if os.path.isfile(args.prediction):
	pheno_in = Single(args.truth, args.prediction)

else:
	pheno_in = All(args.truth, args.prediction)

for k, v in pheno_in.regression().items():
	print(f"{k}: {v}")