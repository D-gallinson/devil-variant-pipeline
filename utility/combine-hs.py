#Combines alignment statistics from multiple CollectHsMetrics files
import os
import pandas as pd
import argparse as ag

parser = ag.ArgumentParser(description="Rename files specified by a CSV")
parser.add_argument("input", help="Path to the input HS directory")
parser.add_argument("output", help="Path to the write the output csv to")
args = parser.parse_args()


input = args.input
files = os.listdir(input)
files = [f"{input}/{file}" for file in files]

metrics = {}
has_cols = False

for file in files:
	with open(file) as handle:
		for line in handle:
			#Search for the line with "METRIC CLASS". The following line is the header and the line following that is the alignment stats
			if line.find("METRICS CLASS") != -1:
				#For the first file, generate the dictionary keys
				if not has_cols:
					metrics = {col: [] for col in handle.readline().rstrip().split()}
					has_cols = True
				else:
					#Eat the stats header for subsequent files beyond the first
					handle.readline()
				rowvals = handle.readline().rstrip().split()
				#Loop through the stats and assign each stat to the corresponding header.
				#NOTE: if there are missing statistic values (except at the end) then all subsequent values will be assigned to the WRONG stats header
				#	   If a column in the summary CSV has odd or outlier values, the original metrics file should be inspected
				for i in range(len(metrics)):
					key = list(metrics.keys())[i]
					if i < len(rowvals):
						metrics[key].append(rowvals[i])
					else:
						#If there are missing stats values, they are appended as NaN in the final columns
						metrics[key].append("NaN")
				break

metrics_df = pd.DataFrame(metrics)
#Change cols to extract additional/fewer columns
cols = ["PCT_SELECTED_BASES", "FOLD_ENRICHMENT", "MEDIAN_TARGET_COVERAGE", \
		"PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_40X"]
relevant_metrics = metrics_df[cols]
relevant_metrics.to_csv(args.output, index=False)