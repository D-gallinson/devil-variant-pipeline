import os
import re
import time
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import statsmodels.formula.api as sm



class Samples:
	base = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data"
	seq_base = f"{base}/data/sequencing"
	batch_ids = [f"{seq_base}/Capture1/rename_key.csv", f"{seq_base}/Capture2/rename_key.csv", f"{seq_base}/Capture3/rename_key.csv", f"{seq_base}/Capture4/rename_key.csv", f"{seq_base}/Capture5/rename_key.csv"]
	prelim = batch_ids[:2]
	dftd_site_arrival = pd.DataFrame({
		"Site": ["Mt William", "Freycinet", "WPP", "Wilmot", "West Takone", "Takone", "Black River", "Arthur River"],
		"arrival_date": [1996, 1999, 2007, 2008, 2013, 2013, 2015, 2019],
		"generations": [12.5, 11, 7, 6.5, 4, 4, 3, 1],
		"years": [25, 22, 14, 13, 8, 8, 6, 2]
	})
	dftd_site_arrival["arrival_date"] = pd.to_datetime(dftd_site_arrival["arrival_date"], format="%Y")


	def __init__(self,
		sample_csv_path=f"{base}/data/pheno_data/master_corrections.csv",
		id_paths=[],
		assume_YOB=True
		):
		dates = ["TrappingDate", "YOB"]
		full_df = pd.read_csv(sample_csv_path)
		if assume_YOB:
			full_df = self.__assume_YOB(full_df)
		full_df[dates] = full_df[dates].apply(pd.to_datetime)
		if id_paths:
			ids = self.__L_ID(id_paths)
			self.sample_df = full_df[full_df["Library number"].isin(ids)].reset_index(drop=True)
		else:
			self.sample_df = full_df
		self.factor_key = {"pheno": [], "key": [], "val": []}


	@staticmethod
	def dynamic_groupby(df, group_col, agg_col, func):
		current_funcs = ["min", "max"]
		groups = df.groupby(group_col)[agg_col]
		if func == "min":
			groups = groups.min()
		elif func == "max":
			groups = groups.max()
		else:
			print(f"Error, func {func} does not exist! Please use one of: {','.join(current_funcs)}")
			return None
		group_dict = {"group": groups.index, "agg": groups.values}
		rows = []
		for i in range(len(group_dict["group"])):
			current_group = group_dict["group"][i]
			current_agg = group_dict["agg"][i]
			rows.append(df.loc[(df[group_col] == current_group) & (df[agg_col] == current_agg)])
		return pd.concat(rows)


	def __str__(self):
		return repr(self.sample_df)


	def add_arrival(self, arrival_cols=["arrival_date"], inplace=True, WPP_tetraploid=True):
		if WPP_tetraploid:
			self.alter_arrival("WPP", 2011)
		arrivals = Samples.dftd_site_arrival[["Site"] + arrival_cols]
		df = self.sample_df.merge(arrivals, how="left", on="Site")
		if not inplace:
			return df
		self.sample_df = df


	def add_col(self, col):
		name = col.name if col.name else "new_col"
		self.sample_df[name] = col


	def alter_arrival(self, site, new_year, current_year=2021):
		i = Samples.dftd_site_arrival[Samples.dftd_site_arrival["Site"] == site].index.values
		Samples.dftd_site_arrival.loc[i, "arrival_date"] = pd.to_datetime(new_year, format="%Y")
		Samples.dftd_site_arrival.loc[i, "years"] = current_year - new_year
		Samples.dftd_site_arrival.loc[i, "generations"] = (current_year - new_year) // 2


	def arrival_age(self, WPP_tetraploid=True):
		if WPP_tetraploid:
			self.alter_arrival("WPP", 2011)
		self.add_arrival(WPP_tetraploid=WPP_tetraploid)
		self.sample_df["arrival_age"] = (self.sample_df["arrival_date"] - self.sample_df["YOB"]).dt.days
		self.sample_df.loc[self.sample_df["arrival_age"] < 0, "arrival_age"] = 0
		self.sample_df = self.sample_df.drop(columns=["arrival_date"])


	def chip_sort(self, df, num_only=False, reset_index=True):
		if num_only:
			sort_df = df.sort_values(by="Microchip", key=lambda x: x.str[-6:])
		else:
			sort_df = df.sort_values(by="Microchip")
		if reset_index:
			sort_df.reset_index(drop=True, inplace=True)
		return sort_df


	def compare_to(self, df_path, cols, sep=",", id_col="Microchip", inplace=True):
		compare_df = pd.read_csv(df_path, sep=sep)
		compare_df = compare_df[[id_col] + cols]
		compare_df.rename(columns={id_col: "Microchip"}, inplace=True)
		combined = self.sample_df.merge(compare_df, how="left", on="Microchip")
		if inplace:
			self.sample_df = combined
		else:
			return combined


	def count(self, col, pattern):
		return len(self.sample_df[self.sample_df[col] == pattern])


	def count_groups(self, col, sorted=True):
		groups = self.sample_df.groupby(col)[col].count()
		return groups.sort_index()


	def date_only(self, col):
		self.sample_df[col] = self.sample_df[col].dt.date


	# === COVARIATE ===
	# Definition: 
	# 			  A nuanced version of DFTD_YOB, this counts the number of generations a devil was born *after* DFTD arrived at
	# 			  its site. Devils born at or before arrival get a value of 0.
	# Arguments:
	# 			  gen_time:		  devil generation time, used to convert YOB days after arrival to YOB generations after arrival
	# 			  WPP_tetraploid: add labels to delineate WPP samples born during the presence of the tetraploid variant
	# Returns:
	# 			  None. Adds a "DFTD_generations" column to sample_df
	# Steps:
	# 			  1) Use 2011 as WPP arrival (if WPP_tetraploid=True) and add arrival dates
	# 			  2) Get the number of days each devil was born after DFTD arrival
	# 			  3) Set devils born before arrival to a generation time of 0
	# 			  4) Convert devils born after arrival to generations after arrival
	def DFTD_generations(self, gen_time=2, WPP_tetraploid=True):
		if WPP_tetraploid:
			self.alter_arrival("WPP", 2011)
		self.add_arrival()
		self.sample_df["DFTD_generations"] = (self.sample_df["YOB"] - self.sample_df["arrival_date"]).dt.days
		self.sample_df.loc[self.sample_df["DFTD_generations"] < 0, "DFTD_generations"] = 0
		self.sample_df["DFTD_generations"] = self.sample_df["DFTD_generations"] / (365 * gen_time)
		self.sample_df.drop(columns="arrival_date", inplace=True)


	# === COVARIATE ===
	# Definition: 
	# 			  Add a categorical covariate to samples of whether a devil was born or after DFTD arrived at its respective
	# 			  site. This is meant to capture devils born in the presence of DFTD selective pressures (post-DFTD) or pre-
	# 			  DFTD selective pressures. Individuals born directly at DFTD arrival are considered pre-DFTD, as their parents
	# 			  experienced no DFTD selective pressures. See the "dftd_site_arrival" DF for arrival times.
	# Arguments:
	# 			  WPP_tetraploid: add labels to delineate WPP samples born during the presence of the tetraploid variant
	# Returns:
	# 			  None. Adds a "born_arrival" column to sample_df
	# Steps:
	# 			  1) Use 2011 as WPP arrival (if WPP_tetraploid=True) and add arrival dates
	# 			  2) Obtain indices for devils born before DFTD arrival and after arrival
	# 			  3) Set these samples to have the proper pre/post-DFTD label and drop the arrival_date col
	# 			  4) If include_tetraploid=True, use "tetraploid" for applicable WPP samples
	# Gotchas:
	# 			  Devils born very near to DFTD arrival (e.g., months after) may not have been affected by DFTD selection.
	# 			  Using "DFTD_generations()" may provide a more nuanced approach which better represents a continuum of
	# 			  selective pressure (e.g., devils born 1 month following arrival are less affected by selection than those
	# 			  born 5 generations after)
	def DFTD_YOB(self, WPP_tetraploid=True, include_tetraploid=True):
		if WPP_tetraploid:
			self.alter_arrival("WPP", 2011)
		self.add_arrival()
		pre_dftd = self.sample_df[self.sample_df["YOB"] <= self.sample_df["arrival_date"]].index
		post_dftd = self.sample_df[self.sample_df["YOB"] > self.sample_df["arrival_date"]].index
		self.sample_df.loc[pre_dftd, "born_arrival"] = "pre-DFTD"
		self.sample_df.loc[post_dftd, "born_arrival"] = "post-DFTD"
		self.sample_df.drop(columns="arrival_date", inplace=True)
		if include_tetraploid:
			self.WPP_tetraploid()


	def duplicates(self):
		cluster = self.sample_df.groupby("Microchip")["Microchip"].count()
		dup_indices = cluster[cluster > 1].index.values
		dups = self.sample_df[self.sample_df["Microchip"].isin(dup_indices)]
		print(f"Duplicate microchips: {dups.shape[0]}/{self.sample_df.shape[0]} ({dups.shape[0]/self.sample_df.shape[0]*100:.2f}%)")
		return dups


	# A sanity check to ensure that entries with "True" in the "Host-tumour pair" column
	# align with what's obtained using the host_tumor_pairs() method below.
	# NOTE: unless all 1056 samples are loaded, this will likely differ from the result
	# 		obtained by host_tumor_pairs(), as it is likely only a single member of the
	# 		pair was sequenced
	def expected_pairs(self):
		first_tissue = self.sample_df["Tissue"].unique()[0]
		# expected_pairs = self.subset("Tissue", first_tissue, inplace=False)
		expected_pairs = self.sample_df[self.sample_df["Host-tumour pair"] == True]
		return expected_pairs


	# Parameter lib_list can be a list of L0123 IDs or a string range as: L0960-L1056
	def extract_libnum(self, lib_list, inplace=True):
		if isinstance(lib_list, str):
			lib_list = lib_list.replace(" ", "")
			lib_list = lib_list.replace("L", "")
			lib_list = lib_list.split("-")
			start = int(lib_list[0])
			end = int(lib_list[1]) + 1
			lib_list = [f"L{i:04d}" for i in range(start, end)]
		extracted = self.sample_df[self.sample_df["Library number"].isin(lib_list)]
		if not inplace:
			return extracted
		self.sample_df = extracted


	def extract_mchips(self, chip_list, inplace=True):
		found = []
		for chip in chip_list:
			last_six = chip[-6:]
			samples = self.sample_df[self.sample_df["Microchip"].str[-6:] == last_six]
			if chip[0] == "H" or chip[0] == "T":
				tissue = "Host" if chip[0] == "H" else "Tumour"
				samples = samples[samples["Tissue"] == tissue]
			found.append(samples)
		extracted = pd.concat(found)
		print(f"Found {extracted.shape[0]} of {len(chip_list)} samples")
		if not inplace:
			return extracted
		self.sample_df = extracted


	# Document fully. This is really only for writing to a file because it converts any date to a str
	def extract_year(self, col, inplace=True, beginning=True):
		dates = self.sample_df[col]
		years = dates[~dates.isna()].astype(str)
		if beginning:
			years = years.str[:4]
		else:
			years = years.str[-4:]
		nans = dates[dates.isna()]
		years = pd.concat([years, nans]).sort_index()
		if not inplace:
			return years
		self.sample_df[col] = years


	def factor_file(self, path, cols=[], silent=False):
		self.update_factors(cols)
		factor_df = pd.DataFrame(self.factor_key)
		if not cols:
			subset = factor_df
		else:
			subset = factor_df[factor_df["pheno"].isin(cols)]
		if subset.empty:
			if not silent:
				print("No factors to print")
			return None
		print(f"Writing factor file to: {path}")
		subset.to_csv(path, index=False, sep="\t")


	def fast_stats(self):
		pairs = self.host_tumor_pairs(inplace=False)["Microchip"]
		singles = self.sample_df.shape[0] - len(pairs)
		multi = len(self.host_tumor_multi()["Microchip"])
		print(f"Samples: {self.sample_df.shape[0]}")
		print(f"Hosts: {self.sample_df[self.sample_df['Tissue'] == 'Host'].shape[0]}")
		print(f"Tumors: {self.sample_df[self.sample_df['Tissue'] == 'Tumour'].shape[0]}")
		print(f"Missing microchips: {self.nan_sum('Microchip', print_flag=False)}")
		print(f"Unique microchips: {len(self.sample_df['Microchip'].unique())}")
		print(f"Singles: {singles}")
		print(f"Pairs (individuals): {len(pairs)}")
		print(f"Pairs: {len(pairs.unique())}")
		print(f"Multi tumor/host: {multi}")


	def get_col(self, col):
		return self.sample_df[col]


	# col can be a string or list of strings for clustering on multiple columns
	def get_groups(self, col):
		return self.sample_df.groupby(col)


	# LIKELY DEPRECATE
	# This considers Pieldivina and Iota to be paired, despite them having no tumor samples
	# Look into extracting pairs only by tumor/host (perhpas get pairs, do something with Host-tumour pair col)
	# __handle_pairs might also be better suited here
	def get_pairs(self, pair_count=None):
		# if pairs_only:
		# 	pairs = self.subset("Host-tumour pair", True, inplace=False)
		pairs = self.sample_df.groupby("Microchip")["Microchip"].count()
		if pair_count:
			pairs = pairs[pairs == pair_count]
		return pairs


	def host_tumor_pairs(self, inplace=True):
		pairs = self.sample_df.groupby("Microchip")["Tissue"].nunique()
		pairs = self.sample_df[self.sample_df["Microchip"].isin(pairs[pairs > 1].index.values)]
		if not inplace:
			return pairs
		print(f"*host_tumor_pairs()* Removing {self.sample_df.shape[0] - pairs.shape[0]} rows out of {self.sample_df.shape[0]}  ({pairs.shape[0]} remaining)")
		self.sample_df = pairs


	def host_tumor_multi(self):
		pairs = self.host_tumor_pairs(inplace=False)
		multi = pairs.groupby(["Microchip", "Tissue"])["Microchip"].count()
		multi_chips = multi[multi > 1].index.get_level_values("Microchip")
		return self.sample_df[self.sample_df["Microchip"].isin(multi_chips)].sort_values(by=["Microchip", "Tissue"])


	def nan_rows(self, col, tissue="both"):
		if tissue == "both":
			subset = self.sample_df
		else:
			subset = self.sample_df[self.sample_df["Tissue"] == tissue]
		return subset[subset[col].isna()]


	def nan_sum(self, col, tissue="both", print_flag=True):
		if tissue == "both":
			subset = self.sample_df
		else:
			subset = self.sample_df[self.sample_df["Tissue"] == tissue]
		total = len(subset)
		nans = subset[col].isna().sum()
		if print_flag and isinstance(nans, np.integer):
			print(f"{col}: {total-nans}/{total} ({(total-nans)/total*100:.2f}%) non-NaNs")
			print(f"{col}: {nans}/{total} ({nans/total*100:.2f}%) NaNs")
		else:
			return nans


	# Definition: 
	# 			  Functions to normalize and transform data with sample_df. These can be strung together to transform and then normalize data, such as
	# 			  the log normal transform: samples.normalize("my_col", method="log"); samples.normalize("my_col"). If the original data within a col
	# 			  wish to be maintained, the new_col argument should be used, with subsequent normalizations/transformations targetting that new col.
	# 			  Normalized cols are named as new cols based on the ordering of "cols" and "new_cols", e.g: cols=["col1", "col2", "col3"] and
	# 			  new_cols=["col1_log", "col2_norm", "col3_lognorm"]
	# Arguments:
	# 			  cols: 	the columns to be normalized <string|list>
	# 			  method: 	method of normalization/transformation <string>
	# 			  inplace: 	whether or not to return data or change the column in place <boolean>
	# 			  new_cols: put the normalized/transformed data into new columns <list>
	# Returns:
	# 			  The normalized/transformed column(s) or None if inplace or new_cols is specified
	# Steps:
	# 			  Depends on the chosen method
	# Gotchas:
	# 			  inplace versus new_cols: if both are specified, new_cols will override inplace=True
	# TODO:
	# 			  Test all methods by hand to ensure the norms/transforms are working.
	# 			  Consider an "undo" method (maintaining raw values with new_cols may be sufficient)
	def normalize(self, cols, method="standard", inplace=True, new_cols=[]):
		methods = ["inv_logit", "log", "logit", "standard", "ranknorm"]
		cols = cols if isinstance(cols, list) else [cols]
		df_cols = self.sample_df[cols]
		if method == "standard":
			means = df_cols.mean()
			stds = df_cols.std()
			norms = (df_cols - means) / stds
		elif method == "inv_logit":
			norms = np.exp(df_cols) / (1 + np.exp(df_cols))
		elif method == "log":
			norms = np.log(df_cols)
		elif method == "logit":
			norms = np.log(df_cols / (1 - df_cols))
		elif method == "ranknorm":
			k = 0.375
			n = df_cols.shape[0]
			norms = df_cols.rank()
			indices = norms.index
			norms = norm.ppf((norms - k) / (n - 2 * k + 1))
			norms = pd.DataFrame(norms, columns=cols, index=indices)
		else:
			print(f"*normalize() ERROR* Method \"{method}\" does not exist! Use: {', '.join(methods)}")
			return None
		if not inplace:
			return norms
		elif new_cols:
			df_cols = new_cols
		else:
			df_cols = cols
		for i in range(len(df_cols)):
			self.sample_df[df_cols[i]] = norms[cols[i]]


	def pair_rows(self, pair_count):
		pairs = self.get_pairs()
		pairs = pairs[pairs == pair_count]
		chips = pairs.index.values
		return self.sample_df[self.sample_df["Microchip"].isin(chips)]


	def plot_col(self, col, outpath, plot="hist", bins=10, xlab=None, ylab=None, title=None):
		plots = ["hist", "bar"]
		data = self.sample_df[col]
		title = col if title == None else title
		xlab = xlab if xlab == None else col
		if plot == "hist":
			ylab = ylab if ylab else "Frequency"
			if bins == "Freedman–Diaconis":
				q25, q75 = np.percentile(data, [0.25, 0.75])
				bin_width = 2 * ((q75 - q25) / len(data)**(-1/3))
				bins = round((data.max() - data.min()) / bin_width)
			plt.hist(data, bins=bins, density=True)
		elif plot == "bar":
			ylab = ylab if ylab else "Count"
			counts = data.value_counts()
			x = counts.index.values
			heights = counts.values
			plt.bar(x, heights)
		else:
			print(f"*plot_col() ERROR* Plot \"{plot}\" does not exist! Use: {', '.join(plots)}")
			return None
		plt.xlabel(xlab)
		plt.ylabel(ylab)
		plt.title(title)
		plt.savefig(outpath)
		plt.close()
		print(f"Plot saved to: {outpath}")


	def print_col(self, col):
		print(self.sample_df[col].to_list())


	def print_cols(self, cols):
		print(self.sample_df[cols])


	def print_factor_key(self):
		if not self.factor_key["key"]:
			print("No factors have been created")
			return None
		factor_df = pd.DataFrame(self.factor_key)
		for pheno in factor_df["pheno"].unique():
			current = factor_df[factor_df["pheno"] == pheno]
			print(pheno)
			print(current[["key", "val"]].to_string(index=False))


	def subset(self, col, pattern, inplace=True, verbose=True):
		subset = self.sample_df[self.sample_df[col] == pattern]
		if inplace and verbose:
			print(f"*subset()* Removing {self.sample_df.shape[0] - subset.shape[0]} rows out of {self.sample_df.shape[0]} ({subset.shape[0]} remaining)")
		if not inplace:
			return subset
		self.sample_df = subset


	def subset_arrival(self, WPP_tetraploid=True, inplace=True):
		tumors = self.sample_df[self.sample_df["Tissue"] == "Tumour"].size > 0
		if tumors:
			tumor_df = self.subset("Tissue", "Tumour", inplace=False)
		if WPP_tetraploid:
			self.alter_arrival("WPP", 2011)
		self.add_arrival()
		subset = self.subset("Tissue", "Host", inplace=False)
		subset = subset.dropna(subset=["YOB", "Site"])
		subset["date_diff"] = (subset["YOB"] - subset["arrival_date"]).dt.days
		final = subset[subset["date_diff"] >= 0]
		if tumors:
			final = pd.concat([final, tumor_df])
			final = final.sort_index()
		final = final.drop(columns=["arrival_date", "date_diff"])
		if not inplace:
			return final
		print(f"*subset_arrival()* Removing {self.sample_df.shape[0] - final.shape[0]} rows out of {self.sample_df.shape[0]} ({final.shape[0]} remaining)")
		self.sample_df = final


	def subset_pairs(self):
		indices = self.get_pairs()
		indices = indices[indices > 1].index
		self.sample_df = self.sample_df[self.sample_df["Microchip"].isin(indices)]


	def subset_non_dup(self, keep="first"):
		start = self.sample_df.shape[0]
		self.sample_df = self.sample_df[~self.sample_df.duplicated(subset=["Microchip"], keep=keep)]
		print(f"*subset_non_dup()* Removing {start - self.sample_df.shape[0]} out of {start} duplicate Microchips ({self.sample_df.shape[0]} remaining)")


	def subset_non_nan(self, col, verbose=True):
		na = self.sample_df[col].isna()
		if isinstance(col, list):
			na = na.any(axis=1)
		if verbose:
			print(f"*subset_non_nan()* Removing {na.values.sum()} rows out of {self.sample_df.shape[0]} ({self.sample_df.shape[0] - na.values.sum()} remaining)")
		self.sample_df = self.sample_df[~na]


	def summarize_pairs(self):
		pairs = self.get_pairs()
		print("PAIR\tCOUNT")
		for pair in pd.unique(pairs):
			print(f"{pair}\t{len(pairs[pairs == pair])}")


	def to_factor(self, col, inplace=True):
		phenotype = self.sample_df[col]
		unique_vals = sorted(phenotype.dropna().unique(), key=str.casefold)
		phenotype = phenotype.fillna("NA")
		factors = [i for i in range(len(unique_vals))]
		phenotype = phenotype.replace(unique_vals, factors)
		self.factor_key["pheno"] += [col for i in range(len(unique_vals))]
		self.factor_key["key"] += list(unique_vals)
		self.factor_key["val"] += factors
		if inplace:
			self.sample_df[col] = phenotype
		else:
			return phenotype


	def ATOMM_pheno(self, outpath, cols, to_pheno=True, to_ATOMM=False):
		tissue = pd.unique(self.sample_df["Tissue"])
		if len(tissue) < 2:
			print(f"*WARNING* the sample DF only contains a single tissue ({tissue[0]}), and thus pairs cannot be generated.")
			print("If this is a TumorSamples object, please specify the argument tissue=\"both\"")
			print("Exiting")
			return None
		host = self.sample_df[self.sample_df["Tissue"] == "Host"]
		if host[cols].isna().any(axis=None):
			print("*WARNING* the sample DF contains NA values in some of the columns.")
			print("If this is unintentional, run: samples.subset_non_nan(my_cols) before writing to a phenotype file")
		self.host_tumor_pairs()
		self.__handle_triplets()
		self.to_vcf_chip()
		tumor = self.chip_sort(self.sample_df[self.sample_df["Tissue"] == "Tumour"], num_only=True)
		self.subset("Tissue", "Host")
		host = self.chip_sort(self.sample_df, num_only=True)
		host.rename(columns={"Microchip": "host_chip"}, inplace=True)
		host["tumor_chip"] = tumor["Microchip"]
		self.sample_df = host
		if to_pheno and not to_ATOMM:
			self.to_pheno(outpath, cols, id_cols=["host_chip", "tumor_chip"], vcf_chip=False)
		if to_ATOMM:
			ids = [i for i in range(1, self.sample_df.shape[0] + 1)]
			intercept = [1 for i in range(1, self.sample_df.shape[0] + 1)]
			ATOMM_df = pd.DataFrame({"host": ids, "patho": ids, "inter": intercept})
			ATOMM_df = pd.concat([ATOMM_df, self.sample_df[cols]], axis=1)
			print(f"Writing ATOMM phenotype file to: {outpath}")
			ATOMM_df.to_csv(outpath, index=False, header=False, sep=" ")


	def to_pheno(self, outpath, cols, id_cols=["Microchip"], vcf_chip=True, sep="\t"):
		if vcf_chip:
			self.to_vcf_chip()
		pheno_df = self.sample_df[id_cols + cols]
		pheno_df = pheno_df.fillna("NA")
		outdir = outpath[:outpath.rfind(".")]
		factor_out = f"{outdir}_FACTOR_KEY.txt"
		self.factor_file(factor_out, cols=cols, silent=True)
		print(f"Writing phenotype file to: {outpath}")
		pheno_df.to_csv(outpath, index=False, sep=sep)


	def to_vcf_chip(self, inplace=True):
		chips = self.sample_df["Microchip"]
		chips = chips.str[-6:]
		tissue = self.sample_df["Tissue"].str[0]
		tissue[~self.sample_df["TumourNumber"].isna()] += self.sample_df["TumourNumber"]
		chips = tissue + "-" + chips	
		if inplace:
			self.sample_df["Microchip"] = chips
		else:
			return chips


	def trapping_age(self, trap_tables=[], max=True):
		if not trap_tables:
			self.sample_df["trap_age"] = (self.sample_df["TrappingDate"] - self.sample_df["YOB"]).dt.days
		else:
			df = pd.concat([pd.read_csv(table, usecols=["Microchip", "TrappingDate"]) for table in trap_tables])
			found_traps = df[df["Microchip"].isin(self.sample_df["Microchip"])]
			found_traps = found_traps[~found_traps.isna()]
			found_traps["TrappingDate"] = pd.to_datetime(found_traps["TrappingDate"])
			if max:
				extreme_trap = found_traps.groupby("Microchip")["TrappingDate"].max()
				col_name = "max_trap_age"
			else:
				extreme_trap = found_traps.groupby("Microchip")["TrappingDate"].min()
				col_name = "min_trap_age"
			extreme_trap.name = "extreme_trap"
			self.sample_df = self.sample_df.merge(extreme_trap, on="Microchip", how="left")
			self.sample_df.loc[self.sample_df["extreme_trap"].isna(), "extreme_trap"] = self.sample_df["TrappingDate"]
			self.sample_df[col_name] = (self.sample_df["extreme_trap"] - self.sample_df["YOB"]).dt.days
			self.sample_df.drop(columns="extreme_trap", inplace=True)


	# Determine total number of times trapped (irrespective of DFTD status)
	# Loop defaulting to the two paths is a bit hacky
	# Also, a TrapDate/TrappingDate general method would be nice to deal with this in the future
	def trap_count(self, path="../../data/pheno_data/originals/captData.csv"):
		default_paths = ["../../data/pheno_data/originals/captData.csv", "../../data/pheno_data/tumorDB.csv"]
		date_cols = ["TrappingDate", "TrapDate"]
		chip_list = self.sample_df["Microchip"]
		counts = []
		for i in range(2):
			df = pd.read_csv(default_paths[i])
			df = df[df["Microchip"].isin(chip_list)]
			trap_col = "TrapDate" if "TrapDate" in df.columns else "TrappingDate"
			count = df.groupby("Microchip")[trap_col].nunique()
			counts.append(count)
			chip_list = self.sample_df.loc[~self.sample_df["Microchip"].isin(count.index), "Microchip"]
		counts = pd.concat(counts)
		counts.name = "trap_count"
		counts = counts.astype("Int64")
		self.sample_df =  self.sample_df.merge(counts, on="Microchip", how="left")


	def undo_factor(self, pheno):
		factor_df = pd.DataFrame(self.factor_key)
		factor_pheno = factor_df[factor_df["pheno"] == pheno]
		sample_pheno = self.sample_df[pheno]
		undone_factors = sample_pheno.replace(factor_pheno["val"].to_list(), factor_pheno["key"].to_list())
		self.sample_df[pheno] = undone_factors
		self.factor_key = factor_df[factor_df["pheno"] != pheno].reset_index(drop=True).to_dict(orient='list')


	def undo_vcf_chip(self, df_path, chip_col="Microchip", header="infer", sep="\t"):
		df = pd.read_csv(df_path, sep=sep, header=header)
		df_six = df[chip_col].str[-6:]
		chips = self.sample_df["Microchip"]
		found = []
		for line in df_six:
			found.append(chips[chips.str[-6:] == line].iloc[0])
		df[chip_col] = found
		return df


	def unique(self, col):
		return self.sample_df[col].unique()


	def update_factors(self, pheno_list):
		for pheno in pheno_list:
			if pheno not in self.factor_key["pheno"]:
				continue
			self.undo_factor(pheno)
			self.to_factor(pheno)


	# === COVARIATE ===
	# Definition: 
	# 			  Add a categorical covariate to WPP samples of whether a devil was born before any DFTD arrived at WPP,
	# 			  during the period where the tetraploid variant was present, or after the normal diploid DFTD arrived.
	# 			  While this method can be used on its own, it is better to use with DFTD_YOB() at default settings.
	# 			  Dates are as follows:
	# 			  pre-DFTD:   YOB < 2007
	# 			  tetraploid: 2007 > YOB < 2011
	# 			  post-DFTD:  YOB > 2011
	# Arguments:
	# 			  None
	# Returns:
	# 			  None. Adds a "born_arrival" column to sample_df if it is not already present
	# Steps:
	# 			  1) Obtain just WPP samples
	# 			  2) Get the indices for pre-DFTD samples, tetraploid samples, and post-DFTD samples
	# 			  3) If the "born_arrival" column is not present, create it with default values of N/A
	# 			  4) Use the previous indices to appropriately add the DFTD label
	def WPP_tetraploid(self):
		WPP_subset = self.subset("Site", "WPP", inplace=False, verbose=False)
		pre_dftd = WPP_subset[WPP_subset["YOB"] < pd.to_datetime(2007, format="%Y")].index
		tetraploid = WPP_subset[(WPP_subset["YOB"] > pd.to_datetime(2007, format="%Y")) & (WPP_subset["YOB"] < pd.to_datetime(2011, format="%Y"))].index
		post_dftd = WPP_subset[WPP_subset["YOB"] > pd.to_datetime(2011, format="%Y")].index
		if not "born_arrival" in self.sample_df.columns:
			print("HYSOH")
			self.sample_df["born_arrival"] = "N/A"
		self.sample_df.loc[pre_dftd, "born_arrival"] = "pre-DFTD"
		self.sample_df.loc[tetraploid, "born_arrival"] = "tetraploid"
		self.sample_df.loc[post_dftd, "born_arrival"] = "post-DFTD"


	def write_pairs(self, outpath):
		pair_df = self.sample_df[["Microchip", "Tissue", "TumourNumber"]].copy()
		pair_df = self.__handle_triplets(pair_df)
		pairs = pair_df.groupby("Microchip")["Microchip"].count()
		pair_loc = pairs[pairs > 1].index
		pair_df = pair_df[pair_df["Microchip"].isin(pair_loc)]
		tumors = pair_df[pair_df["Tissue"] == "Tumour"].sort_values(by="Microchip").reset_index(drop=True)
		tumors = self.to_vcf_chip(tumors)
		hosts = pair_df[pair_df["Tissue"] == "Host"].sort_values(by="Microchip").reset_index(drop=True)
		hosts = self.to_vcf_chip(hosts)
		tumor_host = pd.concat([hosts, tumors], axis=1)
		tumor_host.columns = ["Hosts", "Tumors"]
		tumor_host.to_csv(outpath, index=False, sep="\t")


	def __assume_YOB(self, df, breeding_year="4/1"):
		YOB = df["YOB"].dropna()
		year_indices = YOB[~YOB.str.contains("/")].index
		df.loc[year_indices, "YOB"] = breeding_year + "/" + df.loc[year_indices, "YOB"]
		return df


	def __handle_triplets(self):
		triplets = self.get_pairs(3)
		if triplets.empty:
			return None
		print(f"Found {len(triplets)} tissue duplicates, please select only 1 of the duplicates to retain:")
		remove_loc = []
		for index in triplets.index:
			current_triplet = self.sample_df[self.sample_df["Microchip"] == index]
			if len(current_triplet[current_triplet["Tissue"] == "Tumour"]) > 1:
				dup = current_triplet[current_triplet["Tissue"] == "Tumour"]
			else:
				dup = dup = current_triplet[current_triplet["Tissue"] == "Host"]
			dup_i = list(dup.index.values)
			print(f'{dup[["Microchip", "Library number", "AnimalName", "Tissue", "YOB", "TrappingDate", "Site", "TumourID", "TumourNumber"]]}')
			retain_loc = input(f"Please select the index to be retained ({dup_i[0]} or {dup_i[1]}): ")
			dup_i.remove(int(retain_loc))
			remove_loc += dup_i
		self.sample_df = self.sample_df.loc[~self.sample_df.index.isin(remove_loc)]


	def __isnum(self, x):
		try:
			int(x)
		except:
			return False
		return True


	def __L_ID(self, paths):
		L_IDs = []
		for path in paths:
			if path[-4:] == ".csv":
				fnames = pd.read_csv(path)
				fnames = fnames["original"].to_list()
			else:
				fnames = os.listdir(path)
			L_IDs += list(set([file[:5] for file in fnames]))
		return L_IDs





class TumorSamples(Samples):
	def __init__(self,
		sample_csv_path=f"{Samples.base}/data/pheno_data/master_corrections.csv",
		tumor_csv_path=f"{Samples.base}/data/pheno_data/tumorDB.csv",
		id_paths=[],
		drop_DFTG=5,
		tissue="Host",
		extract_on="Host",
		full_tumor_df=False
	):
		super().__init__(sample_csv_path, id_paths)
		self.drop_DFTG = drop_DFTG
		self.full_tumor_df = full_tumor_df
		dates = ["TrapDate"]
		int_cols = ["TumourNumCertainty", "TumourMerged", "TumourDFTG", "TumourUlcerated", "TumourSecondaryInfection", "TumourBiopsy"]
		str_cols = ["Microchip", "TumourID", "TumourNumber"]
		int_cols = {col: "Int64" for col in int_cols}
		str_cols = {col: str for col in str_cols}
		type_dict = {**str_cols, **int_cols}
		tumor_df_full = pd.read_csv(tumor_csv_path, dtype=type_dict, parse_dates=dates, na_values="n")
		if tissue != "both":
			self.sample_df = self.sample_df[self.sample_df["Tissue"] == tissue]
			if self.sample_df.size == 0:
				print(f"!WARNING! sample_df has no entries! Tissue subsetted to \"{tissue}\", ensure spelling and capitalization are correct")
		self.failed_chips = self.sample_df[~self.sample_df["Microchip"].isin(tumor_df_full["Microchip"])]
		microchips = self.sample_df["Microchip"].unique()
		if extract_on != "both" and tissue == "both":
			microchips = self.sample_df[self.sample_df["Tissue"] == extract_on]["Microchip"].unique()
		if not full_tumor_df:
			self.tumor_df = tumor_df_full[tumor_df_full["Microchip"].isin(microchips)]
		else:
			self.tumor_df = tumor_df_full
		if drop_DFTG and not full_tumor_df:
			self.tumor_df, self.removed_samples, self.removed_tumors = self.__drop_DFTG(self.tumor_df, drop_DFTG)
		self.__init_vars()


	# =============================================================================================================
	# ======================================== tumorDB.csv Updating  START ========================================
	def add_tumor_samples(self, df, dup_cols, DFTG_filter=5, sample_df_only=True, date_diff_tolerance=3, print_found=True):
		df = df.copy()
		tumor_cols = self.tumor_df.columns.values.tolist()
		df_cols = df.columns.values.tolist()
		col_check = [True if col in tumor_cols else False for col in df_cols]
		if not all(col_check):
			col_check = np.array(col_check)
			df_cols = np.array(df_cols)
			print(f"!FAILED! The following columns are not in the Tumour Table: {', '.join(df_cols[~col_check])}")
			print("Please ensure that columns in the incoming DF are present in the Tumour Table")
			return None
		type_dict = dict(self.tumor_df[df_cols].dtypes)
		df = df.astype(type_dict)
		df["Microchip"] = df["Microchip"].astype(str)
		if DFTG_filter:
			df = df[df["TumourDFTG"] >= DFTG_filter]
		df = df.dropna(subset=dup_cols)
		df = df[~df.duplicated(subset=dup_cols)]
		dup_df = self.tumor_df[dup_cols].copy()
		dup_df = dup_df.dropna()
		dup_df["tumorDB"] = True

		merge_dups = pd.concat([dup_df, df])
		marked_dups = merge_dups.duplicated(subset=dup_cols, keep=False)
		final_df = merge_dups[~marked_dups]
		final_df = final_df[final_df["tumorDB"] != True]
		final_df = final_df.drop(columns="tumorDB")
		print(f"Found {final_df.shape[0]} new rows")
		if sample_df_only:
			final_df = final_df[final_df["Microchip"].isin(self.sample_df["Microchip"].unique())]
			print(f"Found {final_df.shape[0]} rows after subsetting on sample_df rows only")
		if date_diff_tolerance:
			maintained_rows = []
			for i in range(final_df.shape[0]):
				current_row = final_df.iloc[i]
				tumor_samples = self.tumor_df[self.tumor_df["Microchip"] == current_row["Microchip"]]
				min_date_diff = abs((tumor_samples["TrapDate"] - current_row["TrapDate"])).min().days
				if min_date_diff >= date_diff_tolerance or tumor_samples.empty:
					maintained_rows.append(i)
			final_df = final_df.iloc[maintained_rows]
			print(f"Found {final_df.shape[0]} rows differing by at least {date_diff_tolerance} days")
			merged_final = pd.concat([self.tumor_df, final_df]).sort_values(by=["Microchip", "TrapDate"])
		if print_found:
			print(final_df.to_string(index=False))
			print(f"Total new rows found: {final_df.shape[0]}")
			print(f"Total rows in tumor_df: {merged_final.shape[0]}")
		return merged_final


	# Definition: 
	# 			  Update NA TumourNumber and TumourID samples in tumorDB.csv by matching these samples with samples from my spreadsheet of sequenced samples.
	# 			  Optionally, any adjacent TrapDates with NA TumourNumbers can be updated. An adjacent sample shares a Microchip with the sample to be updated
	# 			  but differs in TrapDate. For adjacent samples to be updated, three conditions must be met:
	# 			  1) The sample must have multiple trappings (i.e., adjacent samples)
	# 			  2) All adjacent samples must also have NA TumourNumbers
	# 			  3) All adjacent samples must have only a single tumor
	# 			  Updating adjacent samples is based on the assumption that, if a devil was trapped multiple times and only has a single tumor each trapping,
	# 			  then that tumor is the same tumor for all such trappings.
	# Arguments:
	# 			  date_tolerance:  when matching tumorDB samples with my sequenced samples, this is the maximum number of days allowed to differ between the
	# 							   two before I consider the samples to not match
	# 			  update_adjacent: whether or not adjacent TrapDates for a sample should be updated as well (see Definition)
	# 			  outfolder: 	   the output folder to write all files to
	# Returns:
	# 			  None.
	# Steps:
	# 			  1) Remove any of my trailing letters from TumourNumbers
	# 			  2) Obtain the samples which are missing a Microchip or TumourNumber in the tumorDB and loop through each of these samples
	# 			  === LOOP START ===
	# 			  3) If one of the sequenced samples is missing a TrappingDate, save this sample (and reason why it failed) and continue
	# 			  3) Find the corresponding tumorDB sample by Microchip and the TrapDate with the min difference from my sample
	# 			  4) If this min difference is greater than the tolerance, record the sample and failure reason and continue
	# 			  5) Obtain the index of the min date diff then get all samples possessing this date
	# 			  6) Replace the NA tumorDB TumourNumber with my sample's TumourNumber (do the same with TumourID if my sample has a TumourID)
	# 			  7) If update_adjacent=True, ensure that the sample has multiple TrapDates, all TrapDate TumourNumbers are NA, and there is a single tumor
	# 				 per TrapDate. If the sample fails at any point, record the reason of adjacent update failure and continue. Otherwise, update all adjacent
	# 				 samples to have the same TumourNumber
	# 			  === LOOP END ===
	# 			  8) Grab the tumorDB samples which were updated and add a column explaining if adjacent samples were updated
	# 			  9) Grab my samples which failed to be updated and add a column explaining the reason of failure
	# 			  10) Write the updated tumorDB and two explanatory files to outfolder
	# GOTCHAS:
	# 			  - Updating adjacent samples is based on an assumption which is not necessarily true for all cases. However, if it is not done then
	# 				this assumes the single TrapDate for the tumor is the legitimate first observed date, which is also a potentially incorrect assumption
	# 			  - This method should be used in relative isolation from anything else. For example, this method removes my trailing identifiers and may thus
	# 				break any GEMMA analyses by omitting them when writing to_pheno()
	def update_TumourNumber(self, date_tolerance=3, update_adjacent=True, outfolder="../../data/pheno_data"):
		if not self.full_tumor_df:
			print("Please instantiate a TumorSamples object with full_tumor_df=True and run this again")
			exit()
		self.sample_df["TumourNumber"] = self.sample_df["TumourNumber"].str.replace(r"[a-zA-Z]$", "", regex=True)
		_, missing_TumourNumber, _ = self.get_sample_tumors(inplace=False)
		missing_date = 0
		updated = {}
		failed = {}
		for i in range(missing_TumourNumber.shape[0]):
			current = missing_TumourNumber.iloc[i]
			current_index = current.name
			if pd.isnull(current["TrappingDate"]):
				failed[current_index] = "Missing TrappingDate"
				continue
			tumor_samples = self.tumor_df[self.tumor_df["Microchip"] == current["Microchip"]]
			date_diff = abs(current["TrappingDate"] - tumor_samples["TrapDate"]).dt.days
			if date_diff.min() > date_tolerance:
				failed[current_index] = "TrappingDate difference too great"
				continue
			matching_date = tumor_samples.loc[date_diff.idxmin(), "TrapDate"]
			matching_tumor = tumor_samples[tumor_samples["TrapDate"] == matching_date]
			if matching_tumor.shape[0] > 1:
				failed[current_index] = "Multiple tumors at date"
				continue
			self.tumor_df.loc[matching_tumor.index, "TumourNumber"] = current["TumourNumber"]
			if current["TumourID"]:
				self.tumor_df.loc[matching_tumor.index, "TumourID"] = current["TumourID"]
			matching_index = matching_tumor.index[0]
			updated[matching_index] = "update_adjacent=False"
			if update_adjacent:
				if tumor_samples.shape[0] == 1:
					updated[matching_index] = "None. Sample only has a single TrapDate"
					continue
				if not tumor_samples["TumourNumber"].isna().all():
					updated[matching_index] = "None. Some non-NA values found"
					continue
				counts = tumor_samples.groupby("TrapDate")["Microchip"].count()
				if (counts > 1).any():
					updated[matching_index] = "None. Multiple tumors on at least one TrapDate found"
					continue
				self.tumor_df.loc[tumor_samples.index, "TumourNumber"] = current["TumourNumber"]
				updated[matching_index] = f"{tumor_samples.shape[0] - 1} adjacent samples updated"
		updated_samples = self.tumor_df.loc[list(updated.keys()), ["Microchip", "TrapDate", "TumourID", "TumourNumber"]]
		updated_samples["adjacent"] = updated.values()
		failed_samples = self.sample_df.loc[list(failed.keys()), ["Microchip", "TrappingDate", "TumourID", "TumourNumber"]]
		failed_samples["fail_reason"] = failed.values()
		print(f"{updated_samples.shape[0]} successfully updated samples | {failed_samples.shape[0]} failed samples")
		print(f"Writing the following files to {outfolder}")
		print("tumorDB_TumourNumber_update.csv [updated tumorDB]")
		print("TumourNumber_succeeded.csv [successfully updated samples]")
		print("TumourNumber_failed.csv [samples which failed to be updated]")
		self.tumor_df.to_csv(f"{outfolder}/tumorDB_TumourNumber_update.csv", index=False)
		updated_samples.to_csv(f"{outfolder}/TumourNumber_succeeded.csv", index=False)
		failed_samples.to_csv(f"{outfolder}/TumourNumber_failed.csv", index=False)
	# ========================================= tumorDB.csv Updating END ==========================================
	# =============================================================================================================



	# =============================================================================================================
	# ============================================== Phenotype START ==============================================
	# Definition: 
	# 			  Determine cases and controls based on whether a host was infected with DFTD. Because tumorDB.csv represents data collated from all
	# 			  available datasheets, it is assumed that a host not represented in tumorDB was never infected. This method can also utilize samples
	# 			  which were infected but only after the specific trap date (e.g., devil was trapped without DFTD but was infected sometime afterwards).
	# 			  See Margres et al. (2018; DOI: 10.1111/mec.14853) for a similar analysis.
	# Arguments:
	# 			  permit_infected: permit controls to be devils which were uninfected on the sampling capture (i.e., when biopsied and sequenced) but were infected at a later capture
	# 			  age_tolerance:   the required age of the devil in days to be considered a control
	# 			  stats:		   print out stats regarding the final cases and controls
	# Returns:
	# 			  None. Adds a "case" column to the samples DF where 1=case and 0=control
	# Steps:
	# 			  1) Obtain the difference in days from first to last trapping event
	# 			  2) Filter out devils with a single trap event or with fewer days between two events which fail the day_cutoff (these are also saved)
	# 			  3) Utilize the survival proxy or back calculated estimate to find the initial tumor data (see __init_tumor_date() for more details)
	# 			  4) Obtain the date for the most recent trapping event
	# 			  5) Subtract the tumor initial date from the max trap date to obtain devil survival days
	# GOTCHAS:
	# 			  - For duplicate hosts, only the first sample is maintained
	# 			  - Any samples lacking DFTD infection but failing the age_tolerance get NaN values (these are neither cases nor controls)
	# 			  - permit_infected=True indicates these are treated as potential controls such that those failing age_tolerance become NaN
	# 			  - permit_infected=False indicates that these samples will all be considered cases (I believe this is what was done in Margres. et al)
	# 			  - Most of my samples are paired host-pathos, thus this phenotype has very low power with 37 controls (assuming my script is working)
	# TODO:
	# 			  - Ensure permit_infected=False sets all samples infected after my trapDate to cases
	def case_control(self, permit_infected=True, age_tolerance=800, stats=False, chip_tables=[]):
		samples = self.subset("Host-tumour pair", False, inplace=False)
		tumors = self.tumor_df[self.tumor_df["Microchip"].isin(samples["Microchip"])]
		unextracted = self.subset_unextracted(inplace=False)
		unextracted = unextracted[unextracted["Host-tumour pair"] == False]

		if permit_infected:
			scores = tumors.groupby("Microchip")["TumourDFTG"].max()
			scores = tumors[tumors["Microchip"].isin(scores.index)]
			tumor_dates = scores.groupby("Microchip")["TrapDate"].min()
			tumor_dates = tumor_dates.to_frame().merge(self.sample_df[["Microchip", "TrappingDate"]], on="Microchip", how="left")
			tumor_dates.columns = ["Microchip", "tumor_trap", "sample_trap"]
			control_mchips = tumor_dates[tumor_dates["sample_trap"] < tumor_dates["tumor_trap"]]["Microchip"]
			print(control_mchips.to_string(index=False))
			control_mchips = pd.concat([control_mchips, unextracted["Microchip"]])
			controls = self.sample_df[self.sample_df["Microchip"].isin(control_mchips)]
		else:
			controls = unextracted
		controls = controls.drop_duplicates(subset=["Microchip"])
		controls["trap_age"] = (controls["TrappingDate"] - controls["YOB"]).dt.days
		failed_tolerance =  controls[controls["trap_age"] < age_tolerance]
		controls = controls[controls["trap_age"] >= age_tolerance]
		self.sample_df.loc[controls.index, "cases"] = "0"
		self.sample_df.loc[~self.sample_df.index.isin(controls.index), "cases"] = "1"
		self.sample_df.loc[failed_tolerance.index, "cases"] = float("NaN")

		if stats:
			cases = self.sample_df[self.sample_df['cases'] == '1']["Microchip"]
			controls = self.sample_df[self.sample_df['cases'] == '0']["Microchip"]
			cases_dup = len(cases) - len(cases.unique())
			controls_dup = len(controls) - len(controls.unique())
			print("=== CASE-CONTROL STATS ===")
			print(f"Cases: {len(cases.unique())}")
			print(f"Controls: {len(controls.unique())}")
			print(f"Total: {len(cases.unique()) + len(controls.unique())}")
			if cases_dup > 0:
				print(f"Duplicate cases: {cases_dup}")
			if controls_dup > 0:
				print(f"Duplicate controls: {controls_dup}")
			if chip_tables:
				found_dict = {}
				controls = self.sample_df[self.sample_df["cases"] == "0"]["Microchip"]
				for table in chip_tables:
					if controls.empty:
						break
					df = pd.read_csv(table)
					if not "Microchip" in df.columns:
						print(f"Table {table} does not contain a \"Microchip\" column. Skipping this table")
						continue
					found = controls[controls.isin(df["Microchip"])]
					if not found.empty:
						found_dict[table] = found.tolist()
						controls = controls[~controls.isin(found)]
				if found_dict:
					total = sum([len(v) for v in found_dict.values()])
					print(f"Controls found in other tables: {total}")
					for k,v in found_dict.items():
						print(f"\t{os.path.basename(k)}: {len(v)}")
				if not controls.empty:
					print(f"Controls that could not be found: {len(controls)}")
			print("==========================")


	# ANALYSIS PHENOTYPE
	def estimate_age(self, inplace=True, verbose_merge=False, proxy=True, vol_mode="max"):
		tmp_tumor_df = self.tumor_df[["Microchip", "TrapDate", "TumourDepth", "TumourLength", "TumourWidth"]].copy()
		dates_df = self.__init_tumor_date(tmp_tumor_df, proxy, vol_mode=vol_mode)
		YOB = self.sample_df.set_index("Microchip")["YOB"]
		YOB = YOB[~YOB.index.duplicated(keep="first")]
		merged = pd.concat([dates_df, YOB[dates_df.index]], axis=1)
		infection_age = (merged["init_tumor_date"] - merged["YOB"]).dt.days
		infection_age.name = "infection_age"
		infection_age = infection_age.reset_index()
		if verbose_merge:
			merged = merged.drop(columns=["YOB"])
			infection_age = merged.merge(infection_age, how="left", on="Microchip")
		if not inplace:
			return infection_age
		self.sample_df = self.sample_df.merge(infection_age, how="left", on="Microchip")


	# Definition: 
	# 			  Obtain a proxy for devil survival similar to Margres et al. (2018; DOI: 10.1111/mec.14853) using a tumor growth back calculation
	# 			  derived from Wells et al. (2017; DOI: 10.1111/ele.12776). Only devils with 2 trapping events and day_cutoff days between those
	# 			  trapping events are included. After the volume of the largest tumor (from tumors recorded on the earliest trap date) is obtained,
	# 			  the back calculation is used to find the number of days since the tumor was 3 mm^3. This date is then subtracted from the last
	# 			  date the devil was trapped to get survival in days.
	# Arguments:
	# 			  day_cutoff: 	 the minimum number of days a devil must be trapped between two trapping events for inclusion in the analysis <int>
	# 			  proxy:		 True if the proxy (last_trap - first_trap) should be used over the back-calc survival estimate
	# 			  verbose_merge: whether the intermediate calculation values should be added to the samples DF. Mostly for debugging <boolean>
	# Returns:
	# 			  None. Adds at least a "devil_survival_days" column to the samples DF
	# Steps:
	# 			  1) Obtain the difference in days from first to last trapping event
	# 			  2) Filter out devils with a single trap event or with fewer days between two events which fail the day_cutoff (these are also saved)
	# 			  3) Utilize the survival proxy or back calculated estimate to find the initial tumor data (see __init_tumor_date() for more details)
	# 			  4) Obtain the date for the most recent trapping event
	# 			  5) Subtract the tumor initial date from the max trap date to obtain devil survival days
	# GOTCHAS:
	# 			  None
	# TODO:
	# 			  - The calculated dates often do not match up with Margres et al. (2018) table S1.
	# 			  - Attempt to change the calculation using the sum of all tumor volumes on the min trap date.
	# 			  - Use the Compare class to test calculation similarities
	def survival_proxy(self, day_cutoff=40, proxy=True, verbose_merge=False, handle_date="smaller", vol_mode="max"):
		tmp_tumor_df = self.tumor_df[["Microchip", "TrapDate", "TumourDepth", "TumourLength", "TumourWidth"]].copy()
		tmp_tumor_df["TrapDate"] = pd.to_datetime(tmp_tumor_df["TrapDate"])
		
		date_diff = tmp_tumor_df.groupby("Microchip")["TrapDate"].transform(lambda x: x.max() - x.min())
		tmp_tumor_df["date_diff"] = date_diff.dt.days
		self.failed_date_delta = pd.Series(tmp_tumor_df[tmp_tumor_df["date_diff"] < day_cutoff]["Microchip"].unique())
		tmp_tumor_df = tmp_tumor_df[tmp_tumor_df["date_diff"] >= day_cutoff]

		dates_df = self.__init_tumor_date(tmp_tumor_df, proxy, handle_date, vol_mode=vol_mode)
		last_trap = tmp_tumor_df.groupby("Microchip")["TrapDate"].max().loc[dates_df.index]
		devil_survival_days = last_trap - dates_df["init_tumor_date"]
		
		calc_df = pd.DataFrame({
			"Microchip": dates_df.index.values,
			"last_trap": last_trap.dt.date.values,
			"init_tumor_date": dates_df["init_tumor_date"].dt.date.values,
			"devil_survival_days": devil_survival_days.dt.days.values})
		if not proxy:
			calc_df.insert(loc=2, column="volume (cm^3)", value=np.round(dates_df["tumor_volume"].values, 2))
			calc_df.insert(loc=3, column="back_calc", value=dates_df["back_calc"].values)
			calc_df.insert(loc=4, column="first_trap", value=dates_df["min_date"].dt.date.values)
		if not verbose_merge:
			calc_df = calc_df[["Microchip", "devil_survival_days"]]
		self.sample_df = self.sample_df.merge(calc_df, how="left", on="Microchip")


	# Definition: 
	# 			  A modified version of estimate_age() which removes the age effect. This basically looks at how long it took a devil
	# 			  to contract DFTD after it was able to contract the disease. This correction is done by subtracting one of two things:
	# 			  the age of the devil when DFTD arrived OR the min infection_age. The age of the devil on DFTD arrival is subtracted
	# 			  because it is impossible for a devil to get DFTD before the tumor arrives at the devil's site. The min infection_age
	# 			  subtraction is done becasue it is assumed that, based on devil biting behavior, devils younger than this age will not
	# 			  engage in social biting and thus cannot contract DFTD. The value subtracted is the higher of the two, as in: assume a
	# 			  min infection_age of 365d. A devil 800d old on DFTD arrival who gets infected at 850d old has a transmission_age of
	# 			  50d, because it was impossible for the devil to contract DFTD for the first 800d of its life. However, a devil that is
	# 			  200d old on DFTD arrival that contracts DFTD at 375d old has a transmission_age of 10d. Although the devil could not
	# 			  contract DFTD earlier than 200d (due to DFTD not yet existing at its site), the devil was not able to contract DFTD
	# 			  before 365d old due to behavior, and thus the larger min infection_age is subtracted.
	# Arguments:
	# 			  verbose_merge: add all columns used in the transmission_age calculation to sample_df (primarily for debugging)
	# 			  proxy:		 if the initial infection age should be a proxy or back-calc estimate (should only be done for WPP)
	# Returns:
	# 			  None. Adds a "transmission_age" column to sample_df
	# Steps:
	# 			  1) Obtain the infection_age and arrival_date if it was not already added to sample_df
	# 			  2) Calculate arrival - YOB, where positive values indicate a devil born X days prior to DFTD arrival
	# 			  3) For any arrival - YOB less than min infection_age, set the value to 0
	# 			  4) Make a min infection_age Series, setting any value to 0 which is not 0 in the arrival - YOB
	# 				 NOTE: This and 3) ensure that only YOB - arrival OR min infection_age (whichever is larger) is subtracted from infection_age
	# 			  5) Subtract YOB - arrival OR min infection_age from infection age (this removes the non-viable devil infection period)
	# 			  6) Restore the original sample_df or add all calculation columns if verbose_merge
	# 			  7) Add the transmission_age column
	# 			  8) Remove any row with transmission_age <= 0 (including NaN). Some devils were infected before DFTD arrived at their site,
	# 				 possibly due to migration, the arrival_date being an estimate, or recorder error. These samples are removed
	# Gotchas:
	# 			  - DFTD arrival for WPP is assumed to be 2011
	# 			  - Devil 982009104872893 has a negative transmission_age (check for more when I have >312 samples extracted)
	# 			  - This method resets the index of the sample_df (this shouldn't interfere with any other methods)
	# 			  - The correction ALWAYS generates at least a single 0 value. When log transforming, this will be -inf. For now, manually correct
	# 			  	this to a sensible value smaller than all other values
	# TODO:
	# 			  - Make this work with tumors, including determining DFTD arrival for tumor sites and determining if our sequenced tumor is the initial infection
	def transmission(self, verbose_merge=False, proxy=True, vol_mode="max"):
		self.sample_df = self.sample_df.reset_index(drop=True)
		original_df = self.sample_df.copy()
		if not "infection_age" in self.sample_df.columns:
			self.estimate_age(proxy=proxy, vol_mode=vol_mode)
		if not "arrival_date" in self.sample_df.columns:
			self.add_arrival()
		min_age = self.sample_df["infection_age"].min()
		born_before_arrival = (self.sample_df["arrival_date"] - self.sample_df["YOB"]).dt.days.copy()
		born_before_arrival[born_before_arrival < min_age] = 0
		min_series = pd.Series([min_age] * len(born_before_arrival))
		min_series.loc[born_before_arrival[born_before_arrival != 0].index] = 0
		transmission_age = self.sample_df["infection_age"] - born_before_arrival - min_series
		if verbose_merge:
			self.sample_df["born_before_arrival"] = born_before_arrival
			self.sample_df["min_infection_age"] = min_series
		else:
			self.sample_df = original_df
		colname = "transmission_age" if proxy else "transmission_age_est"
		self.sample_df[colname] = transmission_age
		init_N = self.sample_df.shape[0]
		self.negative_transmission = self.sample_df[self.sample_df[colname] < 0]
		self.sample_df = self.sample_df[self.sample_df[colname] >= 0]
		n = self.sample_df.shape[0]
		removed = init_N - n
		negative_N = self.negative_transmission.shape[0]
		print(f"*transmission()* Removing {removed} negative age ({negative_N}) or NaN ({removed - negative_N}) samples out of {init_N} ({n} remaining)")


	# Definition: 
	# 			  Get the tumor growth rate for each devil. Tumor growth is logistic, but most devils have one or two trapping events,
	# 			  making it difficult to estimate the logistic growth rate for an individual. Plus, I don't believe (or know of) a valid
	# 			  way of estimating logistic growth rate with just a few data points; the ideal scenario would be that each devil has
	# 			  tens or hundreds of trappings and I could fit a logistic model to each individual devil. Because this is not the case,
	# 			  I currently estimate growth rate as a linear function, fitting a linear model to each individual devil. For devils with
	# 			  two trappings (the most common scenario), this is equivalent to a basic rise/run slope. For devils with >2 trappings,
	# 			  an ordinary least squares linear model is fit and the date coefficient is used as growth rate. Currently, tumor volume
	# 			  is calculated as the sum of all tumor volumes at a given TrapDate, a measure called tumor load (see Wells et al., 2017;
	# 			  https://doi.org/10.1111/ele.12776). Growth rates are in units of cm^3/day.
	# Arguments:
	# 			  calc_df: return the dataframe containing all intermediate values used in calculating growth rate
	# Returns:
	# 			  Adds a "growth_rate" column to sample_df. Returns the calc_df if calc_df=True
	# Steps:
	# 			  1) Obtain the relevant tumor_df columns and drop any rows with a missing tumor length, width, or depth
	# 			  2) Extract only samples with more than a single trapping event
	# 			  3) Get the volume for each individual tumor and create a new col: "vol"
	# 			  4) Convert the TrapDates to days from first TrapDate (i.e., first TrapDate always becomes 0)
	# 			  5) Sum the volumes per TrapDate to obtain tumor load, convert the series to a DF, and divide by 1000 to convert vol to cm^3
	# 			  6) For each Microchip and TrapDate, get the first days from first TrapDate
	# 			  	 NOTE: this is necessary because multiple tumors on a single TrapDate create duplicate days from first TrapDate
	# 			  7) Iterate through each unique mchip and fit a sample-specific linear regression as volume = date*x + b, saving the date intercept to a list
	# 			  8) Merge the rates with sample_df. If calc_df=True, merge the rates with the vols DF and return this
	# Gotchas:
	# 			  - Some samples with >2 trappings may drop out if they're missing a volume measurement
	# 			  - This method assumes tumor measurements are all in mm
	# 			  - The simplifying linear assumption may be a significant confounder
	# TODO:
	# 			  - Add a parameter to calculate per-tumor growth rate (rather than a summed tumor load)
	# 			  - Make the method amenable to a paired tumor-host analysis
	# 			  - Check the linearity assumption (perhaps relate the growth rates to a starting tumor)
	# 				volume to see if there is a significant association between start_volume and rate
	# 			  - Make this logistic growth rather than linear
	def tumor_growth(self, calc_df=False):
		tmp_tumor_df = self.tumor_df[["Microchip", "TrapDate", "TumourDepth", "TumourLength", "TumourWidth"]].copy()
		tmp_tumor_df = tmp_tumor_df.dropna(subset=["TumourDepth", "TumourLength", "TumourWidth"])
		date_counts = tmp_tumor_df.groupby(["Microchip"])["TrapDate"].nunique()
		date_counts = date_counts[date_counts > 1]
		multi_traps = tmp_tumor_df.loc[tmp_tumor_df["Microchip"].isin(date_counts.index)]
		multi_traps = multi_traps.assign(vol=(multi_traps["TumourDepth"] * multi_traps["TumourLength"] * multi_traps["TumourWidth"]))
		multi_traps["date"] = multi_traps.groupby("Microchip")["TrapDate"].transform(lambda x: (x - x.min()).dt.days)
		vols = multi_traps.groupby(["Microchip", "TrapDate"])["vol"].sum().to_frame()
		vols /= 1000
		vols["date"] = multi_traps.groupby(["Microchip", "TrapDate"])["date"].first()
		chips = vols.index.get_level_values("Microchip").unique()
		rates = []
		for chip in chips:
			sample = vols.loc[chip]
			model = sm.ols(formula="vol ~ date", data=sample).fit()
			rates.append(model.params["date"])
		rate_df = pd.DataFrame({"Microchip": chips, "growth_rate": rates})
		self.sample_df = self.sample_df.merge(rate_df, on="Microchip", how="left")
		if calc_df:
			calc_df = vols.merge(rate_df, on="Microchip")
			return calc_df


	# Definition:
	# 			  This is a convenience function doing one of two things:
	# 			  1) Combine the functionality of __min_cap_volume() and __growth_back_calculation() to find the initial tumor date (back-calc estimate)
	# 			  2) Find the minimum trapping date for which a tumor was observed (survival proxy)
	# Arguments:
	# 			  tumor_df:    a DF with at least the following cols: Microchip, TrapDate, TumourDepth, TumourLength, TumourWidth <pd.DataFrame>
	# 			  proxy_flag:  True if using a basic orixy for the init date tumor (i.e., first date trapped with a tumor present).
	# 						   Otherwise, employ the back-calc to estimate the init tumor date
	# 			  handle_date: how to handle back calc dates which are later than first observed dates (i.e., proxy survival is longer than back calc). 
	# 						   These samples can be dropped (drop), maintain the back calc date (ignore), or use the min observed tumor date (smaller)
	# Returns:
	# 			  a DF with Microchip as the index and one of two column configurations:
	# 			  1) min_obs_date, tumor_volume, min_date, back_calc, init_tumor_date
	# 			  2) init_tumor_date
	# Steps:
	# 			  1) Run __min_cap_volume() to obtain the minimum observed tumor trap date (i.e., survival proxy) and drop the tumor_volume column (not needed for proxy)
	# 			  2) Run __min_cap_volume() again but this time do so for the back calc estimate (i.e., minimum date with tumor volume measurements)
	# 			  3) Divide the volume by 1000 to convert from mm^3 to cm^3 (this will have to be changed if tumor measurements aren't done in mm)
	# 			  BACK CALCULATION
	# 			  4) Run __growth_back_calculation() to obtain days since init tumor (negative days)
	# 			  5) Get the initial tumor date by subtracting the back calc from the min date
	# 			  6) Combine the back calc values and the proxy min date into a single DF
	# 			  7) Find any samples where the min observed date is earlier than the estimated tumor init date
	# 			  8) Save these samples to self.smaller_obs for further manual inspection
	# 			  9) Handle these odd samples by dropping them, ignoring them and utilizing the back calc, or utilizing the min observed date
	# 			  PROXY
	# 			  4) Rename "min_obs_date" (minimum trap date, which considers NA volume rows, obtained from __min_cap_volume) to "init_tumor_date"
	def __init_tumor_date(self, tumor_df, proxy_flag, handle_date="smaller", vol_mode="max"):
		handle_date_list = ["drop", "ignore", "smaller"]
		proxy = self.__min_cap_volume(tumor_df, False, mode=vol_mode)
		proxy = proxy.drop(columns=["tumor_volume"])
		proxy.columns = ["min_obs_date"]
		estimate = self.__min_cap_volume(tumor_df, True, mode=vol_mode)
		estimate["tumor_volume"] /= 1000
		if not proxy_flag:
			back_calc = self.__growth_back_calculation(estimate["tumor_volume"], is_cm=True)
			init_tumor_dates = estimate["min_date"] - pd.to_timedelta(-1 * back_calc, unit="d")
			init_tumor_dates.name = "init_tumor_date"
			back_calc = pd.Series(back_calc, name="back_calc", index=estimate.index)
			final_df = pd.concat([proxy, estimate, back_calc, init_tumor_dates], axis=1)
			final_df.index.rename("Microchip", inplace=True)
			final_df = final_df[~final_df["init_tumor_date"].isna()]
			smaller_obs = final_df[final_df["init_tumor_date"] > final_df["min_obs_date"]].copy()
			smaller_obs["init_date_diff"] = (smaller_obs["init_tumor_date"] - smaller_obs["min_obs_date"]).dt.days
			self.smaller_obs = smaller_obs.reset_index()
			self.smaller_obs = self.smaller_obs.rename(columns={"index": "Microchip"})
			print(f"*WARNING* found {smaller_obs.shape[0]} samples with an earlier observed tumor date than estimated date (try print_smaller_obs() for details)")
			print(f"          These samples will be handled using method={handle_date}")
			if handle_date == "drop":
				final_df = final_df.drop(smaller_obs.index)
			elif handle_date == "ignore":
				pass
			elif handle_date == "smaller":
				final_df.loc[smaller_obs.index, "init_tumor_date"] = smaller_obs["min_obs_date"]
			else:
				print(f"Argument \"{handle_date}\" is not valid for handle_date. Please use: {', '.join(handle_date_list)}")
				return None
		else:
			final_df = proxy.rename(columns={"min_obs_date": "init_tumor_date"})
		return final_df


	# Definition: 
	# 			  A logistic growth back calculation derived from Wells et al. (2017; DOI: 10.1111/ele.12776) which was based on data from WPP.
	# 			  Output is a negative number in days representing days from the min capture date since the tumor was 3 mm^3.
	# 			  The model is not accurate enough to estimate beyond day resolution, and thus rounding to the nearest day 
	# 			  (default round_flag) is advised. This also stores any volumes greater than mmax.
	# Arguments:
	# 			  tumor_volumes: a vector/Series of tumor volumes <pd.Series>
	# 			  mmax: 		 max tumor size in cm^3 (from growth model) <int>
	# 			  init_size: 	 initial tumor size (from growth model) <int>
	# 			  alpha: 		 scale parameter of the logistic growth curve (from growth model) <int>
	# 			  mmax_cutoff:	 how much over mmax can be before discarding the sample (these over values are set to vol=200)
	# 			  is_cm: 		 flag if the incoming volume is in centimeters <boolean>
	# 			  round_flag: 	 flag if the back calculated volumes should be rounded to the nearest int <boolean>
	# Returns:
	# 			  numpy array containing negative numbers representing days since tumor volume = 3 mm^3
	# Steps:
	# 			  1) Convert mm to cm if is_cm=True
	# 			  2) Determine volumes >= mmax, setting those volumes to NaN and storing only those microchips and volumes not already in the over_mmax dict
	# 			  3) Perform the back calculation
	# 			  4) Round the calculations to the nearest int if round_flag=True and return
	# Gotchas:
	# 			  - Formula domain: (0, mmax)
	# 			  - Formula range: (-inf, 0); as volume -> mmax, days -> -inf
	def __growth_back_calculation(self, tumor_volumes, mmax=202, init_size=0.0003, alpha=0.03, mmax_cutoff=100, is_cm=False, round_flag=True):
		tumor_volumes = tumor_volumes.copy()
		if not is_cm:
			tumor_volumes /= 1000
		over = tumor_volumes[(tumor_volumes + 1) >= mmax]
		if over.size > 0:
			vstring = "volumes" if len(over) > 1 else "volume"
			print(f"*WARNING* found {len(over)} {vstring} greater than mmax while performing the back calculation (try print_over_mmax() for details)")
			acceptable_range = over[over <= mmax + mmax_cutoff]
			print(f"          {acceptable_range.size} of these are within the {mmax_cutoff} cutoff, setting these to volume=200 cm^3")
			tumor_volumes[((tumor_volumes + 1) >= mmax) & (tumor_volumes <= mmax + mmax_cutoff)] = 200
			tumor_volumes[(tumor_volumes + 1) >= mmax] = np.nan
		if self.over_mmax["Microchip"]:
			over = over.drop(self.over_mmax["Microchip"])
		if len(over > 0):
			self.over_mmax["Microchip"] += list(over.index.values)
			self.over_mmax["volume"] += list(over.values)
		back_calculation = np.log((mmax / (tumor_volumes.values + 1) - 1) / (mmax - init_size)) / alpha
		if round_flag:
			back_calculation = np.round(back_calculation)
		return back_calculation


	# Definition: 
	# 			  Obtain the tumor volume Series for the minimum capture date with microchips as the index.
	# 			  When a host has multiple tumors, different grouping functions can be used (e.g., max or sum). 
	# 			  If this is being used with the growth back calculation, max should be set for mode. Drops rows with an NA measurement.
	# Arguments:
	#			  tumor_df: a DF with at least the following cols: Microchip, TrapDate, TumourDepth, TumourLength, TumourWidth <pd.DataFrame>
	# 			  dropna:   drop any trap date with NA in depth, length, or width cols. This is necessary for the back calculation
	# 			  mode: 	method of selecting a single volume when multiple tumors are present at the same date <string>
	# Returns: 
	# 			  a DF containing Microchip as the index and cols tumor_volume (in starting units) and min_date (representing the minimum date for a sample)
	# Steps:
	#			  1) get min date for a sample (group on Microchip)
	# 			  2) Iterate through each of these microchips, find all rows with the min date, and generate a DF from these
	# 			  3) Obrain tumor volumes via length x width x depth
	# 			  4) Group these volumes on Microchip and select a single volume for each microchip (i.e., max, mean, sum)
	# Gotchas:
	# 			  - Conditionally drop rows if any NA is found in TumourDepth, TumourLength, or TumourWidth. If this is not done, then volumes with a missing
	# 			  	dimension are NaN. This is fine for the proxy, but will throw an error for the back calculation estimate.
	# TODO:
	# 			  - The loop is an inelegant and inefficient solution, this should be vectorized
	def __min_cap_volume(self, tumor_df, dropna, mode="max"):
		if dropna:
			tumor_df = tumor_df.dropna(subset=["TumourDepth", "TumourLength", "TumourWidth"])
		min_date_groups = tumor_df.groupby("Microchip")["TrapDate"].min()
		min_date_groups.name = "min_date"
		min_dates = {"Microchip": min_date_groups.index, "TrapDate": min_date_groups.values}
		first_tumor_list = []
		
		for i in range(len(min_dates["Microchip"])):
			chip = min_dates['Microchip'][i]
			date = min_dates['TrapDate'][i]
			chip_group = tumor_df[tumor_df["Microchip"] == chip]
			first_tumor_list.append(chip_group[chip_group["TrapDate"] == date])
		
		first_tumor_df = pd.concat(first_tumor_list)
		first_tumor_df["tumor_volume"] = first_tumor_df["TumourDepth"] * first_tumor_df["TumourLength"] * first_tumor_df["TumourWidth"]
		volume_cluster = first_tumor_df.groupby("Microchip")["tumor_volume"]
		if mode == "mean":
			volume_cluster = volume_cluster.mean()
		elif mode == "sum":
			volume_cluster = volume_cluster.sum()
		else:
			volume_cluster = volume_cluster.max()
		return pd.concat([volume_cluster, min_date_groups], axis=1)
	# =============================================== Phenotype END ===============================================
	# =============================================================================================================



	# ============================================== Covariate START ==============================================
	# =============================================================================================================
	def first_infection_tload(self, na_report="col"):
		tumor_df = self.tumor_df.dropna(subset=["TumourLength", "TumourWidth", "TumourDepth"], how="all").copy()
		tumor_df["NA_flag"] = False
		na_indices = tumor_df[(tumor_df["TumourLength"].isna()) | (tumor_df["TumourWidth"].isna()) | (tumor_df["TumourDepth"].isna())].index
		tumor_df.loc[na_indices, "NA_flag"] = True
		tumor_df.loc[na_indices, ["TumourLength", "TumourWidth", "TumourDepth"]] = tumor_df.loc[na_indices, ["TumourLength", "TumourWidth", "TumourDepth"]].fillna(value=1)
		tumor_df["vol"] = tumor_df["TumourLength"] * tumor_df["TumourWidth"] * tumor_df["TumourDepth"]
		tload = tumor_df.groupby(["Microchip", "TrapDate"])["vol"].sum().round(1)
		tload.name = "init_tload"
		tumor_df = tumor_df.merge(tload, on=["Microchip", "TrapDate"], how="left")
		min_date_tload = tumor_df.loc[tumor_df.groupby("Microchip")["TrapDate"].idxmin(), ["Microchip", "NA_flag", "init_tload"]]
		cols = ["Microchip", "init_tload"]
		if na_report == "str":
			min_date_tload["init_tload"] = min_date_tload["init_tload"].astype(str)
			min_date_tload.loc[min_date_tload["NA_flag"], "init_tload"] = ">=" + min_date_tload.loc[min_date_tload["NA_flag"], "init_tload"]
		elif na_report == "col":
			cols.append("NA_flag")
		min_date_tload = min_date_tload[cols]
		self.sample_df = self.sample_df.merge(min_date_tload, on="Microchip", how="left")


	# Definition: 
	# 			  Determine if the sequenced tumor is the first tumor to infect the devil. The largest tumor at the earliest TrapDate is considered
	# 			  to be the first tumor to infect the devil. If there is only a single tumor at the min TrapDate, it is automatically the first
	# 			  infected tumor. In the absence of any tumor measurement data, the method looks to see if there is a single tumor at the earliest
	# 			  TrapDate. If there is, it is considered the first infected. If there are multiple tumors, then the sample is considered a multi
	# 			  sample. Because the method compares TumourNumbers, a missing TumourNumber results in a value of NA. Of course, if the tumor
	# 			  cannot be found in tumorDB, or it is dropped due to TumourDFTG<5, it automatically gets NA. Output labels are as follows:
	# 			  yes = 	 first infected tumor
	# 			  no  = 	 not the first infected tumor
	# 			  NaN = 	 missing Microchip or TumourNumber; first infected tumor is indeterminate 
	# 			  multiple = missing tumor measurements and multiple tumors at min TrapDate; first infected tumor is indeterminate 
	# Arguments:
	#			  None.
	# Returns: 
	# 			  None. Adds a "first_infected" column
	# Steps:
	#			  1) Ensure get_sample_tumors() has not yet been called
	# 			  2) Drop any tumors with NA for all dimension cols
	# 			  === VOLUME METHOD ===
	# 			  3) For any tumors with at least one measurement, fill NAs with the value 1
	# 			  4) Obtain tumor volume in mm^3
	# 			  5) Obtain the largest volume for a given TrapDate by clustering on Microchip and TrapDate and using loc of the index
	# 			  6) Get the min TrapDate vie clustering on Microchip and using loc of the min index. The remaining tumors are the first infected
	# 				 based on their volumes
	# 			  === MIN DATE METHOD ===
	# 			  7) Find any samples which failed to get a first infected tumor (either due to not being present in tumorDB or lacking measurements)
	# 			  8) Count the number of tumors for a given TrapDate by clustering on Microchip and TrapDate. Merge with failed DF
	# 			  9) Get the loc of the min TrapDates and use the previously obtained counts to determine which of these min dates have 1 tummor
	# 			  10) Add these single tumor min dates to the volume first infected samples
	# 			  === END ===
	# 			  11) Merge the first infected tumor numbers with sample_df and store indices where first infected = NA
	# 			  12) Compare where my sampled TumourNumber (scrubbed for letter IDs) equals the first infected number
	# 			  13) Fill the NAs back in using the stored indices (converts the bool series to a float)
	# 			  14) Fill in the value "2.0" where any Microchip with >1 tumor and no measurements at the min TrapDate were found
	# 			  15) Replace the floats with appropriate strings
	# 			  16) Merge with sample_df
	# Gotchas:
	# 			  - Tumors with NA in length, width, and depth are dropped
	# 			  - Tumors with at least one measurement in length, width, or depth are maintained. The missing measurements are given a value
	# 				of 1
	# 			  - The largest tumor, even if it is only bigger by 1 mm^3, is considered the first infected
	# 			  - TrapDates other than the min are ignored. Thus, if the min TrapDate has a 1 mm^3 tumor and a TrapDate one week later
	# 				has a different tumor 200 cm^3, the 1 mm^3 tumor is the first infected (this scenario should be very rare)
	# 			  - None of my currently sequenced tumors have a label of multiple
	# TODO:
	# 			  - Test the "multiple" label by making up data and adding it to one of the "failed_single" min TrapDates
	def pre_infection(self):
		if self.sample_tumor_flag:
			print("This method cannot be used after running get_sample_tumors(). Please run this method before extracting sample tumors")
			exit()
		na = self.tumor_df[self.tumor_df[["TumourLength", "TumourWidth", "TumourDepth"]].isna().all(axis=1)]
		tumor_df = self.tumor_df.dropna(subset=["TumourLength", "TumourWidth", "TumourDepth"], how="all").copy()
		tumor_df[["TumourLength", "TumourWidth", "TumourDepth"]] = tumor_df[["TumourLength", "TumourWidth", "TumourDepth"]].fillna(value=1)
		tumor_df["vol"] = tumor_df["TumourLength"] * tumor_df["TumourWidth"] * tumor_df["TumourDepth"]
		max_tumor_i = tumor_df.groupby(["Microchip", "TrapDate"])["vol"].idxmax()
		tumor_df = tumor_df.loc[max_tumor_i]
		min_trap_i = tumor_df.groupby("Microchip")["TrapDate"].idxmin()
		tumor_df = tumor_df.loc[min_trap_i]
		failed = self.sample_df[~self.sample_df["Microchip"].isin(tumor_df["Microchip"])]
		failed_tumor_df = self.tumor_df[self.tumor_df["Microchip"].isin(failed["Microchip"])]
		failed_tumor_count = failed_tumor_df.groupby(["Microchip", "TrapDate"])["TumourNumber"].count()
		failed_tumor_count.name = "tumor_count"
		failed_tumor_df = failed_tumor_df.merge(failed_tumor_count, on="Microchip", how="left")
		failed_min_trap_i = failed_tumor_df.groupby("Microchip")["TrapDate"].idxmin()
		failed_tumor_df = failed_tumor_df.loc[failed_min_trap_i]
		failed_single = failed_tumor_df[failed_tumor_df["tumor_count"] == 1]
		failed_multi = failed_tumor_df.loc[failed_tumor_df["tumor_count"] > 1, "Microchip"]
		full_df = pd.concat([tumor_df, failed_single])
		full_df = full_df.rename(columns={"TumourNumber": "first_number"})
		merged = self.sample_df.merge(full_df[["Microchip", "first_number"]], on="Microchip", how="left")
		na_vals = merged[merged["first_number"].isna()].index
		first_infection = merged["TumourNumber"].str.replace(r"[a-zA-Z]$", "", regex=True) == merged["first_number"]
		first_infection[na_vals] = np.nan
		merged["first_infection"] = first_infection
		if not failed_multi.empty:
			merged.loc[merged["Microchip"].isin(failed_multi), "first_infection"] = 2.0
		merged["first_infection"] = merged["first_infection"].replace({0: "no", 1: "yes", 2: "multiple"})
		self.sample_df = self.sample_df.merge(merged[["Microchip", "TumourNumber", "first_infection"]], on=["Microchip", "TumourNumber"], how="left")


	# COVARIATE
	def tumor_count(self):
		num_tumors = self.tumor_df.groupby(["Microchip"])["TumourNumber"].unique().str.len()
		num_tumors = num_tumors.to_frame().reset_index()
		num_tumors = num_tumors.rename(columns={"TumourNumber": "tumor_count"})
		self.sample_df = self.sample_df.merge(num_tumors, how="left", on="Microchip")


	def coinfection(self):
		if self.sample_tumor_flag:
			print("This method cannot be used after running get_sample_tumors(). Please run this method before extracting sample tumors")
			exit()
		# print(self.tumor_df.groupby("Microchip")["TrapDate"].min())
		date_clusters = self.tumor_df.groupby(["Microchip", "TrapDate"])["TrapDate"].count()
		# print(date_clusters.groupby("Microchip").loc["TrapDate"].min())
		print(date_clusters.reset_index(0))
	# =============================================== Covariate END ===============================================
	# =============================================================================================================



	# ============================================== Subsetting START =============================================
	# =============================================================================================================
	# Definition: 
	# 			  Subset the tumor_df to include just the tumors sequenced in the sample_df. Ideally, every sample_df tumor sample will exactly match a single
	# 			  tumor_df tumor in 4 columns: Microchip, TumourNumber, TumourID, and TrappingDate. I consider a tumor's Microchip and TumourNumber to be
	# 			  its primary identifiers such that a failure to match either of these indicates the tumor is missing from tumor_df. Further validation,
	# 			  and importantly for phenotypes selection of the proper TrappingDate, is provided by matching TrappingDate and TumourID. However, in
	# 			  some cases a mismatch of one or both of these occurs. For example, it is possible that a tumor matches its TumourID but has a TrappingDate
	# 			  differing by hundreds of days. This method was designed to dynamically extract tumors considering these inconsistencies such that
	# 			  tumors may be extracted just by matching closest dates irrespective of TumourID mismatches, extract based on TumourID matches irrespective
	# 			  of TrappingDate mismatches, or extract based on the TrappingDate where TumourIDs match. Any tumors with missing TumourIDs (either in the
	# 			  tumor_df or sample_df) are extracted based solely upon TrappingDate matches.
	# Arguments:
	# 			  date_tolerance: The number of days TrappingDate mismatches can differ before dropping the sample. "None" means samples will not
	# 							  be dropped based on date differences whereas "0" means any sample with differing TrappingDates will be dropped
	# 			  id_match:		  Indicates any samples with matching TumourIDs will be retained irrespective of differences in TrappingDates
	# 			  closest_date:	  Indicates TrappingDate differences will be considered only for samples with the closest dates, irrespective of
	# 							  TumourID matches. Setting this to False will use the TrappingDate of the matching TumourID (if a match exists),
	# 							  even if this date is not closest to the sample_df TrappingDate
	# 			  specific_tumor: Indicates if *just* the single matching tumor should be extracted (this is likely less useful for most phenotypes)
	# 			  stats: 		  Print extraction stats
	# 			  inplace:		  Replace the tumor_df with the extracts (returns them otherwise)
	# Returns:
	# 			  Extracted tumors if inplace=False
	# Steps:
	# 			  1) Get just the samples in sample_df (and remove trailing A/B TumourNumber IDs) which have a tumor_df entry and intialize a bunch of lists
	# 			  2) Loop through every one of these sample_df samples
	# 			  === LOOP START ===
	# 			  3) Obtain the tumor_df entry (or entries) which matches the current sample_df entry's Microchip and TumourNumber
	# 			  	3a) If no tumor_df entry is found, save the sample_df index in not_found_indices and continue to the next sample_df sample
	# 			  	3b) If the current sample_df entry lacks a TrappingDate, print its mchip and continue (this should not occur) <- THIS SHOULD MAYBE RAISE AN EXCEPTION
	# 			  4) Obtain the absolute value of the difference in days between the tumor_df and sample_df match TrappingDate
	# 			  5) Attempt to match the sample_df TumourID (if it exists) with the tumor_df TumourID
	# 			  6) Determine the index of tumor_df TrappingDate which is closest to the sample_df TrappingDate
	# 				 If specific_tumor=True, make this single index what's extracted from the tumor_df. Otherwise, make all matching mchip/TumourNumbers the extracted indices
	# 			  7) If the tumor_df and sample_df TrappingDates are different, save the index of both the sample_df and a single tumor_df to separate lists.
	# 				 If closest_date=False, the matching TumourID sample TrappingDate from tumor_df is used (and if specific_tumor=True, this becomes the extracted tumor).
	# 			  8) If the difference in TrappingDates is within date_tolerance, the matching tumor(s) is/are extracted
	# 			  9) If id_match=True and a matching TumourID was found, the tumor(s) is/are extracted (overrides date_tolerance)
	# 			  10) Otherwise, no tumor is extracted from the tumor_df
	# 			  === LOOP END ===
	# 			  11) Extracted tumors are pulled from tumor_df via index, unextracted tumors are pulled from sample_df via index, and the date mismatches (even for tumors which
	# 			  	  were extracted) are pulled both from sample_df and tumor_df and are then merged. If closest_date=True, the *single* tumor_df sample, and corresponding date_diff,
	# 			  	  represent the tumor_df closest in TrappingDate to the sample_df TrappingDate. Otherwise, if a matching TumourID is found, the tumor_df sample and date_diff
	# 				  represent the matching TumourID tumor_df sample
	# 			  12) Print the failed date_tolerance extracts and myriad extract stats if stats=True
	# GOTCHAS:
	# 			  - If a perfectly matching TrappingDate is found in tumor_df, this is considered a true match even if the TumourIDs differ
	# 			  - Toggling closest_date=False can result in failed extracts if id_match=False which would otherwise be successful. This only
	# 				occurs when the tumor_df has multiple TrappingDates and a matching TumourID with a date is further than the sample_df date:
	# 				=== closest_tumor=True, id_match=False ===
	# 				Removing 18 samples with differences in dates greater than the tolerance (date_tolerance=50)
    # 					  Microchip TumourNumber TumourID_sample TrapDate_sample TumourID_DB TrapDate_DB  date_diff
	# 				982009106029805            1             610      2012-03-01         708  2012-03-06          5
	# 				=== closest_tumor=False, id_match=False ===
	# 				Removing 33 samples with differences in dates greater than the tolerance (date_tolerance=50)
    # 					  Microchip TumourNumber TumourID_sample TrapDate_sample TumourID_DB TrapDate_DB  date_diff
	# 				982009106029805            1             610      2012-03-01         610  2011-11-10        112
	# TODO:
	# 			  - See Step 3b about missing sample_df TrappingDates raising an exception
	# 			  - This is a big, unwieldy mess. Try and clean it up some
	# 			  - Optimize, the loop is inefficient (time tests to see if this is really inefficient and thus necessary)
	def get_sample_tumors(self, date_tolerance=50, id_match=False, closest_date=True, specific_tumor=False, stats=False, inplace=True):
		self.date_tolerance = date_tolerance
		sample_subset = self.sample_df[self.sample_df["Microchip"].isin(self.tumor_df["Microchip"])].copy()
		sample_subset["TumourNumber"] = sample_subset["TumourNumber"].str.replace(r"[a-zA-Z]$", "", regex=True)
		sample_tumors = sample_subset[["Microchip", "TumourNumber", "TumourID", "TrappingDate"]]
		found_indices = []
		not_found_indices = []
		date_mismatch_tumor = []
		date_mismatch_sample = []
		mismatch_vals = []
		nan_traps = []
		failed_dates = 0
		for i in range(sample_tumors.shape[0]):
			current_sample = sample_tumors.iloc[i]
			found_tumor = self.tumor_df.loc[(self.tumor_df["Microchip"] == current_sample["Microchip"]) & (self.tumor_df["TumourNumber"] == current_sample["TumourNumber"])]
			if found_tumor.size == 0:
				not_found_indices.append(current_sample.name)
				continue
			if pd.isnull(current_sample["TrappingDate"]):
				print(current_sample["Microchip"])
				continue
			date_diff = abs(current_sample["TrappingDate"] - found_tumor[["TrapDate"]])
			diff_val = date_diff["TrapDate"].min().days
			tumor_id_match = found_tumor[found_tumor["TumourID"] == current_sample["TumourID"]]
			closest_date_index = [date_diff.idxmin()[0]]
			if specific_tumor:
				found_index = closest_date_index
			else:
				found_index = found_tumor.index.tolist()
			if diff_val > 0:
				if tumor_id_match.size > 0 and not closest_date:
					closest_date_index = tumor_id_match.index.tolist()
					if specific_tumor:
						found_index = closest_date_index
					diff_val = abs(tumor_id_match["TrapDate"] - current_sample["TrappingDate"]).iloc[0].days
				date_mismatch_tumor += closest_date_index
				date_mismatch_sample.append(current_sample.name)
				mismatch_vals.append(diff_val)
			if not date_tolerance or diff_val <= date_tolerance:
				found_indices += found_index
			elif id_match and tumor_id_match.size > 0:
				found_indices += found_index
			else:
				failed_dates += 1
		extracted = self.tumor_df.loc[found_indices]
		unextracted_tumors = self.sample_df.loc[not_found_indices]
		mismatch_tumor = self.tumor_df.loc[date_mismatch_tumor][["Microchip", "TumourNumber", "TumourID", "TrapDate"]]
		mismatch_sample = sample_subset.loc[date_mismatch_sample][["Microchip", "TumourNumber", "TumourID", "TrappingDate"]]
		tumor_date_mismatches = mismatch_sample.merge(mismatch_tumor, on=["Microchip", "TumourNumber"])
		tumor_date_mismatches = tumor_date_mismatches.rename(columns={"TumourID_x": "TumourID_sample", "TrappingDate": "TrapDate_sample", "TumourID_y": "TumourID_DB", "TrapDate": "TrapDate_DB"})
		tumor_date_mismatches["date_diff"] = mismatch_vals
		if date_tolerance and inplace == True:
			print(f"Removing {failed_dates} tumor_df samples with differences in dates greater than the tolerance (date_tolerance={date_tolerance})")
		if stats:
			mismatch_days = tumor_date_mismatches["date_diff"]
			not_extracted = unextracted_tumors["Microchip"].unique()[:3]
			total_mismatch = tumor_date_mismatches[tumor_date_mismatches["TumourID_sample"] != tumor_date_mismatches["TumourID_DB"]]
			to_extract_n = len(sample_subset["Microchip"].unique())
			extracted_n = len(extracted["Microchip"].unique())
			print(f"To extract: {to_extract_n}")
			print(f"Extracted: {extracted_n}")
			print(f"Not extracted: {to_extract_n - extracted_n}")
			print(f"Failed TumourNumbers: {len(unextracted_tumors['Microchip'].unique())}")
			print(f"NaN TrappingDate (sample_df): {len(nan_traps)}")
			print(f"Number of matching trap dates: {extracted_n - len(mismatch_days)}")
			print(f"Number of date mismatches: {len(mismatch_days)}")
			print(f"Date mismatches failing tolerance (day_diff > {date_tolerance}): {len(mismatch_days[mismatch_days > date_tolerance])}")
			print(f"Number of date mismatches with no tumourID: {len(total_mismatch)}")
			print(f"Mean date mismatch (days): {round(sum(mismatch_days)/len(mismatch_days))}")
			print(f"Median date mismatch (days): {np.median(mismatch_days)}")
			print(f"Failed mchips (first 3): {','.join(not_extracted)}")
			print("""NOTE: Adding up failed TumourNumbers + date mismatches may not equal \"Not Extracted\" \n      if a devil with two tumors had one tumor extract and the other fail""")
		if not inplace:
			return extracted, unextracted_tumors, tumor_date_mismatches
		self.sample_tumor_flag = True
		self.tumor_df = extracted
		self.unextracted_tumors = unextracted_tumors
		self.tumor_date_mismatches = tumor_date_mismatches


	def subset_extracted(self, inplace=True):
		extracted_mchips = self.tumor_df["Microchip"].unique()
		subset = self.sample_df[self.sample_df["Microchip"].isin(extracted_mchips)]
		if not inplace:
			return subset
		print(f"*subset_extracted()* Removing {self.sample_df.shape[0] - subset.shape[0]} rows out of {self.sample_df.shape[0]} ({subset.shape[0]} remaining)")
		self.sample_df = subset


	def subset_trap_n(self, trap_N, greater=True, inplace=True):
		unique_traps = self.tumor_df.groupby("Microchip")["TrapDate"].nunique()
		if not greater:
			unique_traps = unique_traps[unique_traps == trap_N].index
		else:
			unique_traps = unique_traps[unique_traps > trap_N].index
		df = self.tumor_df[self.tumor_df["Microchip"].isin(unique_traps)]
		print(f"*subset_trap_n()* Removing {self.tumor_df.shape[0] - df.shape[0]} rows out of {self.tumor_df.shape[0]} ({df.shape[0]} remaining)")
		if not inplace:
			return df
		self.tumor_df = df


	def subset_unextracted(self, inplace=True):
		unextracted_mchips = self.tumor_df["Microchip"].unique()
		subset = self.sample_df[~self.sample_df["Microchip"].isin(unextracted_mchips)]
		if not inplace:
			return subset
		print("*WARNING* This will break all functionality of the TumorSamples class!")
		self.sample_df = subset
	# =============================================== Subsetting END ==============================================
	# =============================================================================================================



	# =============================================== Reporting START =============================================
	# =============================================================================================================
	def fast_stats_tumor(self):
		tissue = self.sample_df["Tissue"].unique()
		if(len(tissue) > 1):
			tissue = "Tumors and Devils"
		else:
			tissue = f"{tissue[0]}s"
		no_db = self.sample_df[~self.sample_df["Microchip"].isin(self.tumor_df["Microchip"].unique())].shape[0]
		in_db = self.sample_df[self.sample_df["Microchip"].isin(self.tumor_df["Microchip"].unique())].shape[0]
		unique_samples = len(self.sample_df["Microchip"].unique())
		unique_tumor_finds = len(self.tumor_df["Microchip"].unique())
		in_tumor_db = self.sample_df[self.sample_df["Microchip"].isin(self.tumor_df["Microchip"].unique())].shape[0]
		sample_n = self.sample_df.shape[0]
		chip_date_cluster = self.tumor_df.groupby(["Microchip", "TrapDate"])["TrapDate"].count()
		trap_cluster = chip_date_cluster.groupby("Microchip").count()
		traps = self.tumor_df.groupby("Microchip")["TrapDate"].count()
		tumor_percent = (in_tumor_db/sample_n)*100
		tumor_percent_uniq = (unique_tumor_finds/unique_samples)*100

		print(f"===== {tissue} =====")
		print(f"Found in tumor DB: {tumor_percent:.1f}% ({in_tumor_db}/{sample_n})")
		print(f"Found in tumor DB (unique mchips): {tumor_percent_uniq:.1f}% ({unique_tumor_finds}/{unique_samples})")
		print(f"Failed microchips: {self.failed_chips.shape[0]}")
		if self.drop_DFTG:
			print(f"Individual tumors removed (min DFTG={self.drop_DFTG}): {self.removed_tumors['Microchip'].unique().shape[0]}")
			print(f"Samples removed with DFTG < {self.drop_DFTG}: {self.removed_samples['Microchip'].unique().shape[0]}")
		print(f"1 trap: {len(trap_cluster[trap_cluster == 1])}")
		print(f"More than 1 trap: {len(trap_cluster[trap_cluster > 1])}")

		if tumor_percent > 100 or tumor_percent_uniq > 100:
			print("!WARNING! Tumor DB extraction rate is above 100%, is the tumor DF out of sync? Try sync_tumor_df()")


	def get_date_mismatches(self, failed_tolerance=True):
		if self.tumor_date_mismatches.empty:
			print("Please run get_sample_tumors() before using this method")
			exit()
		if failed_tolerance:
			return self.tumor_date_mismatches[self.tumor_date_mismatches["date_diff"] > self.date_tolerance]
		return self.tumor_date_mismatches


	def get_DFTG_removed(self, remove_type="sample", cols=["Microchip", "TrapDate", "TumourNumber", "TumourDFTG"]):
		if not cols:
			cols = self.tumor_df.columns
		if remove_type == "sample":
			return self.removed_samples[cols]
		elif remove_type == "tumor":
			return self.removed_tumors[cols]
		else:
			print(f"{remove_type} not recognized as a remove_type! Please use: \"sample\" or \"tumor\"")
			return None


	def get_multi_sub_40(self):
		if len(self.failed_date_delta) == 0:
			print("Please run survival_proxy() before using this method")
			return None
		single_trap = self.subset_trap_n(1, greater=False, inplace=False)
		return self.failed_date_delta[~self.failed_date_delta.isin(single_trap["Microchip"].unique())]


	def get_unextracted_tumors(self):
		if self.unextracted_tumors.empty:
			print("Please run get_sample_tumors() before using this method")
			exit()
		return self.unextracted_tumors


	def no_tumor_entry(self):
		hosts = self.sample_df[self.sample_df["Tissue"] == "Host"]
		hosts_not_in = hosts[~hosts["Microchip"].isin(self.tumor_df["Microchip"].unique())]
		return hosts_not_in


	def print_over_mmax(self):
		if not self.over_mmax["Microchip"]:
			print("No back calculations have been done yet")
			return None
		print("Microchip\tVolume (cm^3)")
		N = len(self.over_mmax["Microchip"])
		for i in range(N):
			print(f"{self.over_mmax['Microchip'][i]}\t{self.over_mmax['volume'][i]:.3f}")
		print(f"Length: {N}")


	def print_smaller_obs(self):
		if not self.smaller_obs:
			print("No back calculations have been done yet")
			return None
		print(self.smaller_obs.to_string(index=False))
	# ================================================ Reporting END ==============================================
	# =============================================================================================================



	# ================================================ Utility START ==============================================
	# =============================================================================================================
	# Definition: 
	# 			  Sync the tumor DF with the sample DF. This is only needed when an update is made to the sample DF (i.e., a subsetting operation)
	# Arguments:
	# 			  None: Pass
	# Returns:
	# 			  None. Updates the tumor DF so its samples match those in the sample DF
	# Steps:
	# 			  1) Obtain the unique sample DF mchips
	# 			  2) Save the number of tumor DF entries before syncing
	# 			  3) Update the tumor DF to include only those samples with mchips matching the sample DF's mchips
	# 			  4) Print the number of tumor samples removed and the number of remaining tumor samples
	def sync_tumor_df(self):
		sample_mchip = self.sample_df["Microchip"].unique()
		before = self.tumor_df.shape[0]
		self.tumor_df = self.tumor_df[self.tumor_df["Microchip"].isin(sample_mchip)]
		after = self.tumor_df.shape[0]
		print(f"*sync_tumor_df()* Removed {before-after} tumor samples ({after} tumor samples remaining)")


	def __drop_DFTG(self, df, min_score):
		df = df.copy()
		df.loc[df["TumourDFTG"].isna(), "TumourDFTG"] = 0
		df.loc[df["TumourNumber"].isna(), "TumourNumber"] = "n"
		passed_indices = []
		removed_indices = []
		for chip in df["Microchip"].unique():
			sample = df[df["Microchip"] == chip]
			for tumor_num in sample["TumourNumber"].unique():
				tumor = sample[sample["TumourNumber"] == tumor_num]
				max_DFTG = tumor["TumourDFTG"].max()
				if max_DFTG >= min_score:
					passed_indices += tumor.index.values.tolist()
				else:
					removed_indices += tumor.index.values.tolist()
		removed_tumors = df.loc[removed_indices]
		cleaned_df = df.loc[passed_indices]
		cleaned_df["TumourNumber"] = cleaned_df["TumourNumber"].replace("n", np.nan)
		df_chips = df["Microchip"].unique()
		cleaned_chips = cleaned_df["Microchip"].unique()
		bool_mask = np.in1d(df_chips, cleaned_chips)
		removed_sample_chips = df_chips[~bool_mask]
		removed_samples = df[df["Microchip"].isin(removed_sample_chips)]
		return cleaned_df, removed_samples, removed_tumors


	# Definition: 
	# 			  Convenience method to initialize variables
	# Variables:
	# 			  survival_proxy():
	# 			  	failed_date_delta
	# 			  	over_mmax
	# 				smaller_obs
	# 			  get_sample_tumors():
	# 				date_tolerance
	# 				sample_tumor_flag
	# 			  	tumor_date_mismatches
	# 				unextracted_tumors
	def __init_vars(self):
		self.failed_date_delta = []
		self.over_mmax = {"Microchip": [], "volume": []}		
		self.smaller_obs = None
		self.date_tolerance = None
		self.sample_tumor_flag = False
		self.tumor_date_mismatches = pd.DataFrame([])
		self.unextracted_tumors = pd.DataFrame([])
	# ================================================= Utility END ===============================================
	# =============================================================================================================



	# ============================================== In Progress START ============================================
	# =============================================================================================================	
	# Might get rid of this
	def update_dftd_score(self, csv_paths, date_tolerance=3, DFTG_min=5):
		if self.drop_DFTG:
			print("Updating DFTD scores (TumourDFTG) can only be done if a TumorSamples object is instantiated with drop_DFTG=False!")
			exit()
		if not self.full_tumor_df:
			print("Updating DFTD scores (TumourDFTG) can only be done if a TumorSamples object is instantiated with full_tumor_df=True!")
			exit()
		extracted_tumors = self.tumor_df[self.tumor_df["Microchip"].isin(self.sample_df["Microchip"].unique())]
		_, dropped_samples, _ = self.__drop_DFTG(extracted_tumors, DFTG_min)
		search_df = [pd.read_csv(path) for path in csv_paths]
		score_cols = ["DFTDScore", "DFTDscore", "Score"]
		for i in range(len(search_df)):
			df = search_df[i]
			colname = [col for col in score_cols if col in df.columns][0]
			search_df[i] = df.rename(columns={colname: "DFTDScore"})
			search_df[i] = search_df[i][["Microchip", "TrappingDate", "DFTDScore"]]
		search_df = pd.concat(search_df)
		search_df = search_df[search_df["Microchip"].isin(dropped_samples["Microchip"])]
		search_df["TrappingDate"] = pd.to_datetime(search_df["TrappingDate"])
		found = {}
		failed = []
		failed_dates = []
		for i in range(dropped_samples.shape[0]):
			tumor = dropped_samples.iloc[i]
			found_samples = search_df[search_df["Microchip"] == tumor["Microchip"]]
			if found_samples.empty:
				failed.append(tumor.name)
				continue
			date_diff = abs(found_samples["TrappingDate"] - tumor["TrapDate"])
			closest_date = min(date_diff).days
			if closest_date > date_tolerance:
				failed_dates.append(tumor.name)
				continue
			closest_date_i = date_diff.idxmin()
			num_samples = found_samples.loc[closest_date_i]
			num_samples = found_samples[found_samples["TrappingDate"] == num_samples["TrappingDate"]].shape[0]
			found[closest_date_i] = num_samples
		found_score = search_df.loc[found.keys()]
		found_score["number_of_samples"] = found.values()
		found_score = found_score.drop_duplicates(subset="Microchip")
		no_chip = dropped_samples.loc[failed]
		no_date = dropped_samples.loc[failed_dates]
		print(found_score.shape[0])
		print(no_chip.shape[0])
		print(no_date.shape[0])
		print(len(dropped_samples["Microchip"].unique()))
	# =============================================== In Progress END =============================================
	# =============================================================================================================





class Compare:
	def __init__(self, df, extract_path, sep="\t", id="Microchip"):
		extract_df = pd.read_csv(extract_path, sep=sep)
		extract_df = extract_df.rename(columns={"mchip": "Microchip"})
		merged_df = df.merge(extract_df, how="inner", on=id)
		self.id = id
		self.df = merged_df
		self.matched = pd.Series()
		self.unmatched = pd.Series()


	def compare_strict(self, col1, col2):
		self.matched = self.df[self.df[col1] == self.df[col2]][[self.id, col1, col2]]
		self.unmatched = self.df[self.df[col1] != self.df[col2]][[self.id, col1, col2]]
		self.comp_cols = [col1, col2]


	def df_print(self, cols):
		print(self.df[cols])


	def diff(self, n=1, biggest=True):
		first = self.unmatched[self.comp_cols[0]]
		second = self.unmatched[self.comp_cols[1]]
		abs_diffs = sorted((first - second).abs(), reverse=biggest)
		return abs_diffs[0:n]


	def diff_stats(self):
		first = self.unmatched[self.comp_cols[0]]
		second = self.unmatched[self.comp_cols[1]]
		abs_diffs = (first - second).abs()
		mean = abs_diffs.mean()
		median = abs_diffs.median()
		return mean, median


	def show_diff(self, print_only=True):
		first = self.unmatched[self.comp_cols[0]]
		second = self.unmatched[self.comp_cols[1]]
		abs_diffs = (first - second).abs()
		diff_frame = pd.DataFrame({
			"Microchip": self.unmatched["Microchip"],
			self.comp_cols[0]: first,
			self.comp_cols[1]: second,
			"diff": abs_diffs
			})
		if print_only:
			print(diff_frame)
		else:
			return diff_frame


	def stats(self):
		if self.matched.empty and self.unmatched.empty:
			print("Please run a comparison method before printing stats!")
			return None
		mean, median = self.diff_stats()
		print(f"Matched: {self.matched.shape[0]}")
		print(f"Unmatched: {self.unmatched.shape[0]}")
		print("----- UNMATCHED -----")
		print(f"Biggest difference: {self.diff()[0]}")
		print(f"Smallest difference: {self.diff(biggest=False)[0]}")
		print(f"Mean difference: {mean:.1f}")
		print(f"Median difference: {median}")


	def subset_non_nan(self, col):
		self.df = self.df[~self.df[col].isna()]