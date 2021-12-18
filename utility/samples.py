import os
import re
import time
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



class Samples:
	base = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data"
	seq_base = f"{base}/data/sequencing"
	batch_ids = [f"{seq_base}/Capture1/rename_key.csv", f"{seq_base}/Capture2/rename_key.csv", f"{seq_base}/Capture3/rename_key.csv", f"{seq_base}/Capture4/rename_key.csv", f"{seq_base}/Capture5/rename_key.csv"]
	prelim = batch_ids[:2]
	dftd_site_arrival = pd.DataFrame({
		"Site": ["Freycinet", "WPP", "Takone", "Black River", "Arthur River"],
		"arrival_date": [1999, 2007, 2013, 2015, 2019],
		"generations": [11, 7, 4, 3, 1],
		"years": [22, 14, 8, 6, 2]
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


	def add_arrival(self, arrival_cols=["arrival_date"]):
		arrivals = Samples.dftd_site_arrival[["Site"] + arrival_cols]
		self.sample_df = self.sample_df.merge(arrivals, how="left", on="Site")


	def add_col(self, col):
		name = col.name if col.name else "new_col"
		self.sample_df[name] = col


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


	def extract_mchips(self, chip_list):
		found = []
		for chip in chip_list:
			last_six = chip[-6:]
			samples = self.sample_df[self.sample_df["Microchip"].str[-6:] == last_six]
			if chip[0] == "H" or chip[0] == "T":
				tissue = "Host" if chip[0] == "H" else "Tumour"
				samples = samples[samples["Tissue"] == tissue]
			found.append(samples)
		return pd.concat(found)


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
		pairs = len(self.host_tumor_pairs(inplace=False)["Microchip"])
		singles = self.sample_df.shape[0] - pairs
		multi = len(self.host_tumor_multi()["Microchip"])
		print(f"Samples: {self.sample_df.shape[0]}")
		print(f"Hosts: {self.sample_df[self.sample_df['Tissue'] == 'Host'].shape[0]}")
		print(f"Tumors: {self.sample_df[self.sample_df['Tissue'] == 'Tumour'].shape[0]}")
		print(f"Missing microchips: {self.nan_sum('Microchip', print_flag=False)}")
		print(f"Unique microchips: {len(self.sample_df['Microchip'].unique())}")
		print(f"Singles: {singles}")
		print(f"Pairs: {pairs}")
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
		print(f"Removing {self.sample_df.shape[0] - pairs.shape[0]} rows out of {self.sample_df.shape[0]}")
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
		methods = ["inv_logit", "log", "logit", "standard"]
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


	def subset(self, col, pattern, inplace=True):
		subset = self.sample_df[self.sample_df[col] == pattern]
		if not inplace:
			return subset
		self.sample_df = subset


	def subset_pairs(self):
		indices = self.get_pairs()
		indices = indices[indices > 1].index
		self.sample_df = self.sample_df[self.sample_df["Microchip"].isin(indices)]


	def subset_non_nan(self, col, verbose=True):
		na = self.sample_df[col].isna()
		if isinstance(col, list):
			na = na.any(axis=1)
		if verbose:
			print(f"Removing {na.values.sum()} rows out of {self.sample_df.shape[0]}")
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


	def ATOMM_pheno(self, outpath, cols, to_pheno=True):
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
		if to_pheno:
			self.to_pheno(outpath, cols, id_cols=["host_chip", "tumor_chip"], vcf_chip=False)


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
				print(f"!WARNING! sample_df has no entries! Tissue subsetted to \"{tissue}\", ensure spelling and capitalization is correct")
		microchips = self.sample_df["Microchip"].unique()
		if extract_on != "both" and tissue == "both":
			microchips = self.sample_df[self.sample_df["Tissue"] == extract_on]["Microchip"].unique()
		if not full_tumor_df:
			self.tumor_df = tumor_df_full[tumor_df_full["Microchip"].isin(microchips)]
		else:
			self.tumor_df = tumor_df_full
		self.drop_DFTG = drop_DFTG
		if drop_DFTG:
			self.tumor_df, self.removed_samples, self.removed_tumors = self.__drop_DFTG(self.tumor_df, drop_DFTG)
		self.over_mmax = {"Microchip": [], "volume": []}
		self.failed_date_delta = []


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

	
	# ANALYSIS PHENOTYPE
	def estimate_age(self, inplace=True, verbose_merge=False):
		tmp_tumor_df = self.tumor_df[["Microchip", "TrapDate", "TumourDepth", "TumourLength", "TumourWidth"]].copy()
		dates_df = self.__init_tumor_date(tmp_tumor_df)
		YOB = self.sample_df[self.sample_df["Tissue"] == "Host"].set_index("Microchip")["YOB"]
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


	def fast_stats_tumor(self):
		tissue = self.sample_df["Tissue"].unique()
		if(len(tissue) > 1):
			tissue = "Tumors and devils"
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
		print(f"Individual tumors removed (min DFTG={self.drop_DFTG}): {self.removed_tumors['Microchip'].unique().shape[0]}")
		print(f"Samples removed with DFTG < {self.drop_DFTG}: {self.removed_samples['Microchip'].unique().shape[0]}")
		print(f"1 trap: {len(trap_cluster[trap_cluster == 1])}")
		print(f"More than 1 trap: {len(trap_cluster[trap_cluster > 1])}")

		if tumor_percent > 100 or tumor_percent_uniq > 100:
			print("!WARNING! Tumor DB extraction rate is above 100%, is the tumor DF out of sync? Try sync_tumor_df()")


	# === ANALYSIS PHENOTYPE ===
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
	# 			  The calculated dates often do not match up with Margres et al. (2018) table S1.
	# 			  Attempt to change the calculation using the sum of all tumor volumes on the min trap date.
	# 			  Use the Compare class to test calculation similarities
	def survival_proxy(self, day_cutoff=40, proxy=True, verbose_merge=False):
		tmp_tumor_df = self.tumor_df[["Microchip", "TrapDate", "TumourDepth", "TumourLength", "TumourWidth"]].copy()
		tmp_tumor_df["TrapDate"] = pd.to_datetime(tmp_tumor_df["TrapDate"])
		
		date_diff = tmp_tumor_df.groupby("Microchip")["TrapDate"].transform(lambda x: x.max() - x.min())
		tmp_tumor_df["date_diff"] = date_diff.dt.days
		self.failed_date_delta = pd.Series(tmp_tumor_df[tmp_tumor_df["date_diff"] < day_cutoff]["Microchip"].unique())
		tmp_tumor_df = tmp_tumor_df[tmp_tumor_df["date_diff"] >= day_cutoff]

		dates_df = self.__init_tumor_date(tmp_tumor_df, proxy)
		last_trap = tmp_tumor_df.groupby("Microchip")["TrapDate"].max().loc[dates_df.index]
		devil_survival_days = last_trap - dates_df["init_tumor_date"]
		
		calc_df = pd.DataFrame({
			"Microchip": dates_df.index.values,
			"volume (cm^3)": np.round(dates_df["tumor_volume"].values, 2),
			"last_trap": last_trap.dt.date.values,
			"init_tumor_date": dates_df["init_tumor_date"].dt.date.values,
			"devil_survival_days": devil_survival_days.dt.days.values})
		if not proxy:
			calc_df.insert(loc=2, column="back_calc", value=dates_df["back_calc"].values)
			calc_df.insert(loc=3, column="first_trap", value=dates_df["min_date"].dt.date.values)
		if not verbose_merge:
			calc_df = calc_df[["Microchip", "devil_survival_days"]]
		self.sample_df = self.sample_df.merge(calc_df, how="left", on="Microchip")


	def get_multi_sub_40(self):
		if len(self.failed_date_delta) == 0:
			print("Please run survival_proxy() before using this method")
			return None
		single_trap = self.get_trap(1)
		return self.failed_date_delta[~self.failed_date_delta.isin(single_trap)]


	def get_trap(self, trap_N):
		trap_nums = self.tumor_df.groupby(["Microchip", "TrapDate"])["TrapDate"].count().groupby("Microchip").count()
		return trap_nums[trap_nums == trap_N].index


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


	def get_sample_tumors(self, stats=False, date_tolerance=100):
		sample_subset = self.sample_df[self.sample_df["Microchip"].isin(self.tumor_df["Microchip"])]
		sample_tumors = sample_subset[["Microchip", "TumourNumber", "TrappingDate"]]
		found_indices = []
		date_mismatch = {}
		nan_traps = []
		for i in range(sample_tumors.shape[0]):
			current_sample = sample_tumors.iloc[i]
			found_index = self.tumor_df.loc[(self.tumor_df["Microchip"] == current_sample["Microchip"]) & (self.tumor_df["TumourNumber"] == current_sample["TumourNumber"])]
			if found_index.shape[0] > 1:
				date_diff = current_sample["TrappingDate"] - found_index[["TrapDate"]]
				diff_val = abs(date_diff["TrapDate"].min().days)
				if diff_val > 0:
					date_mismatch[current_sample["Microchip"]] = diff_val
				found_index = [date_diff.idxmin()[0]]
				if np.isnan(found_index):
					nan_traps.append(current_sample["Microchip"])
					continue
			else:
				found_index = list(found_index.index.values)
			found_indices += found_index
		extracted = self.tumor_df.loc[found_indices]
		for k,v in date_mismatch.items():
			print(f"{k}: {v}")
		if stats:
			mismatch_days = date_mismatch.values()
			not_extracted = sample_subset[~sample_subset["Microchip"].isin(extracted["Microchip"])]
			not_extracted = not_extracted["Microchip"].unique()[:3]
			print(f"To extract: {sample_subset.shape[0]}")
			print(f"Extracted: {extracted.shape[0]}")
			print(f"Not extracted: {sample_subset.shape[0] - extracted.shape[0]}")
			print(f"NaN TrappingDate (sample_df): {len(nan_traps)}")
			print(f"Number of date mismatches: {len(mismatch_days)}")
			print(f"Mean date mismatch (days): {round(sum(mismatch_days)/len(mismatch_days))}")
			print(f"Failed mchips: {','.join(not_extracted)}")
		return extracted


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


	def subset_extracted(self, inplace=True):
		extracted_mchips = self.tumor_df["Microchip"].unique()
		subset = self.sample_df[self.sample_df["Microchip"].isin(extracted_mchips)]
		if not inplace:
			return subset
		self.sample_df = subset


	def subset_unextracted(self, inplace=True):
		unextracted_mchips = self.tumor_df["Microchip"].unique()
		subset = self.sample_df[~self.sample_df["Microchip"].isin(unextracted_mchips)]
		if not inplace:
			return subset
		print("*WARNING* This will break all functionality of the TumorSamples class!")
		self.sample_df = subset


	def sync_tumor_df(self):
		sample_mchip = self.sample_df["Microchip"].unique()
		before = self.tumor_df.shape[0]
		self.tumor_df = self.tumor_df[self.tumor_df["Microchip"].isin(sample_mchip)]
		after = self.tumor_df.shape[0]
		print(f"Removed {before-after} tumor samples ({after} tumor samples remaining)")


	# === ANALYSIS PHENOTYPE ===
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
	# 			  None: Template
	# Returns:
	# 			  None. Adds a "transmission_age" column to sample_df
	# Steps:
	# 			  1) Obtain the infection_age and arrival_date if it was not already added to sample_df
	# 			  2) Calculate arrival - YOB, where positive values indicate a devil born X days prior to DFTD arrival
	# 			  3) For any arrival - YOB less than min infection_age, set the value to 0
	# 			  4) Make a min infection_age Series, setting any value to 0 which is not 0 in the arrival - YOB
	# 				 NOTE: This and 3) ensure that only YOB - arrival OR min infection_age (whichever is larger) is subtracted from infection_age
	# 			  5) Subtract YOB - arrival OR min infection_age from infection age (this removes the non-viable devil infection period)
	# 			  6) Restore the original sample_df (which may not have had estimate_age or add_survival)
	# Gotchas:
	# 			  Devil 982009104872893 has a negative transmission_age (check for more when I have >312 samples extracted)
	# 
	# 			  The correction ALWAYS generates at least a single 0 value. When log transforming, this will be -inf. For now, manually correct
	# 			  this to a sensible value smaller than all other values
	# 
	# 			  Tumor transmission_age values are also often nonsense. This is expected, ensure these are ignored for tumors (only hosts can have infection/transmission age)
	def transmission(self):
		original_df = self.sample_df.copy()
		if not "infection_age" in self.sample_df.columns:
			self.estimate_age()
		if not "arrival_date" in self.sample_df.columns:
			self.add_arrival()
		min_age = self.sample_df["infection_age"].min()
		born_before_arrival = (self.sample_df["arrival_date"] - self.sample_df["YOB"]).dt.days.copy()
		born_before_arrival[born_before_arrival < min_age] = 0
		min_series = pd.Series([min_age] * len(born_before_arrival))
		min_series.loc[born_before_arrival[born_before_arrival != 0].index] = 0
		transmission_age = self.sample_df["infection_age"] - born_before_arrival - min_series
		self.sample_df = original_df
		self.sample_df["transmission_age"] = transmission_age


	# ANALYSIS PHENOTYPE
	def tumor_count(self):
		num_tumors = self.tumor_df.groupby(["Microchip"])["TumourNumber"].unique().str.len()
		num_tumors = num_tumors.to_frame().reset_index()
		num_tumors = num_tumors.rename(columns={"TumourNumber": "tumor_count"})
		self.sample_df = self.sample_df.merge(num_tumors, how="left", on="Microchip")


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
		df_chips = df["Microchip"].unique()
		cleaned_chips = cleaned_df["Microchip"].unique()
		bool_mask = np.in1d(df_chips, cleaned_chips)
		removed_sample_chips = df_chips[~bool_mask]
		removed_samples = df[df["Microchip"].isin(removed_sample_chips)]
		return cleaned_df, removed_samples, removed_tumors


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
	# 			  Formula domain: (0, mmax)
	# 			  Formula range: (-inf, 0); as volume -> mmax, days -> -inf
	def __growth_back_calculation(self, tumor_volumes, mmax=202, init_size=0.0003, alpha=0.03, mmax_cutoff=100, is_cm=False, round_flag=True):
		tumor_volumes = tumor_volumes.copy()
		if not is_cm:
			tumor_volumes /= 1000
		over = tumor_volumes[(tumor_volumes + 1) >= mmax]
		if over.size > 0:
			vstring = "volumes" if len(over) > 1 else "volume"
			print(f"*WARNING* found {len(over)} {vstring} greater than mmax while performing the back calculation (try print_over_mmax() for details)")
			acceptable_range = over[over <= mmax + mmax_cutoff]
			print(f"{acceptable_range.size} of these are within the {mmax_cutoff} cutoff, setting these to volume=200 cm^3")
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
	# 			  Conditionally drop rows if any NA is found in TumourDepth, TumourLength, or TumourWidth. If this is not done, then volumes with a missing
	# 			  dimension are NaN. This is fine for the proxy, but will throw an error for the back calculation estimate.
	# TODO:
	# 			  The loop is an inelegant and inefficient solution, this should be vectorized
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


	# Definition:
	# 			  This is a convenience function doing one of two things:
	# 			  1) Combine the functionality of __min_cap_volume() and __growth_back_calculation() to find the initial tumor date (back-calc estimate)
	# 			  2) Find the minimum trapping date for which a tumor was observed (survival proxy)
	# Arguments:
	# 			  tumor_df:   a DF with at least the following cols: Microchip, TrapDate, TumourDepth, TumourLength, TumourWidth <pd.DataFrame>
	# 			  proxy_flag: True if using a basic for the init date (i.e., first date trapped with a tumor present). Otherwise, employ the back-calc to estimate the init tumor date
	# Returns:
	# 			  a DF with Microchip as the index and one of two column configurations:
	# 			  1) tumor_volume, min_date, back_calc, init_tumor_date
	# 			  2) tumor_volume, init_tumor_date
	# Steps:
	# 			  1: For __min_cap_volume(), dropna should be true when proxy_flag is false
	# 			  2) Run __min_cap_volume() to obtain the minimum trap date and max tumor volume for that date
	# 			  BACK CALCULATION
	# 			  3) Run __growth_back_calculation() to obtain days since init tumor (negative days)
	# 			  4) Get the initial tumor date by subtracting the back calc from the min date
	# 			  5) Combine everything into a DF
	# 			  PROXY
	# 			  3) Rename "min_date" (minimum trap date obtained from __min_cap_volume) to "init_tumor_date"
	def __init_tumor_date(self, tumor_df, proxy_flag):
		dropna_flag = not proxy_flag
		volumes = self.__min_cap_volume(tumor_df, dropna_flag)
		volumes["tumor_volume"] /= 1000
		if not proxy_flag:
			back_calc = self.__growth_back_calculation(volumes["tumor_volume"], is_cm=True)
			init_tumor_dates = volumes["min_date"] - pd.to_timedelta(-1 * back_calc, unit="d")
			init_tumor_dates.name = "init_tumor_date"
			back_calc = pd.Series(back_calc, name="back_calc", index=volumes.index)
			final_df = pd.concat([volumes, back_calc, init_tumor_dates], axis=1)
		else:
			final_df = volumes.rename(columns={"min_date": "init_tumor_date"})
		return final_df



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