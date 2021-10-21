import os
import re
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



class Samples:
	base = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data"
	batch_ids = [f"{base}/Capture1_6-11-21/rename_key.csv", f"{base}/Capture2_7-29-21/rename_key.csv", f"{base}/Capture3/rename_key.csv", f"{base}/Capture4/rename_key.csv", f"{base}/Capture5/rename_key.csv"]
	prelim = batch_ids[:2]


	def __init__(self,
		sample_csv_path="/work_bgfs/d/dgallinson/data/pheno_data/master_corrections.csv",
		id_paths=[]
		):
		full_df = pd.read_csv(sample_csv_path)
		if id_paths:
			ids = self.__L_ID(id_paths)
			self.sample_df = full_df[full_df["Library number"].isin(ids)].reset_index(drop=True)
		else:
			self.sample_df = full_df
		self.factor_key = {"pheno": [], "key": [], "val": []}


	def __str__(self):
		return repr(self.sample_df)


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


	def count_pairs(self, pair_N):
		pairs = self.get_pairs()
		return len(pairs[pairs == pair_N])


	def extract_year(self, col, replace=True):
		dates = self.sample_df[col]
		years = dates[~dates.isna()].astype(str)
		if (years.str.find(".") != -1).any():
			years = years.str[:4]
		else:
			years = years.str[-4:]
		nans = dates[dates.isna()]
		years = pd.concat([years, nans]).sort_index()
		if not replace:
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
		subset.to_csv(path, index=False, sep="\t")


	def fast_stats(self):
		print(f"Samples: {self.sample_df.shape[0]}")
		print(f"Hosts: {self.sample_df[self.sample_df['Tissue'] == 'Host'].shape[0]}")
		print(f"Tumors: {self.sample_df[self.sample_df['Tissue'] == 'Tumour'].shape[0]}")
		print(f"Missing microchips: {self.nan_sum('Microchip', print_flag=False)}")
		print(f"Singles: {self.count_pairs(1)}")
		print(f"Pairs: {self.count_pairs(2)}")
		print(f"Triplets: {self.count_pairs(3)}")


	# col can be a string or list of strings for clustering on multiple columns
	def get_groups(self, col):
		return self.sample_df.groupby(col)


	def get_pairs(self, pair_count=None):
		pairs = self.sample_df.groupby("Microchip")["Microchip"].count()
		if pair_count:
			pairs = pairs[pairs == pair_count]
		return pairs


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
			print(f"{col}: {nans}/{total} ({nans/total*100:.2f}%) NaNs")
		else:
			return nans


	# TEST SOME OF THE NORMED VALUES BY HAND TO ENSURE THIS IS WORKING
	def normalize(self, cols, method="standard", inplace=True):
		methods = ["standard", "log"]
		cols = cols if isinstance(cols, list) else [cols]
		df_cols = self.sample_df[cols]
		if method == "standard":
			means = df_cols.mean()
			stds = df_cols.std()
			norms = (df_cols - means) / stds
		elif method == "log":
			norms = np.log(df_cols)
		else:
			print(f"*normalize() ERROR* Method \"{method}\" does not exist! Use: {', '.join(methods)}")
			return None
		if not inplace:
			return norms
		for col in cols:
			self.sample_df[col] = norms[col]



	def pair_rows(self, pair_count):
		pairs = self.get_pairs()
		pairs = pairs[pairs == pair_count]
		chips = pairs.index.values
		return self.sample_df[self.sample_df["Microchip"].isin(chips)]


	def plot_col(self, col, outpath, plot="hist"):
		plots = ["hist"]
		data = self.sample_df[col]
		if plot == "hist":
			q25, q75 = np.percentile(data, [0.25, 0.75])
			bin_width = 2 * ((q75 - q25) / len(data)**(-1/3))
			bins = round((data.max() - data.min()) / bin_width)
			plt.hist(data, bins=bins, density=True)
		else:
			print(f"*plot_col() ERROR* Plot \"{plot}\" does not exist! Use: {', '.join(plots)}")
			return None
		plt.ylabel(col)
		plt.savefig(outpath)
		plt.close()
		print(f"Plot saved to: {outpath}")


	def print_col(self, col):
		print(self.sample_df[col].to_list())


	def print_factor_key(self):
		if not self.factor_key["key"]:
			print("No factors have been created")
			return None
		factor_df = pd.DataFrame(self.factor_key)
		for pheno in factor_df["pheno"].unique():
			current = factor_df[factor_df["pheno"] == pheno]
			print(pheno)
			print(current[["key", "val"]].to_string(index=False))


	def subset(self, col, pattern):
		self.sample_df = self.sample_df[self.sample_df[col] == pattern]


	def subset_pairs(self):
		indices = self.sample_df.groupby("Microchip")["Microchip"].count()
		indices = indices[indices > 1].index
		self.sample_df = self.sample_df[self.sample_df["Microchip"].isin(indices)]


	def subset_non_nan(self, col):
		na = self.sample_df[col].isna()
		if isinstance(col, list):
			na = na.any(axis=1)
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


	def ATOMM_pheno(self, outpath, cols):
		self.subset_pairs()
		self.__handle_triplets()
		self.to_vcf_chip()
		tumor = self.chip_sort(self.sample_df[self.sample_df["Tissue"] == "Tumour"], num_only=True)
		self.subset("Tissue", "Host")
		host = self.chip_sort(self.sample_df, num_only=True)
		host.rename(columns={"Microchip": "host_chip"}, inplace=True)
		host["tumor_chip"] = tumor["Microchip"]
		self.sample_df = host
		self.to_pheno(outpath, cols, id_cols=["host_chip", "tumor_chip"], vcf_chip=False)


	def to_pheno(self, outpath, cols, id_cols=["Microchip"], vcf_chip=True):
		if vcf_chip:
			self.to_vcf_chip()
		pheno_df = self.sample_df[id_cols + cols]
		pheno_df = pheno_df.fillna("NA")
		outdir = outpath[:outpath.rfind(".")]
		factor_out = f"{outdir}_FACTOR_KEY.txt"
		self.factor_file(factor_out, cols=cols, silent=True)
		pheno_df.to_csv(outpath, index=False, sep="\t")


	def to_vcf_chip(self, inplace=True):
		chips = self.sample_df["Microchip"]
		chips = chips.str[-6:]
		tissue = self.sample_df["Tissue"].str[0]
		tissue[tissue == "T"] += self.sample_df["TumourNumber"]
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


	def __handle_triplets(self):
		triplets = self.get_pairs(3)
		if triplets.empty:
			return None
		print(f"Found {len(triplets)} triplets, please select only 1 of the duplicates to retain:")
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
		sample_csv_path="/work_bgfs/d/dgallinson/data/pheno_data/master_corrections.csv",
		tumor_csv_path="/work_bgfs/d/dgallinson/data/pheno_data/originals/RTumourTableRodrigo.csv",
		id_paths=[],
		tissue="Host"
		):
		super().__init__(sample_csv_path, id_paths)
		tumor_df_full = pd.read_csv(tumor_csv_path, dtype={"TumourNumber": str})
		if tissue != "both":
			self.sample_df = self.sample_df[self.sample_df["Tissue"] == tissue]
		microchips = self.sample_df["Microchip"].unique()
		self.tumor_df = tumor_df_full[tumor_df_full["Microchip"].isin(microchips)]
		self.over_mmax = {"Microchip": [], "volume": []}
		self.under_40 = []

	
	# ANALYSIS PHENOTYPE
	def estimate_age(self):
		tmp_tumor_df = self.tumor_df[["Microchip", "TrapDate", "TumourDepth", "TumourLength", "TumourWidth"]].copy()
		tmp_tumor_df["TrapDate"] = pd.to_datetime(tmp_tumor_df["TrapDate"])
		dates_df = self.__init_tumor_date(tmp_tumor_df)
		YOB = self.sample_df.set_index("Microchip")["YOB"]
		YOB = YOB[~YOB.index.duplicated(keep="first")]
		YOB = pd.to_datetime(YOB)
		merged = pd.concat([dates_df, YOB[dates_df.index]], axis=1)
		infection_age = (merged["init_tumor_date"] - merged["YOB"]).dt.days
		infection_age.name = "infection_age"
		infection_age = infection_age.reset_index()
		self.sample_df = self.sample_df.merge(infection_age, how="left", on="Microchip")


	def fast_stats_tumor(self):
		unique_samples = len(self.sample_df["Microchip"].unique())
		found_in_tumor = len(self.tumor_df["Microchip"].unique())
		chip_date_cluster = self.tumor_df.groupby(["Microchip", "TrapDate"])["TrapDate"].count()
		trap_cluster = chip_date_cluster.groupby("Microchip").count()

		print(f"Found in tumor DB: {(found_in_tumor/unique_samples)*100:.1f}% ({found_in_tumor}/{unique_samples})")
		print(f"1 trap: {len(trap_cluster[trap_cluster == 1])}")
		print(f"More than 1 trap: {len(trap_cluster[trap_cluster > 1])}")


	# ANALYSIS PHENOTYPE
	# Survival proxy following Margres et at. (2018; DOI: 10.1111/mec.14853) and using a growth back calculation derived from Wells et al. (2017; DOI: 10.1111/ele.12776)
	# OVERALL: obtain the difference in days between the first and last trapping trip and throw out those <40 days (including single-trip events). Calculate volume
	# via width x depth x length and put this into the back-calculation to get days since the tumor was 3 mm^3.
	# DETAILED
	# 1) Obtain the difference in days from first to last trapping event
	# 2) Throw out days < 40 (single trap events get tossed) and any rows with NA in width, depth, or length (also save the Microchips of tossed out samples)
	# 3) Group on Microchip and get the min date per Microchip group. If the devil had multiple tumors (most did), then there will be multiple min dates
	# 	3a) Within a Microchip group, get the tumorDB row(s) corresponding to the min date
	# 		*NOTE: herein lies my most awful loop, surely an elegant solution exists to this
	# 4) Calculate tumor volume
	# 	4a) Group on Microchip and get a single number per group (i.e., max or mean)
	# 5) Growth back calculation (also save volumes >= mmax-1)
	# 6) Get date of tumor start (at volume of 3 mm^3) and devil survival in days
	# 	6a) Group on Microchip and get the min and max dates
	# 	6b) Convert the back calculated value to a positive value in days and subtract that from the start (min) date to obtain tumor start date
	# 	6c) Subtract tumor start date from last trap (max) date to obtain survival in days
	# 7) Collate all of the info from step 6 and add it to the relevant Microchips in sample_df (missing values get NaN)
	def survival_proxy(self, day_cutoff=40):
		tmp_tumor_df = self.tumor_df[["Microchip", "TrapDate", "TumourDepth", "TumourLength", "TumourWidth"]].copy()
		tmp_tumor_df["TrapDate"] = pd.to_datetime(tmp_tumor_df["TrapDate"])
		
		date_diff = tmp_tumor_df.groupby("Microchip")["TrapDate"].transform(lambda x: x.max() - x.min())
		tmp_tumor_df["date_diff"] = date_diff.dt.days
		self.under_40 = pd.Series(tmp_tumor_df[tmp_tumor_df["date_diff"] < day_cutoff]["Microchip"].unique())
		tmp_tumor_df = tmp_tumor_df[tmp_tumor_df["date_diff"] >= day_cutoff]
		
		dates_df = self.__init_tumor_date(tmp_tumor_df)
		last_trap = tmp_tumor_df.groupby("Microchip")["TrapDate"].max().loc[dates_df.index]
		devil_survival_days = last_trap - dates_df["init_tumor_date"]
		
		calc_df = pd.DataFrame({
			"Microchip": dates_df.index.values,
			"volume (cm^3)": np.round(dates_df["tumor_volume"].values, 2),
			"back_calc": dates_df["back_calc"].values,
			"first_trap": dates_df["min_date"].dt.date.values,
			"last_trap": last_trap.dt.date.values,
			"init_tumor_date": dates_df["init_tumor_date"].dt.date.values,
			"devil_survival_days": devil_survival_days.dt.days.values})
		self.sample_df = self.sample_df.merge(calc_df, how="left", on="Microchip")


	def get_multi_sub_40(self):
		if len(self.under_40) == 0:
			print("Please run survival_proxy() before using this method")
			return None
		single_trap = self.get_trap(1)
		return self.under_40[~self.under_40.isin(single_trap)]


	def get_trap(self, trap_N):
		trap_nums = self.tumor_df.groupby(["Microchip", "TrapDate"])["TrapDate"].count().groupby("Microchip").count()
		return trap_nums[trap_nums == trap_N].index


	# This fails to get the exact tumors we sequenced but gets close. Due to mismatching info
	# between the saple CSV and tumor CSV, further resolution is impossible. Once I obtain the
	# proper phenotype data I will update this method (likely to grab based on TrapDate)
	def get_sample_tumors(self):
		sample_subset = self.sample_df[self.sample_df["Microchip"].isin(self.tumor_df["Microchip"])]
		sample_tumors = sample_subset[["Microchip", "TumourNumber"]]
		found_indices = []
		for i in range(sample_tumors.shape[0]):
			current_sample = sample_tumors.iloc[i]
			found_index = self.tumor_df.loc[(self.tumor_df["Microchip"] == current_sample["Microchip"]) & (self.tumor_df["TumourNumber"] == current_sample["TumourNumber"])]
			found_index = list(found_index.index.values)
			found_indices += found_index
		return self.tumor_df.loc[found_indices]


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


	# ANALYSIS PHENOTYPE
	def tumor_count(self):
		num_tumors = self.tumor_df.groupby(["Microchip"])["TumourNumber"].unique().str.len()
		num_tumors = num_tumors.to_frame().reset_index()
		num_tumors = num_tumors.rename(columns={"TumourNumber": "tumor_count"})
		self.sample_df = self.sample_df.merge(num_tumors, how="left", on="Microchip")


	# Definition: 
	# 			  A logistic growth back calculation derived from Wells et al. (2017; DOI: 10.1111/ele.12776) which was based on data from WPP.
	# 			  Output is a negative number in days representing days from the min capture date since the tumor was 3 mm^3.
	# 			  The model is not accurate enough to estimate beyond day resolution, and thus rounding to the nearest day 
	# 			  (default round_flag) is advised. This also stores any volumes greater than mmax.
	# Arguments:
	# 			  tumor_volumes: a vector/Series of tumor volumes
	# 			  mmax: 		 max tumor size in cm^3 (from growth model)
	# 			  init_size: 	 initial tumor size (from growth model)
	# 			  alpha: 		 scale parameter of the logistic growth curve (from growth model)
	# 			  is_cm: 		 flag if the incoming volume is in centimeters
	# 			  round_flag: 	 flag if the back calculated volumes should be rounded to the nearest int
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
	def __growth_back_calculation(self, tumor_volumes, mmax=202, init_size=0.0003, alpha=0.03, is_cm=False, round_flag=True):
		tumor_volumes = tumor_volumes.copy()
		if not is_cm:
			tumor_volumes /= 1000
		over = tumor_volumes[(tumor_volumes + 1) >= mmax]
		if over.size > 0:
			vstring = "volumes" if len(over) > 1 else "volume"
			print(f"*WARNING* found {len(over)} {vstring} greater than mmax while performing the back calculation (try print_over_mmax() for details)")
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
	# 			  If this is being used with the growth back calculation, max should set for mode. Drops rows with an NA measurement.
	# Arguments:
	#			  tumor_df: a DF with at least the following cols: Microchip, TrapDate, TumourDepth, TumourLength, TumourWidth
	# 			  mode: 	method of selecting a single volume when multiple tumors are present at the same date
	# Returns: 
	# 			  a DF containing Microchip as the index and cols tumor_volume (in starting units) and min_date (representing the minimum date for a sample)
	# Steps:
	#			  1) get min date for a sample (group on Microchip)
	# 			  2) Iterate through each of these microchips, find all rows with the min date, and generate a DF from these
	# 			  3) Obrain tumor volumes via length x width x depth
	# 			  4) Group these volumes on Microchip and select a single volume for each microchip (i.e., max, mean, sum)
	# Gotchas:
	# 			  Drops rows if any NA is found in TumourDepth, TumourLength, or TumourWidth
	# TODO:
	# 			  The loop is an inelegant and inefficient solution, this should be vectorized
	def __min_cap_volume(self, tumor_df, mode="max"):
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
	# 			  This is  a convenience function combining the functionality of __min_cap_volume()
	# 			  and __growth_back_calculation() to find the initial tumor date.
	# Arguments:
	# 			  tumor_df: a DF with at least the following cols: Microchip, TrapDate, TumourDepth, TumourLength, TumourWidth
	# Returns:
	# 			  a DF with Microchip as the index and the following cols: tumor_volume, min_date, back_calc, init_tumor_date
	# Steps:
	# 			  1) Run __min_cap_volume() and __growth_back_calculation
	# 			  2) Get the initial tumor date by subtracting the back calc from the min date
	# 			  3) Combine everything into a DF
	def __init_tumor_date(self, tumor_df):
		volumes = self.__min_cap_volume(tumor_df)
		volumes["tumor_volume"] /= 1000
		back_calc = self.__growth_back_calculation(volumes["tumor_volume"], is_cm=True)
		init_tumor_dates = volumes["min_date"] - pd.to_timedelta(-1 * back_calc, unit="d")
		init_tumor_dates.name = "init_tumor_date"
		back_calc = pd.Series(back_calc, name="back_calc", index=volumes.index)
		return pd.concat([volumes, back_calc, init_tumor_dates], axis=1)



class Compare:
	def __init__(self, df, cols):
		self.df = df
		self.cols = cols
		self.matched = pd.DataFrame([])
		self.unmatched = pd.DataFrame([])


	def compare_strict(self):
		self.matched = self.df[self.df[self.cols[0]] == self.df[self.cols[1]]]
		self.unmatched = self.df[self.df[self.cols[0]] != self.df[self.cols[1]]]


	def diff(self, n=1, biggest=True):
		first = self.unmatched[self.cols[0]]
		second = self.unmatched[self.cols[1]]
		abs_diffs = sorted((first - second).abs(), reverse=biggest)
		return abs_diffs[0:n]


	def diff_stats(self):
		first = self.unmatched[self.cols[0]]
		second = self.unmatched[self.cols[1]]
		abs_diffs = (first - second).abs()
		mean = abs_diffs.mean()
		median = abs_diffs.median()
		return mean, median


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