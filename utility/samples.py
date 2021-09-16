import os
import re
import time
import pandas as pd
import numpy as np



class Samples:
	def __init__(self, sample_csv_path, id_paths=[]):
		full_df = pd.read_csv(sample_csv_path)
		if id_paths:
			ids = self.__L_ID(id_paths)
			self.sample_df = full_df[full_df["Library number"].isin(ids)].reset_index(drop=True)
		else:
			self.sample_df = full_df
		self.factor_key = {"key": [], "val": []}


	def __str__(self):
		return repr(self.sample_df)


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


	def pair_rows(self, pair_count):
		pairs = self.get_pairs()
		pairs = pairs[pairs == pair_count]
		chips = pairs.index.values
		return self.sample_df[self.sample_df["Microchip"].isin(chips)]


	def print_col(self, col):
		print(self.sample_df[col].to_list())


	def print_factor_key(self):
		if not self.factor_key["key"]:
			print("No factors have been created")
			return None
		for i in range(len(self.factor_key["key"])):
			print(f"{self.factor_key['key'][i]}\t{self.factor_key['val'][i]}")


	def subset(self, col, pattern):
		self.sample_df = self.sample_df[self.sample_df[col] == pattern]


	def susbset_non_nan(self, col):
		self.sample_df = self.sample_df[~self.sample_df[col].isna()]


	def summarize_pairs(self):
		pairs = self.get_pairs()
		print("PAIR\tCOUNT")
		for pair in pd.unique(pairs):
			print(f"{pair}\t{len(pairs[pairs == pair])}")


	def to_factor(self, col):
		phenotype = self.sample_df[col]
		unique_vals = np.sort(phenotype.dropna().unique())
		phenotype = phenotype.fillna("NA")
		factors = [i for i in range(len(unique_vals))]
		phenotype = phenotype.replace(unique_vals, factors)
		self.factor_key["key"] += list(unique_vals)
		self.factor_key["val"] += factors
		return phenotype


	def to_pheno(self, outpath, cols, vcf_chip=True, extract_year=True, auto_factor=True):
		if extract_year:
			year_cols = ["YOB", "TrappingDate"]
			for year_col in year_cols:
				self.extract_year(year_col)
		microchips = self.sample_df["Microchip"]
		if vcf_chip:
			microchips = self.to_vcf_chip()
			microchips.name = "Microchip"
		output_df_list = [microchips]
		if not isinstance(cols, list):
			cols = [cols]
		gen_factor_flag = False
		for col in cols:
			numeric_flag = self.__isnum(self.sample_df[col].dropna().iloc[0])
			if not numeric_flag and auto_factor:
				new_col = self.to_factor(col)
				gen_factor_flag = True
			else:
				new_col = self.sample_df[col]
			output_df_list.append(new_col)
		pheno_df = pd.concat(output_df_list, axis=1, keys=[col.name for col in output_df_list])
		pheno_df = pheno_df.fillna("NA")
		if gen_factor_flag:
			outdir = outpath[:outpath.rfind(".")]
			factor_out = f"{outdir}_FACTOR_KEY.txt"
			factor_key_df = pd.DataFrame(self.factor_key)
			factor_key_df.to_csv(factor_out, index=False, sep="\t")
		pheno_df.to_csv(outpath, index=False, sep="\t")


	def to_vcf_chip(self):
		chips = self.sample_df["Microchip"]
		chips = chips.str[-6:]
		tissue = self.sample_df["Tissue"].str[0]
		tissue[tissue == "T"] += self.sample_df["TumourNumber"]
		chips = tissue + "-" + chips	
		return chips


	def unique(self, col):
		return self.sample_df[col].unique()


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
	def __init__(self, sample_csv_path, tumor_csv_path, id_paths=[], tissue="Host"):
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
		if "devil_survival_days" not in self.sample_df.columns:
			print("Please run survival_proxy() before using this method")
			return None



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
		# STEP 1
		tmp_tumor_df["TrapDate"] = pd.to_datetime(tmp_tumor_df["TrapDate"])
		date_diff = tmp_tumor_df.groupby("Microchip")["TrapDate"].transform(lambda x: x.max() - x.min())
		tmp_tumor_df["date_diff"] = date_diff.dt.days
		# STEP 2
		self.under_40 = pd.Series(tmp_tumor_df[tmp_tumor_df["date_diff"] < day_cutoff]["Microchip"].unique())
		tmp_tumor_df = tmp_tumor_df[tmp_tumor_df["date_diff"] >= day_cutoff]
		volume_cluster = self.__min_cap_volume(tmp_tumor_df)
		# STEP 5
		back_calculation = self.__growth_back_calculation(volume_cluster)
		# STEP 6
		chip_ids = volume_cluster.index
		# STEP 6a
		chip_groups = tmp_tumor_df.groupby("Microchip")["TrapDate"]
		first_trap = chip_groups.min().loc[chip_ids]
		last_trap = chip_groups.max().loc[chip_ids]
		# STEP 6b
		tumor_start_date = first_trap - pd.to_timedelta(-1 * back_calculation, unit="d")
		# STEP 6c
		devil_survival_days = last_trap - tumor_start_date
		# STEP 7
		calc_df = pd.DataFrame({
			"Microchip": chip_ids,
			"volume (cm^3)": np.round(volume_cluster.values, 2),
			"back_calc": back_calculation,
			"first_trap": first_trap.dt.date.values,
			"last_trap": last_trap.dt.date.values,
			"tumor_start_date": tumor_start_date.dt.date.values,
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


	# ANALYSIS PHENOTYPE
	def tumor_count(self):
		num_tumors = self.tumor_df.groupby(["Microchip"])["TumourNumber"].unique().str.len()
		num_tumors = num_tumors.to_frame().reset_index()
		num_tumors = num_tumors.rename(columns={"TumourNumber": "tumor_count"})
		self.sample_df = self.sample_df.merge(num_tumors, how="left", on="Microchip")


	# A logistic growth back calculation derived from Wells et al. (2017; DOI: 10.1111/ele.12776) which was based on data from WPP.
	# Output is a negative number in days representing days from the min capture date since the tumor was 3 mm^3.
	# The model is not accurate enough to estimate beyond day resolution, and thus rounding to the nearest day (default round_flag) is advised.
	# This also stores any volumes greater than mmax
	# Domain: (0, inf)
	# Range: (-inf, 0); as volume -> 201, days -> -inf
	def __growth_back_calculation(self, tumor_volumes, mmax=202, init_size=0.0003, alpha=0.03, is_cm=False, round_flag=True):
		if not is_cm:
			tumor_volumes /= 1000
		over = tumor_volumes[(tumor_volumes + 1) >= mmax]
		if self.over_mmax["Microchip"]:
			over = over.drop(self.over_mmax["Microchip"])
		if len(over > 0):
			self.over_mmax["Microchip"] += list(over.index.values)
			self.over_mmax["volume"] += list(over.values)
		back_calculation = np.log((mmax / (tumor_volumes.values + 1) - 1) / (mmax - init_size)) / alpha
		if round_flag:
			back_calculation = np.round(back_calculation)
		return back_calculation


	# Obtain the tumor volume Series for the minimum capture date with microchips as the index. When a host has multiple tumors, different grouping
	# functions can be used (e.g., max or sum). If this is being used with the growth back calculation, max should set for mode. Drops rows with an NA measurement.
	# This inelegant method obtains each microchip group's minimum date and then loops through each microchip number, slices rows on the microchip,
	# and extracts all rows equal to that microchip's minimum date. Something utilizing Pandas would be nicer.
	def __min_cap_volume(self, tumor_df, mode="max"):
		tumor_df = tumor_df.dropna(subset=["TumourDepth", "TumourLength", "TumourWidth"])
		min_dates = tumor_df.groupby("Microchip")["TrapDate"].min()
		min_dates = {"Microchip": min_dates.index, "TrapDate": min_dates.values}
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
		return volume_cluster



# pd.set_option('precision', 0)

samples = "/work_bgfs/d/dgallinson/data/pheno_data/master_corrections.csv"
tumors = "/work_bgfs/d/dgallinson/data/pheno_data/originals/RTumourTableRodrigo.csv"
base = "/shares_bgfs/margres_lab/Devils/BEE_Probe_Data"

batch_ids = [f"{base}/Capture1_6-11-21/rename_key.csv", f"{base}/Capture2_7-29-21/rename_key.csv", f"{base}/Capture3/rename_key.csv", f"{base}/Capture4/rename_key.csv", f"{base}/Capture5/NVS109A_Margres_CaptureSeq5_R1"]
b2_id = [f"{base}/Capture2_7-29-21/rename_key.csv"]

# b2 = Samples(samples)
# b2.subset("Tissue", "Host")
# print(b2.pair_rows(2))
# b2.subset("Tissue", "Host")
# trap_year = b2.extract_year("TrappingDate", replace=False).dropna()
# trap_year = trap_year.astype(int)
# print(len(trap_year[trap_year > 2018]))
# b2.to_pheno("../../phenotype_YOB.txt", ["YOB", "Site"], auto_factor=False)

tumor = TumorSamples(samples, tumors, tissue="Tumour")
# matches = tumor.get_sample_tumors()

# print(matches[matches["Microchip"] == "982000356426081"][["Microchip", "TumourNumber"]])
# print(matches[matches["TumourNumber"] != "1"][["Microchip", "TumourNumber"]])
# tumor.to_pheno("../../survival_proxy.txt", ["first_trap", "last_trap", "volume (cm^3)", "back_calc", "tumor_start_date", "devil_survival_days"], auto_factor=False)