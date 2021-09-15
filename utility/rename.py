#Rename all files in a directory (or recursively through subdirectories) based on a CSV specifying an ID and rename column

import pandas as pd
import argparse as ag
import os
import re


class Rename:
	def __init__(self, root, name_path, id_col_name, rename_col_name, tissue_col_name, tumor_id_col_name, full_names, interactive_mode, undo_flag):
		self.interactive_mode = interactive_mode
		rename_df = pd.read_csv(name_path, dtype=str)
		
		if not undo_flag:
			if full_names:
				renames = rename_df[rename_col_name]
			else:
				renames = self.__gen_names(rename_df, rename_col_name, tissue_col_name, tumor_id_col_name)
			tIDs = rename_df[tumor_id_col_name]
			self.rename_list = self.__gen_rename_tID(renames, tIDs)
		else:
			self.rename_list = rename_df["original"]
			id_col_name = "rename"
		self.target_list = rename_df[id_col_name]


	def rename(self, path):
		rename_list = self.__rename_list(self.target_list, self.rename_list, path)
		self.__rename_files(rename_list, path)


	def rename_recur(self, path):
		rename_list = self.__rename_list(self.target_list, self.rename_list, path)
		for root, dirs, files in os.walk(path):
			if root[-1] != "/":
				root += "/"
			for dir in dirs:
				rename_list += self.__rename_list(self.target_list, self.rename_list, root + dir)
		self.__rename_files(rename_list, path)


	def undo(self, csv_path):
		name_df = self.read_csv(csv_path)
		targets = name_df["rename"]
		rename = name_df["original"]


	def __gen_names(self, df, rename_col, tissue_col, tumor_id_col):
		chips = df[rename_col].str[-6:]
		ids = df[tissue_col].str[0]
		names = ids + "-" + chips
		return names


	def __gen_rename_tID(self, renames, tIDs):
		dict = {"ID": renames.to_list(), "tID": tIDs.to_list()}
		df = pd.DataFrame(dict)
		tumors = df[df["ID"].str[0] == "T"].index.values
		df.loc[tumors, "ID"] = df.loc[tumors]["ID"].str[0] + df.loc[tumors]["tID"] + df.loc[tumors]["ID"].str[1:]
		return df["ID"]


	def __rename_list(self, targets, renames, path):
		if path[-1] != "/":
			path += "/"
		files = sorted([file for file in os.listdir(path) if not os.path.isdir(path + file)])
		rename_list = []
		for file in files:
			for i in range(len(targets)):
				target = targets[i]
				rename = renames[i]
				if target in file:
					src = path + file
					dest = file.replace(target, rename)
					dest = path + dest
					rename_list.append([src, dest])
					break
		return rename_list


	def __rename_files(self, rename_list, path):
		if not rename_list:
			print("Found 0 files to rename, exiting.")
			exit()

		if self.interactive_mode:
			display_names = input(f"Found {len(rename_list)} files to rename. Show renaming list? [(Y)es, (N)o, (H)ead] ").lower()
			if display_names == "y" or display_names == "h":
				head = True if display_names == "h" else False
				self.__print_flist(rename_list, head)
				rename_flag = True if input("Continue with renaming? [y/n] ").lower() == "y" else False
		else:
			print("********************************************************************************")
			self.__gen_name_key(rename_list, path)
			print("********************************************************************************\n")
			print(f"Found {len(rename_list)} files to rename.")
			self.__print_flist(rename_list, False)
			rename_flag = True

		if rename_flag:
			for pair in rename_list:
				os.rename(pair[0], pair[1])
			print(f"Renamed {len(rename_list)} files")

		else:
			print("No files renamed")


	def __gen_name_key(self, rename_list, path):
		name_list = []
		for name in rename_list:
			original = name[0][name[0].rfind("/")+1:]
			rename = name[1][name[1].rfind("/")+1:]
			name_list.append([original, rename])
		name_df = pd.DataFrame(name_list, columns=["original", "rename"])
		name_df.to_csv(f"{path}/rename_key.csv", index=False)
		print(f"Name key written to: {path}/rename_key.csv")
		print(f"---Run the following command to undo renaming---")
		print(f"python3 rename.py --undo {path} {path}/rename_key.csv")


	def __print_flist(self, rename_list, head):
		end = 5 if head else len(rename_list)
		for i in range(end):
			pair = rename_list[i]
			print(f"{pair[0]} -------> {pair[1]}")



parser = ag.ArgumentParser(description="Rename files specified by a CSV")
parser.add_argument("path", help="Directory path for files to be renamed")
parser.add_argument("rename_csv", help="Path to the CSV to guide renaming")
parser.add_argument("-ic", "--id-col", dest="id_col_name", default="Library number", help="Name of the ID column in the CSV file [default=Library number]")
parser.add_argument("-rc", "--rename-col", dest="rename_col_name", default="Microchip", help="Name of the rename column in the CSV file [default=Microchip]")
parser.add_argument("-tc", "--tissue-col", dest="tissue_col_name", default="Tissue", help="Name of the tissue column in the CSV file [default=Tissue]")
parser.add_argument("-tid", "--tid-col", dest="tumor_id_col_name", default="TumourNumber", help="Name of the tumor ID column in the CSV file [default=TumourNumber]")
parser.add_argument("-f", "--full-names", dest="full_names_flag", action="store_true", help="Used when the names are properly formatted: T#ID_Microchip[-6:]")
parser.add_argument("-i", "--interactive", dest="interactive_flag", action="store_true", help="Use the tool in interactive mode")
parser.add_argument("-r", "--recursive", dest="recur_flag", action="store_true", help="Flag to search subdirectories recursively")
parser.add_argument("-u", "--undo", dest="undo_flag", action="store_true", help="Undo renaming. Supply the parent directory containing all renamed files and the path to the rename_key.csv")
args = parser.parse_args()

rObj = Rename(
	args.path,
	args.rename_csv,
	args.id_col_name,
	args.rename_col_name,
	args.tissue_col_name,
	args.tumor_id_col_name,
	args.full_names_flag,
	args.interactive_flag,
	args.undo_flag
)


if args.recur_flag:
	rObj.rename_recur(args.path)
else:
	rObj.rename(args.path)