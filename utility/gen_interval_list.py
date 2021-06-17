import pandas as pd


df = pd.read_table("data/intervals/mSarHar1.11_scaffolds.csv")
df = df.iloc[0:-1]

assigned_chr = list(df["Assigned-Molecule"][:8])
unassigned_chr = list(df["GenBank-Accn"][8:])
chr = assigned_chr + unassigned_chr

cols = ["chr", "start", "stop", "+", "name"]
interval_df = df[["Sequence-Length", "Sequence-Name"]]
tmp_dict = {
	"chr": chr,
	"start": ["1" for i in range(len(interval_df))],
	"plus": ["+" for i in range(len(interval_df))]
}

interval_df = pd.concat([interval_df, pd.DataFrame(tmp_dict)], axis=1)
# interval_df = interval_df[["chr", "start", "Sequence-Length", "plus", "Sequence-Name"]]
interval_df = interval_df[["chr", "start", "Sequence-Length"]]
interval_df.to_csv("data/intervals.list", index=False, header=False, sep="\t")