#############################################
# === DEFINITION ===
# Calculate the MAF for a haploid file of 0s
# and 1s (outputs statistics regarding MAF)
# 
# === USAGE ===
# Rscript hap_maf.R <input_path> <start_col>
# start_col: 1-indexed integer representing
#            the first col with SNP data
#############################################

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
path <- args[1]
start_col <- args[2]

# Read in the haploid matrix, considering SNPs from col start_col onwards
haploid_mat <- read.table(path, sep = " ")
haploid_mat <- haploid_mat[, start_col:ncol(haploid_mat)]
cat("Columns with SNP data:", ncol(haploid_mat), "\n")

# Obtain the MAF for each variant site, then output overall site-wide stats
samples <- ncol(haploid_mat)
sample_maf <- rowSums(haploid_mat)
sample_maf <- sample_maf / samples
sample_maf[sample_maf > 0.5] <- 1 - sample_maf[sample_maf > 0.5]
zeros <- length(sample_maf[sample_maf == 0])
print(summary(sample_maf))
cat("Sample size:", length(sample_maf), "\n")
cat("Number of 0 MAFs: ",
    zeros,
    " (", round((zeros / length(sample_maf)) * 100, 2), "%)",
    "\n",
    sep = "")