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

library(argparse)
library(tibble)



get_mac <- function(hap_mat){
  samples <- ncol(haploid_mat)
  midpoint <- samples / 2
  site_mac <- rowSums(haploid_mat)
  site_mac[site_mac > midpoint] <- samples - site_mac[site_mac > midpoint]
  return(site_mac)
}

get_maf <- function(hap_mat){
  samples <- ncol(haploid_mat)
  site_maf <- rowSums(haploid_mat)
  site_maf <- site_maf / samples
  site_maf[site_maf > 0.5] <- 1 - site_maf[site_maf > 0.5]
  return(site_maf)
}

maf_filter <- function(hap_mat, chr_vec, maf_vec, maf){
  remove_i <- which(maf_vec < maf)
  cat("Number of SNPs:", nrow(hap_mat), "\n")
  cat("Removing:", length(remove_i), "\n")
  hap_mat <- cbind(chr = chr_vec, hap_mat)
  hap_mat <- hap_mat[-remove_i, ]
  return(hap_mat)
}

mac_filter <- function(hap_mat, chr_vec, mac_vec, mac){
  remove_i <- which(mac_vec < mac)
  cat("Number of SNPs:", nrow(hap_mat), "\n")
  cat("Removing:", length(remove_i), "\n")
  hap_mat <- cbind(chr = chr_vec, hap_mat)
  hap_mat <- hap_mat[-remove_i, ]
  return(hap_mat)
}

write_output <- function(filtered_mat, outpath){
  ids <- seq(nrow(filtered_mat))
  output_mat <- add_column(filtered_mat, ids, .after = "chr")
  write.table(output_mat, outpath, col.names = F, row.names = F, quote = F)
}

# CLI arguments
parser <- ArgumentParser()
parser$add_argument("haploid_path", help="Path to the input haploid sequence file")
parser$add_argument("-c", "--col", type="integer", default=3, help="First column containing SNP information [default=3]")
parser$add_argument("-chr", "--chr_col", type="integer", default=1, help="Column number containing the chr IDs [default=1]")
parser$add_argument("--maf", type="double", default=-1, help="Filter the haploid file by removing any site with a minor allele frequency (MAF) lower than this value [default=-1]")
parser$add_argument("--mac", type="integer", default=-1, help="Filter the haploid file by removing any site with a minor allele count (MAC) lower than this value [default=-1]")
parser$add_argument("-o", "--outpath", default="hap_filter.txt", help="Output path if filtering is conducted [default=hap_filter.txt]")
parser$add_argument("-v", "--verbose", type="logical", default=T, help="Print the commands to console")
args <- parser$parse_args()

if(args$verbose){
  cat("Parameters:\n")
  cat("\t", args$haploid_path, 
      "\n\t--col ", args$col, "\n", 
      "\t--chr_col ", args$chr_col, "\n",
      sep = "")
  if(args$maf > 0){
    cat("\t--maf", args$maf, "\n")
  } else if (args$mac > 0){
    cat("\t--mac", args$mac, "\n")
  }
  if (args$maf > 0 || args$mac > 0){
    cat("\t--outpath", args$outpath, "\n")
  }
  cat("\n")
}

# Read in the haploid matrix, considering SNPs from col start_col onwards
haploid_mat <- read.table(args$haploid_path, sep = " ")
chr_vec <- haploid_mat[, args$chr_col]
haploid_mat <- haploid_mat[, args$col:ncol(haploid_mat)]
cat("Columns with SNP data:", ncol(haploid_mat), "\n")

# If --maf or --mac are specified, filter the sequence file before obtaining MAF stats
if(args$maf > 0){
  site_maf <- get_maf(haploid_mat)
  haploid_mat <- maf_filter(haploid_mat, chr_vec, site_maf, args$maf)
} else if(args$mac > 0){
  site_mac <- get_mac(haploid_mat)
  haploid_mat <- mac_filter(haploid_mat, chr_vec, site_mac, args$mac)
}

if(args$maf > 0 || args$mac > 0){
  write_output(haploid_mat, args$outpath)
  haploid_mat <- haploid_mat[,2:ncol(haploid_mat)]
  cat("\n")
}

# Get the MAF of each site and the MAF=0 sites
site_maf <- get_maf(haploid_mat)
zeros <- length(site_maf[site_maf == 0])
zero_indices <- which(site_maf == 0)

# Output overall site-wide stats
cat("===== MAF STATS =====\n")
print(summary(site_maf))
cat("Number of SNPs:", length(site_maf), "\n")
cat("Number of 0 MAFs: ",
    zeros,
    " (", round((zeros / length(site_maf)) * 100, 2), "%)",
    "\n",
    sep = "")
if(zeros != 0){
  cat("First 3 rows with MAF=0: ", paste(zero_indices[1:3], collapse = ","),
  "\n",
  sep = "")
}