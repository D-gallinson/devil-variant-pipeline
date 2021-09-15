library(dplyr)

############################################################
# This script takes as input the following paths in order:
# .INFO file
# .sample file
# .site file
# output/path
# Three CSV files equivalent to R's summary() are then 
# generated in output/path
############################################################



################### FUNCTIONS #######################
# Function to generate quantiles within all cols of a matrix
quantiles <- function(df, quant){
  quarts <- data.frame(matrix(nrow=1, ncol=ncol(df)))
  colnames(quarts) <- colnames(df)
  for (i in 1:ncol(df)){
    quarts[, i] <- quantile(df[, i], quant, na.rm=TRUE)
  }
  return(quarts)
}

# Convenience function to add the full path and unique tissue_folder ID to a file name
output_path <- function(fname, dir="tables", ext="csv"){
  outname <- paste(fname, output_id, ext, sep=".")
  if(fname == ""){
    outname <- substr(outname, 2, nchar(outname))
  }
  full_out <- paste(output, dir, outname, sep="/")
  return(full_out)
}

# All functions below are to handle rm vals in apply
na_mean <- function(x){
  return(mean(x, na.rm=TRUE))
}

na_min <- function(x){
  return(min(x, na.rm=TRUE))
}

na_max <- function(x){
  return(max(x, na.rm=TRUE)) 
}

na_median <- function(x){
  return(median(x, na.rm=TRUE))
}

# Main processing, this generates and saves the summary table to a .csv
table_summary <- function(df, out){
  minima <- data.frame(t(apply(df, 2, na_min)))
  rownames(minima) <- "Min."
  
  first_qu <- quantiles(df, 0.25)
  rownames(first_qu) <- "1st Qu."
  
  median <- data.frame(t(apply(df, 2, na_median)))
  rownames(median) <- "Median"
  
  mean <- data.frame(t(apply(df, 2, na_mean)))
  rownames(mean) <- "Mean"
  
  third_qu <- quantiles(df, 0.75)
  rownames(third_qu) <- "3rd Qu."
  
  maxima <- data.frame(t(apply(df, 2, na_max)))
  rownames(maxima) <- "Max."
  
  summary_stats <- rbind(minima, first_qu, median, mean, third_qu, maxima)
  write.csv(summary_stats, out)
}

################### VARIABLES ###################
# Arguments must be in the order: output_id, output
args <- commandArgs(trailingOnly=TRUE)
output_id <- args[1]
output <- args[2]

# File names
base_name <- paste(output, "table", output_id, sep="/")
info_file <- paste(base_name, "INFO", sep=".")
sample_file <- paste(base_name, "sample", sep=".")
site_file <- paste(base_name, "site", sep=".")

info_file <- output_path("", dir="data", ext="INFO")
sample_file <- output_path("", dir="data", ext="sample")
site_file <- output_path("", dir="data", ext="site")

################### LOADED DATAFRAMES ###################
# Load the .INFO file, skipping the CHROM, POS, REF, and ALT cols
info <- read.csv(info_file, sep="\t", na.strings="?")
info <- info[,5:ncol(info)]

# Load the .sample file, skipping the CHROM and POS cols
sample <- read.csv(sample_file, sep="\t", na.strings="?")
sample <- sample[,2:ncol(sample)]

# Load the .site file and generate a MAF col. Extract the relevant cols
site <- read.csv(site_file, sep="\t", na.strings="?")
site$MAF <- site %>% select(REF_FREQ, ALT_FREQ) %>% apply(1, function(x) min(x))
site$MAC <- site %>% select(REF_COUNT, ALT_COUNT) %>% apply(1, function(x) min(x))
site <- site %>% select(QUAL, MEAN_DEPTH, F_MISS, MAF, MAC)

################### SUMMARY TABLES ###################
table_summary(info, output_path("info_summary"))
table_summary(sample, output_path("sample_summary"))
table_summary(site, output_path("site_summary"))

################### BOTTOM 10 MAC ###################
low_ten <- site[site$MAC <= 10,]
low_ten <- as.data.frame(table(low_ten$MAC))
colnames(low_ten) <- c("MAC", "count")
n_samples <- dim(sample)[1]
low_ten$MAC <- as.numeric(as.character(low_ten$MAC))
low_ten$MAF <- signif(low_ten$MAC/n_samples, 2)
low_ten <- low_ten[,c(1,3,2)]
write.csv(low_ten, output_path("bottom_mac.csv"), row.names=F)

################### SNP DENSITY ###################
density_df <- read.table(output_path("", dir="data", ext="snpden"), header = T)
density_freq <- as.data.frame(table(density_df$SNP_COUNT))
colnames(density_freq) <- c("#SNPs", "count")
write.csv(density_freq, output_path("SNP_density"), row.names=F)

# ======================== SNP Density in bp ========================
# These functions generate basic summary stats of the average base pairs
# between SNPs. This assumes the VCF file is position sorted and finds
# the distance between subsequent SNPs on a per-chromosome basis.
# 4 files are generated from this:
# 1) per_chrom.bp.snpden = bp snpden and SNP count per chromosome
# 2) all_chroms.bp.snpden = quartiles/avg over all chroms
# 3) big_chroms.bp.snpden = only chroms 1-6 and X quartiles/avg
# 4) little_chroms.bp.snpden = chrom Y and all unplaced scaffolds quartiles/avg
distance <- function(pos_vec){
  previous_snps <- pos_vec[1:length(pos_vec)-1]
  next_snps <- pos_vec[2:length(pos_vec)]
  diff <- next_snps - previous_snps
  return(round(sum(diff)/length(diff)))
}


chrom_diff <- function(pos_df){
  chroms <- unique(pos_df$CHROM)
  chroms_N <- length(chroms)
  chrom_df <- data.frame(chrom = numeric(chroms_N),
                         mean_SNP_distance = numeric(chroms_N),
                         SNP_count = numeric(chroms_N))
  i <- 1
  for(chrom in chroms){
    chrom_SNPs <- pos_df[pos_df$CHROM == chrom, "POS"]
    avg_distance <- distance(chrom_SNPs)
    chrom_df[i, "chrom"] <- chrom
    chrom_df[i, "mean_SNP_distance"] <- avg_distance
    chrom_df[i, "SNP_count"] <- length(chrom_SNPs)
    i <- i + 1
  }
  
  return(chrom_df)
}


gen_stats <- function(pos_df){
  per_chrom <- chrom_diff(pos_df)
  big_chrom <- per_chrom[1:7, c("mean_SNP_distance", "SNP_count")]
  little_chrom <- per_chrom[8:nrow(per_chrom), c("mean_SNP_distance", "SNP_count")]
  write.csv(per_chrom, output_path("per_chrom.bp.snpden"), row.names=F)
  sink(output_path("big_chroms.bp", ext="snpden"))
  print(summary(big_chrom))
  sink(output_path("little_chroms.bp", ext="snpden"))
  print(summary(little_chrom))
  sink(output_path("all_chroms.bp", ext="snpden"))
  print(summary(per_chrom[, c("mean_SNP_distance", "SNP_count")]))
  sink()
}


pos_df <- read.table(paste(output, "tmp.pos", sep="/"), header=T)
gen_stats(pos_df)