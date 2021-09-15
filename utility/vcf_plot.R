library(ggplot2)
library(dplyr)
library(karyoploteR)
library(VariantAnnotation)
library(gwasvcf)

#############################################################
# This script generates plots for the following VCF metrics:
# ---info---
# DP, FS, QD, SOR, MQ, MQRankSum, ReadPosRankSum
# 
# ---sample---
#  sites, mean depth, missingness frequency
# 
# ---site---
# mean depth, QUAL, MAF, missingness frequency
# 
# 
#############################################################


################### FUNCTIONS ###################
# Convenience naming function
output_path <- function(fname, dir="plots", ext="pdf"){
  outname <- paste(fname, output_id, ext, sep=".")
  if(fname == ""){
    outname <- substr(outname, 2, nchar(outname))
  }
  full_out <- paste(output, dir, outname, sep="/")
  return(full_out)
}

# Arguments must be in the order: output_id, output
args <- commandArgs(trailingOnly=TRUE)
output_id <- args[1]
output <- args[2]

base_name <- paste(output, output_id, sep="/")
info_file <- paste(base_name, "INFO", sep=".")
sample_file <- paste(base_name, "sample", sep=".")
site_file <- paste(base_name, "site", sep=".")

info_file <- output_path("", dir="data", ext="INFO")
sample_file <- output_path("", dir="data", ext="sample")
site_file <- output_path("", dir="data", ext="site")

info <- read.csv(info_file, sep="\t", na.strings="?")
sample <- read.csv(sample_file, sep="\t", na.strings="?")
site <- read.csv(site_file, sep="\t", na.strings="?")

# -----INFO column plots-----
pdf(output_path("info_plots"))

info_DP <- ggplot(info, aes(DP)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_DP + ggtitle("INFO DP") + theme_light() + xlim(0, 300)

info_FS <- ggplot(info, aes(FS)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_FS + ggtitle("INFO FS") + theme_light()

info_QD <- ggplot(info, aes(QD)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_QD + ggtitle("INFO QD") + theme_light()

info_SOR <- ggplot(info, aes(SOR)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_SOR + ggtitle("INFO SOR") + theme_light()

info_MQ <- ggplot(info, aes(MQ)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_MQ + ggtitle("INFO MQ") + theme_light()

info_MQRankSum <- ggplot(info, aes(MQRankSum)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_MQRankSum + ggtitle("INFO MQRankSum") + theme_light()

info_ReadPosRankSum <- ggplot(info, aes(ReadPosRankSum)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_ReadPosRankSum + ggtitle("INFO ReadPosRankSum") + theme_light()

dev.off()


# -----Sample-Level Plots-----
pdf(output_path("sample_plots"))

# N_sites
depth <- ggplot(sample, aes(N_SITES)) + geom_histogram(fill = "blue", colour = "black", alpha = 0.3)
depth + ggtitle("Number of Sites Per Sample") + theme_light()

# mean_depth
depth <- ggplot(sample, aes(MEAN_DEPTH)) + geom_histogram(fill = "blue", colour = "black", alpha = 0.3)
depth + ggtitle("Mean Depth Per Sample") + theme_light()

# frequency_miss
depth <- ggplot(sample, aes(F_MISS)) + geom_histogram(fill = "blue", colour = "black", alpha = 0.3)
depth + ggtitle("Missingness Frequency Per Sample (lower is better)") + theme_light()

dev.off()


# -----Site-Level Plots-----
pdf(output_path("site_plots"))

# mean_depth
site_depth <- ggplot(site, aes(MEAN_DEPTH)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
site_depth + ggtitle("Mean Depth (DP) Per Site") + theme_light() + xlim(0, 100)

# quality
quality <- ggplot(site, aes(QUAL)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
quality + ggtitle("Site Quality") + theme_light() + xlim(0, 300)

# MAF
site$maf <- site %>% dplyr::select(REF_FREQ, ALT_FREQ) %>% apply(1, function(x) min(x))
freq <- ggplot(site, aes(maf)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
freq + ggtitle("Minor Allele Frequency") + theme_light()

# MAC
site$mac <- site %>% dplyr::select(REF_COUNT, ALT_COUNT) %>% apply(1, function(x) min(x))
binwidth <- floor(dim(sample)[1] / 30)
count <- ggplot(site, aes(mac)) + geom_histogram(fill = "blue", colour = "black", alpha = 0.3, bindwidth = binwidth)
count + ggtitle(paste("Minor Allele Count (binwidth = ", binwidth, ")")) + theme_light()

# missingness
missing_site <- ggplot(site, aes(F_MISS)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
missing_site + ggtitle("Frequency of Missing Samples Per Site (F_MISS - lower is better)") + theme_light()

# -----SNP Density-----
pdf(output_path("SNP_density"))

# density_file <- paste(output, paste(output_id, "snpden", sep="."), sep="/")
density_file <- output_path("", dir="data", ext="snpden")
density_df <- read.table(density_file, header = T)
density <- ggplot(density_df, aes(SNP_COUNT)) + 
  geom_bar(fill = "blue", colour = "black", alpha = 0.3, bindwidth = 1, bins = 50)
density + ggtitle("SNP Density over 125kb intervals (binwidth = 1)") + theme_light() + xlim(0, 50)

dev.off()


# -----Karyotype Density-----
# pdf(paste(output, "karyotype_density.pdf", sep="/"))
pdf(output_path("karyotype_density"))

density <- read.table(paste(output, "tmp.grange.density", sep="/"), header=TRUE)
mSarHar <- toGRanges(paste(output, "tmp.chroms", sep="/"))
snps <- toGRanges(density)
kp <- plotKaryotype(genome = mSarHar)
kpPlotDensity(kp, data = snps)

dev.off()