library(ggplot2)
library(dplyr)

sample_type <- "SNPs-host"
variants <- "/work_bgfs/d/dgallinson/outputs/intermediates/8_joint-variants"
stat_out <- "/work_bgfs/d/dgallinson/outputs/results/variants/stats"

input_file_info <- paste(variants, paste(sample_type, ".INFO", sep=""), sep="/")
sample <- paste(variants, "vcftools_metrics", sample_type, sep="/")

info_summary_out <- paste(stat_out, paste(sample_type, "-info.summary", sep=""), sep="/")
sample_summary_out <- paste(stat_out, paste(sample_type, "-sample.summary", sep=""), sep="/")

pdf("depth.pdf")
mean_depth_df <- read.csv(paste(sample, ".ldepth.mean", sep = ""), sep = "\t")
nrow(mean_depth_df)
dev.off()
q()

pdf(paste(stat_out, paste(sample_type, "-info.pdf", sep=""), sep="/"))

# -----INFO column plots-----
info_df <- read.csv(input_file_info, sep = "\t", na.strings = "?")

info_DP <- ggplot(info_df, aes(DP)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_DP + ggtitle("INFO DP") + theme_light() + xlim(0, 300)

info_FS <- ggplot(info_df, aes(FS)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_FS + ggtitle("INFO FS") + theme_light()

info_QD <- ggplot(info_df, aes(QD)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_QD + ggtitle("INFO QD") + theme_light()

info_SOR <- ggplot(info_df, aes(SOR)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_SOR + ggtitle("INFO SOR") + theme_light()

info_MQ <- ggplot(info_df, aes(MQ)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_MQ + ggtitle("INFO MQ") + theme_light()

info_MQRankSum <- ggplot(info_df, aes(MQRankSum)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_MQRankSum + ggtitle("INFO MQRankSum") + theme_light()

info_ReadPosRankSum <- ggplot(info_df, aes(ReadPosRankSum)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
info_ReadPosRankSum + ggtitle("INFO ReadPosRankSum") + theme_light()

info_summary <- summary(info_df[,5:ncol(info_df)])
capture.output(info_summary, file = info_summary_out)

dev.off()


pdf(paste(stat_out, paste(sample_type, "-sample.pdf", sep=""), sep="/"))

# -----Sample Level Plots-----
# --depth
depth_df <- read.csv(paste(sample, ".idepth", sep = ""), sep = "\t")
depth <- ggplot(depth_df, aes(MEAN_DEPTH)) + geom_histogram(fill = "blue", colour = "black", alpha = 0.3)
depth + ggtitle("Mean Depth Per Sample") + theme_light()
depth_summary <- summary(depth_df$MEAN_DEPTH)

# --site-mean-depth
mean_depth_df <- read.csv(paste(sample, ".ldepth.mean", sep = ""), sep = "\t")
site_depth <- ggplot(mean_depth_df, aes(MEAN_DEPTH)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
site_depth + ggtitle("Mean Depth Per Site") + theme_light() + xlim(0, 300)
site_depth_summary <- summary(mean_depth_df$MEAN_DEPTH)

# --freq2
freq_df <- read.csv(paste(sample, ".frq", sep = ""), sep = "\t", col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "A1", "A2"))
freq_df$maf <- freq_df %>% select(A1, A2) %>% apply(1, function(x) min(x))
freq <- ggplot(freq_df, aes(maf)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
freq + ggtitle("Minor Allele Frequency") + theme_light()
freq_summary <- summary(freq_df$maf)

# --site-quality
quality_df <- read.csv(paste(sample, ".lqual", sep = ""), sep = "\t")
quality <- ggplot(quality_df, aes(QUAL)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
quality + ggtitle("Site Quality (GQ)") + theme_light()
quality_summary <- summary(quality_df$QUAL)

# --missing-indv
missing_indv_df <- read.csv(paste(sample, ".imiss", sep = ""), sep = "\t")
missing_indv <- ggplot(missing_indv_df, aes(F_MISS)) + geom_histogram(fill = "blue", colour = "black", alpha = 0.3)
missing_indv + ggtitle("Frequency of Missing Samples Per Individual") + theme_light()
missing_indv_summary <- summary(missing_indv_df$F_MISS)

# --missing-site
missing_site_df <- read.csv(paste(sample, ".lmiss", sep = ""), sep = "\t")
missing_site <- ggplot(missing_site_df, aes(F_MISS)) + geom_density(fill = "blue", colour = "black", alpha = 0.3)
missing_site + ggtitle("Frequency of Missing Samples Per Site") + theme_light()
missing_site_summary <- summary(missing_site_df$F_MISS)

dev.off()

# capture.output("Mean Depth Per Sample (--depth)\n", file = sample_summary_out)
capture.output(depth_summary, file = sample_summary_out, append = T)

# capture.output("\n\nMean Depth Per site (--site-mean-depth)", file = sample_summary_out, append = T)
capture.output(site_depth_summary, file = sample_summary_out, append = T)

# capture.output("\n\nMinor Allele Frequency (--freq2)\n", file = sample_summary_out, append = T)
capture.output(freq_summary, file = sample_summary_out, append = T)

# capture.output("\n\Site Quality (GQ) (--site-quality)\n", file = sample_summary_out, append = T)
capture.output(quality_summary, file = sample_summary_out, append = T)

# capture.output("\n\Frequency of Missing Samples Per Individual (--missing-indv)\n", file = sample_summary_out, append = T)
capture.output(missing_indv_summary, file = sample_summary_out, append = T)

# capture.output("\n\nFrequency of Missing Samples Per Site (--missing-site)\n", file = sample_summary_out, append = T)
capture.output(missing_site_summary, file = sample_summary_out, append = T)


# pdf(paste(stat_out, "DENSITY_SNPs-host.pdf", sep="/"))
# hist(snp_df$DP[snp_df$DP < 300], main = "DP Histogram", xlab = "DP", breaks = seq(from=0, to=300, by=5))
# hist(snp_df$FS, main = "FS Histogram", xlab = "FS", breaks = seq(from=0, to=60, by=2))
# 
# density_cols <- c("QD", "SOR","MQ", "MQRankSum", "ReadPosRankSum")
# 
# for(col in density_cols){
#   d <- density(snp_df[[col]], na.rm = T)
#   plot(d, main = paste(col, "Probability Density Function", sep=" "), xlab = col)
# }
# dev.off()