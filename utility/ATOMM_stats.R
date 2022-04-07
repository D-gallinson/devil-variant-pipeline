library(ggplot2)
library(ggsci)
library(reshape2)
library(stringr)
library(dplyr)
library(cowplot)
library(grid)
library(gridExtra)

# Single command line argument which points to the output folder
args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
diff_cutoff <- 1

# Sort the README on run (this is out of order due to the array job)
readpath <- paste(dir, "/README.txt", sep="")
readme <- read.csv(readpath, sep="\t")
readme <- readme[order(readme$Run, decreasing=F),]
write.table(readme, readpath, sep="\t", row.names=F, quote=F)

# Grab paths to all estimate.txt files
estimate_paths <- list.files(dir, pattern="estimate.txt", full.names=T, recursive=T)
dfs <- list()

# Read all estimate.txt files into a list of DFs
i <- 1
for(file in estimate_paths){
  run <- str_replace(basename(dirname(file)), "run", "")
  dfs[[i]] <- read.table(file, header=F)
  dfs[[i]]$Run <- as.numeric(run)
  i <- i + 1
}

# Combine the estimate DFs into a plottable DF and write the estimate DF
estimates <- bind_rows(dfs)
colnames(estimates) <- c("Host", "Pathogen", "Interaction", "Noise", "Run")
estimate_plot <- reshape2::melt(estimates[,1:4])
write.table(estimates, paste(dir, "/run_estimates.txt", sep=""), sep="\t", quote=F, row.names=F)

# Generate a boxplot for host, pathogen, interaction, and noise
ATOMM_box <- ggplot(estimate_plot, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  scale_fill_jama() +
  ylab("Proportion") + ggtitle("ATOMM Host-Pathogen Heritability") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        legend.position="none",
        text=element_text(size=22),
        plot.title = element_text(size = 20))

# Save the plot
ggsave(paste(dir, "/heritability_boxes.png", sep=""),
       ATOMM_box,
       height = 7.28,
       width = 7.48)


# Obtain the model fitting stats (function init value and end after fitting)
run_out <- list.files(dir, pattern="*.out", full.names=T, recursive=T)
dfs <- list()

# Read in each .out file, grabbing run#, init vals, init func val, final func val, and diff
i <- 1
for(file in run_out){
  func <- readLines(file, n=1000)
  run <- as.numeric(str_replace(func[2], "Run ", ""))
  init <- substr(func[3], 6, nchar(func[3]))
  start <- as.numeric(strsplit(str_squish(func[18]), " ")[[1]][3])
  end <- as.numeric(strsplit(str_squish(func[length(func)-5]), " ")[[1]][3])
  diff <- round(start - end, 3)
  dfs[[i]] <- data.frame(Run=run, Init=init, Start=start, End=end, Diff=diff)
  i <- i + 1
}

# Combine DFs into a single DF and write to output dir
model_fit <- bind_rows(dfs)
write.table(model_fit, paste(dir, "/model_fits.txt", sep=""), sep="\t", row.names=F, quote=F)

# Plot a histogram of the func diffs
func_plot <- ggplot(model_fit, aes(x=Diff)) +
  geom_histogram(color="black", fill="skyblue3") +
  xlab("Start - Final") +
  ggtitle("Difference From Model Initialization to Final Fit")

ggsave(paste(dir, "/fit_histogram.png", sep=""),
       func_plot)


# Merge the estimate and model DFs, expanding the init values
joint <- merge(estimates, model_fit, on="Run")
joint$Init <- gsub("\\[|\\]", "", joint$Init)
joint$Init <- str_replace_all(joint$Init, " ", "")
inits <- data.frame(str_split(joint$Init, ",", simplify=T))
inits <- mutate_all(inits, function(x) as.numeric(as.character(x)))
colnames(inits) <- c("Host_init", "Patho_init", "Inter_init")
inits$Noise_init <- 1 - rowSums(inits)
joint <- cbind(joint, inits)
joint_out <- joint[, !(names(joint) %in% "Init")]
write.table(joint_out, paste(dir, "/all_stats.csv", sep=""), sep=",", row.names=F, quote=F)

# Convenience plotting function
scatters <- function(df, x, y){
  plot <- ggplot(df, aes_string(x=x, y=y, size="Diff")) +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    theme_bw()
  return(plot)
}

# Attempt (perhaps in vain) to determine if model fitting is sufficient
# by comparing the init val to the heritability estimate scaled by model fit
host <- scatters(joint, "Host_init", "Host")
patho <- scatters(joint, "Patho_init", "Pathogen")
inter <- scatters(joint, "Inter_init", "Interaction")
noise <- scatters(joint, "Noise_init", "Noise")

plot <- plot_grid(host,
                  patho,
                  inter,
                  noise,
                  ncol=2)

ggsave(paste(dir, "/fit_vs_estimate.png", sep=""),
       plot)


# Make a final heritability estimate boxplot, excluding model fits below diff_cutoff
final_boxes <- joint[joint$Diff >= diff_cutoff,
                     c("Host", "Pathogen", "Interaction", "Noise")]
N <- nrow(final_boxes)
final_boxes <- melt(final_boxes)

ATOMM_box <- ggplot(final_boxes, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  scale_fill_jama() +
  ylab("Proportion") + ggtitle(paste("ATOMM Host-Pathogen Heritability (N=", N, ")", sep="")) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        legend.position="none",
        text=element_text(size=22),
        plot.title = element_text(size = 20))

# Save the plot
ggsave(paste(dir, "/heritability_boxes_trimmed.png", sep=""),
       ATOMM_box,
       height = 7.28,
       width = 7.48)