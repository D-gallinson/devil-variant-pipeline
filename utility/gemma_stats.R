library(ggplot2)
library(rstan)
library(stringr)
library(reshape2)

downsample <- function(df, samples=-1){
  N <- nrow(df)
  if(samples == -1){
    if(N < 1e5){
      samples <- N
    } else{
      samples <- 5e4
    }
  }
  iter <- round(nrow(df) / samples)
  return(df[seq(1, N, iter),])
}

logit <- function(x){
  return(log(x / (1 - x)))
}

to_sim_mat <- function(df, chains=1){
  arr <- array(unlist(df), c(nrow(df), chains, ncol(df)),
               dimnames=list(NULL, NULL, colnames(df)))
  return(arr)
}

to_CI <- function(x){
  return(quantile(x, c(0.025, 0.975)))
}

violin <- function(df, color, ylab, ylim=c(0, NA)){
  ggplot(df, aes(x=variable, y=value)) +
    geom_violin(fill=color) +
    stat_summary(fun=mean, geom="point", size=2, color="black") +
    stat_summary(fun=to_CI, geom="line", size=0.25, color="black") +
    ylim(ylim) +
    ylab(ylab) +
    theme(axis.title.x=element_blank())
}

# Read in input/output CLI args
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]

# Adjust based on the GEMMA run
iters <- 10
transform <- T

# Prepare the output fname
if(substr(output, nchar(output), nchar(output)) == "/"){
  output <- substr(output, 1, nchar(output)-1)
}
out_dir <- str_extract(output, regex("([^/]+$)"))
output <- paste(output, out_dir, sep="/")

# Read in prefix.hyp.txt and downsample iterations
gemma_hyp_full <- read.csv(input, sep="\t")
model_params <- colnames(gemma_hyp_full)
gemma_hyp_full$iter <- seq(1, nrow(gemma_hyp_full)*iters, iters)
stats_df <- gemma_hyp_full[, c("pve", "pge", "n_gamma")]

logit_cols <- c("h", "pve", "rho", "pge")
log_cols <- c("pi", "n_gamma")
if(transform){
  gemma_hyp_full[, logit_cols] <- logit(gemma_hyp_full[, logit_cols])
  gemma_hyp_full[, log_cols] <- log(gemma_hyp_full[, log_cols])
  output <- paste(output, "transform", sep=".")
}

gemma_hyp <- downsample(gemma_hyp_full)

# Generate DF to hold summary stats
stat_len <- 4
stat_df <- data.frame(pve = numeric(stat_len),
                      pge = numeric(stat_len),
                      n_gamma = numeric(stat_len),
                      row.names = c("Mean", "2.5%", "97.5%", "Median"))

# rstan posterior distribution stats
sim_mat <- to_sim_mat(gemma_hyp_full[,1:length(gemma_hyp_full)-1])
fit_stats <- monitor(sim_mat, warmup=0)

# Generate plots and stats per parameter
pdf(paste(output, "convergence_plots.pdf", sep="."))
for(param in model_params){
  y_label <- param
  if(param %in% logit_cols && transform){
    y_label <- paste("logit(", param, ")", sep="")
  }else if(param %in% log_cols && transform){
    y_label <- paste("log(", param, ")", sep="")
  }
  
  plot <- ggplot(gemma_hyp, aes_string(x="iter", y=param)) + 
    geom_point() +
    ggtitle(paste("GEMMA Hyperparameter:", param)) + labs(x="Iter", y=y_label)
  print(plot)
  
  if(param %in% c("pve", "pge", "n_gamma")){
    col <- stats_df[, param]
    stat_df["Mean", param] <- mean(col)
    CI <- quantile(col, c(0.025, 0.975))
    stat_df["2.5%", param] <- CI["2.5%"]
    stat_df["97.5%", param] <- CI["97.5%"]
    stat_df["Median", param] <- median(col)
  }
}
dev.off()
write.csv(stat_df, paste(output, "stats.csv", sep="."))
write.csv(fit_stats, paste(output, "convergence_stats.csv", sep="."))

# Generate PVE and PGE violin plots
violin_df <- read.csv(input, sep="\t")
violin_df <- violin_df[, c("pve", "pge")]
colnames(violin_df) <- c("PVE", "PGE")
violin_df <- melt(violin_df)

plot <- violin(violin_df, "skyblue3", "Proportion Variance Explained")
ggsave(
  paste(output, "violin_plots.png", sep="."),
  plot=plot
)
