library(ggplot2)
library(rstan)

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

to_sim_mat <- function(df, chains=1){
  arr <- array(unlist(df), c(nrow(df), chains, ncol(df)),
               dimnames=list(NULL, NULL, colnames(df)))
  return(arr)
}

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]

# Adjust based on the GEMMA run
iters <- 10

# Read in prefix.hyp.txt and downsample iterations
gemma_hyp_full <- read.csv(input, sep="\t")
model_params <- colnames(gemma_hyp_full)
gemma_hyp_full$iter <- seq(1, nrow(gemma_hyp_full)*iters, iters)
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
pdf(paste(output, "convergence_plots.pdf", sep=""))
for(param in model_params){
  plot <- ggplot(gemma_hyp, aes_string(x="iter", y=param)) + 
                   geom_point() +
                   ggtitle(paste("GEMMA Hyperparameter:", param)) + xlab("Iter") + ylab(param)
  print(plot)
  
  if(param %in% c("pve", "pge", "n_gamma")){
    col <- gemma_hyp_full[, param]
    stat_df["Mean", param] <- mean(col)
    CI <- quantile(col, c(0.025, 0.975))
    stat_df["2.5%", param] <- CI["2.5%"]
    stat_df["97.5%", param] <- CI["97.5%"]
    stat_df["Median", param] <- median(col)
  }
}
dev.off()
write.csv(stat_df, paste(output, "stats.csv", sep=""))
write.csv(fit_stats, paste(output, "fit_stats.csv", sep=""))