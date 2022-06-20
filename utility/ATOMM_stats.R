library(stringr)
library(dplyr)
library(argparse)
library(readr)


# 3 command line arguments:
# arg[1] = points to the output folder
# arg[2] = processing mode ("prior" or "permutation")
# arg[3] = within a model .out file, number of lines from the final line with the last MLE
#          NOTE: open a single .out file and count up from the bottom to get this value
# args <- commandArgs(trailingOnly=TRUE)
# dir <- args[1]
# mode <- args[2]
# MLE_final_line <- as.numeric(args[3])

# CLI arguments
parser <- ArgumentParser()
parser$add_argument("run_dir", help="Path to the directory containing all model runs")
parser$add_argument("log_dir", help="Path to the directory containing all log files for a model")
parser$add_argument("-m", "--mode", default="prior", help="Processing mode, either \"prior\" or \"permutation\" [default=prior]")
parser$add_argument("-l", "--line", type="integer", default=7, help="Number of lines from the bottom of the .out file containing the final MLE (indexed at line 1 as start) [default=7]")
parser$add_argument("-r", "--readme", default="README.txt", help="Path to the README file [default=README.txt]")
parser$add_argument("-o", "--outpath", default="stats.csv", help="Output path for the stats CSV [default=stats.csv]")
args <- parser$parse_args()

if(args$mode == "prior"){
  # Sort the README on run (this is out of order due to the array job)
  readme <- read.csv(args$readme, sep="\t")
  readme <- readme[order(readme$Run, decreasing=F),]
  write.table(readme, args$readme, sep="\t", row.names=F, quote=F)
  
  
  # Read in each .out log file, grabbing run#, init vals, init model val, final model val
  i <- 1
  run_out <- list.files(args$log_dir, pattern="*.out", full.names=T, recursive=T)
  dfs <- list()
  
  for(file in run_out){
    model <- readLines(file, n=10000)
    run <- as.numeric(str_replace(model[2], "Run ", ""))
    inits <- substr(model[3], 1, nchar(model[3]))
    inits <- str_extract(inits, "\\[.*\\]")
    inits <- substr(inits, 2, nchar(inits)-1)
    inits <- strsplit(inits, split=", ")
    init_host <- as.double(inits[[1]][1])
    init_patho <- as.double(inits[[1]][2])
    init_inter <- as.double(inits[[1]][3])
    init_noise <- 1 - sum(init_host, init_patho, init_inter)
    start <- as.numeric(strsplit(str_squish(model[18]), " ")[[1]][3])
    end <- as.numeric(strsplit(str_squish(model[length(model) - args$line]), " ")[[1]][3])
    dfs[[i]] <- data.frame(run=run, 
                           host_i=init_host, pathogen_i=init_patho, 
                           interaction_i=init_inter, noise_i=init_noise,
                           init_MLE=start, final_MLE=end)
    i <- i + 1
  }
  
  # Combine DFs into a single DF and write to output dir
  model_df <- bind_rows(dfs)
  
  # Grab paths to all estimate.txt files
  estimate_paths <- list.files(args$run_dir, pattern="estimate.txt", full.names=T, recursive=T)
  dfs <- list()
  
  # Read all estimate.txt files into a list of DFs
  i <- 1
  for(file in estimate_paths){
    # run <- str_replace(basename(dirname(file)), "run", "")
    run <- parse_number(file)
    dfs[[i]] <- read.table(file, header=F)
    dfs[[i]]$Run <- as.numeric(run)
    i <- i + 1
  }
  
  # Combine the estimate DFs into a DF
  estimates <- bind_rows(dfs)
  colnames(estimates) <- c("host_PVE", "pathogen_PVE", "interaction_PVE", "noise_PVE", "run")
  
  
  # Merge the model Df and estimates DF
  final_df <- merge(model_df, estimates)
  final_df <- final_df[order(final_df$final_MLE, decreasing=F),]
  write.table(final_df, args$outpath, sep=",", row.names=F, quote=F)
} else {
  # GENERAL
  # Of the 5 inits per permutation, get the best run (i.e., min MLE), remove all other init runs,
  # and replace all SLURM .out files with the best init.out file 
  
  # Delete all original log files (these are not helpful in permutation mode)
  log_files <- list.files(args$log_dir, recursive=F, full.names=T)
  unlink(log_files)
  
  # Get the full path to all perm directories
  perm_dirs <- list.files(args$run_dir, pattern="perm*", full.names=T, include.dirs=T, recursive=F)
  failed_runs <- c()
  for(sub_dirs in perm_dirs){
    # Loop through all init subdirs within the current perm dir
    outfiles <- list.files(sub_dirs, pattern="*.out", full.names=T, recursive=T)
    best_MLE <- double(length(outfiles))
    init_files <- character(length(outfiles))
    i <- 1
    for(outfile in outfiles){
      # Get the final MLE for each of the init runs for a given perm
      init_files[i] <- basename(outfile)
      model <- readLines(outfile, n=10000)
      if(length(model) < 12){
        failed_runs <- c(failed_runs, outfile)
        next
      }
      best_MLE[i] <- as.numeric(strsplit(str_squish(model[length(model) - args$line]), " ")[[1]][3])
      i <- i + 1
    }
    # Get the first init# with the best (i.e., lowest) MLE
    best_MLE <- best_MLE[!is.na(best_MLE)]
    best_init_i <- which(best_MLE == min(best_MLE))[1]
    best_init <- tools::file_path_sans_ext(init_files[best_init_i])

    # Move all run files out of the best MLE init# dir and delete all init dirs
    to_delete <- list.files(sub_dirs, full.names=T, include.dirs=T, recursive=F)
    keep_path <- paste(sub_dirs, best_init, sep="/")
    keep_files <- list.files(keep_path, full.names=T)
    to <- paste(sub_dirs, basename(keep_files), sep="/")
    file.rename(keep_files, to)
    unlink(to_delete, recursive=T, force=T)
    
    # Move the best init#.out file to the logs dir
    log_file <- paste(sub_dirs, "/", best_init, ".out", sep="")
    new_log <- paste(args$log_dir, "/", basename(dirname(log_file)), ".out", sep="")
    file.rename(log_file, new_log)
  }

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Save runs which failed to a CSV
  # MAYBE CHANGE? MIGHT DO NEW PERM MODE WHERE THIS NO LONGER OCCURS
  # failed_runs <- data.frame(failed_path=failed_runs)
  # write.table(failed_runs, paste(dir, "failed_runs.csv", sep="/"), sep=",", quote=F, row.names=F)
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  # GENERAL
  # Get a DF of all relevant stats from the permutation tests
  estimate_paths <- list.files(args$run_dir, pattern="estimate.txt", full.names=T, recursive=T)
  dfs <- list()

  # Read all estimate.txt files into a list of DFs
  i <- 1
  for(file in estimate_paths){
    run <- str_replace(basename(dirname(file)), "perm", "")
    dfs[[i]] <- read.table(file, header=F)
    dfs[[i]]$Run <- as.numeric(run)
    i <- i + 1
  }

  estimates <- bind_rows(dfs)
  colnames(estimates) <- c("Host", "Pathogen", "Interaction", "Noise", "Run")
  estimates <- estimates[, c(5, 1:4)]

  # Obtain the model fitting stats (function init value and end after fitting)
  run_out <- list.files(args$log_dir, pattern="*.out", full.names=T, recursive=T)
  dfs <- list()

  # Read in each .out file, grabbing run#, init vals, init func val, final func val, and diff
  i <- 1
  for(file in run_out){
    func <- readLines(file, n=10000)
    run <- as.numeric(str_replace(func[2], "Run ", ""))
    end <- as.numeric(strsplit(str_squish(func[length(func)-5]), " ")[[1]][3])
    dfs[[i]] <- data.frame(Run=run, MLE=end)
    i <- i + 1
  }
  model_fit <- bind_rows(dfs)

  # Combine the DFs and sort by run
  df <- merge(estimates, model_fit)
  df <- df[order(df$Run),]
  write.table(df, args$outpath, sep=",", row.names=F, quote=F)
}