library(argparse)


inv_logit <- function(x){
  return(exp(x) / (1 + exp(x)))
}

log_norm <- function(x){
  log_t <- log(x)
  return((x - mean(x, na.rm=T)) / sd(x, na.rm=T))
}

logit <- function(x){
  return(log(x / (1 - x)))
}

log_norm_undo <- function(x){
  x <- exp(x)
  return(exp((x * sd(x)) + mean(x)))
  return((x + mean(x, na.rm=T)) * sd(x, na.rm=T))
}


parser <- ArgumentParser()
parser$add_argument("pheno_file", help="Path to the input phenotype file (this file will be overwritten by this script)")
parser$add_argument("-t", "--transform", default="", help="Transform the input phenotype data with: inv_logit, logit, or log [default=None]")
parser$add_argument("-x", "--xval", action="store_true", default=F, help="Xval mode to randomly subset the phenotype file to test/training [default=False]")
parser$add_argument("-s", "--split", type="integer", default=20, help="Percentage to split into the testing set (e.g., -s 40 means 40% test, 60% training) [default=20]")
args <- parser$parse_args()

# These commands suck, try flags
# args <- commandArgs(trailingOnly=TRUE)
# pheno_input <- args[1]
# transform <- args[2]
# test_percent <- 20

df <- read.csv(args$pheno_file, sep="\t")
pheno <- df[,2]

if(args$xval == F){
  if(args$transform == "inv_logit"){
    pheno <- inv_logit(pheno)   
  } else if(args$transform == "logit"){
    pheno <- logit(pheno)
  } else if(args$transform == "log"){
    pheno <- log(pheno)
  } else{
    print("!!!! Error in transform function, please use: inv_logit, logit, or log !!!!")
    q()
  }
  write.csv(pheno, args$pheno_file)
}else{
  test_percent <- round(nrow(df) * (args$split/100))
  na_rows <- sample.int(nrow(df), test_percent)
  df[na_rows, ] <- "NA"
  write.csv(df, args$pheno_file)
}