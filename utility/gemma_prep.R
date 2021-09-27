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


args <- commandArgs(trailingOnly=TRUE)
pheno_input <- args[1]
transform <- args[2]

df <- read.csv(pheno_input, sep="\t")
pheno <- df[,2]

if(is.na(transform)){
  if(transform == "inv_logit"){
    pheno <- inv_logit(pheno)   
  } else if(transform == "logit"){
    pheno <- logit(pheno)
  } else if(transform == "log"){
    pheno <- log(pheno)
  } else{
    print("!!!! Error in transform function, please use: inv_logit, logit, or log !!!!")
  }
  write.csv(pheno, pheno_input)
}

print(pheno)
print(log_norm_transform(pheno))
print(log_norm_undo(log_norm_transform(pheno)))