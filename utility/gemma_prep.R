library(argparse)


parser <- ArgumentParser()
parser$add_argument("pheno_file", help="Path to the input phenotype file")
parser$add_argument("outpath", help="Path to output the test_training matrix file")
parser$add_argument("-s", "--split", type="integer", default=20, help="Percentage to split into the testing set (e.g., -s 40 means 40% test, 60% training) [default=20]")
args <- parser$parse_args()

# Read in the phenotype file and obtain the indices of
# all rows without an NA phenotype. Randomize these indices
df <- read.table(args$pheno_file, sep = "\t", skip = 1)
non_na <- rownames(df[!is.na(df$V2),])
non_na <- sample(non_na)
n <- length(non_na)

# Obtain the base number of test set samples per run (div)
# and the runs which will have div+1 samples (rem)
split <- 100 / args$split
div <- floor(n / split)
rem <- n %% split

# Generate each run's test set sample size. The runs which have
# div+1 samples are randomized
sizes <- rep(div, split)
longer <- sample(c(rep(1, rem), rep(0, split - rem)))
col_sizes <- sizes + longer

# Generate a list of vectors such that each vector represents
# the phenotype indices which will be converted to NA (i.e., the 
# test set) for a given run
start <- 1
i <- 1
test_cols <- list()
for(size in col_sizes){
  nums <- non_na[start:(start + size - 1)]
  if(size == div && var(col_sizes) != 0){
    nums <- c(nums, -1)
  }
  test_cols[[i]] <- nums
  start <- start + size
  i <- i + 1
}

# Convert the vectors to a df where columns represent the
# phenotype inidices to be NAed in a GEMMA run (e.g., 10
# columns means 10 GEMMA runs with a 10/90 test/training split
# whereby col_1 represents the indices which will be NAed in
# the first GEMMA run)
out_df <- as.data.frame(t(do.call(rbind, test_cols)))
out_df[out_df == -1] <- ""
write.table(out_df, args$outpath, col.names = F, sep = "\t", row.names = F, quote = F)