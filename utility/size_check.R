#############################################
# USAGE:
# Rscript size_check.R <input_path> [unit]
# input_path: path to input folder
# unit: B, K, M, or G for unit of bytes
#       (default=M)
#############################################

args <- commandArgs(trailingOnly=TRUE)
path <- args[1]
pathlen <- nchar(path)

# Remove trailing fslash if present
if (substr(path, pathlen, pathlen) == "/") {
  path <- substr(path, 1, nchar(path)-1)
}

# Get file sizes in dir path
files <- list.files(path, full.names=T, recursive=T, pattern="*.fastq.gz")
sizes <- numeric(length(files))
i <- 1
for (f in files) {
  sizes[i] <- file.size(f)
  i <- i + 1
}

# Convert to specified units (default MB)
if (is.na(args[2]) == F) {
  units <- list("B"=1,
                "K"=1024,
                "M"=1024**2,
                "G"=1024**3)
  sizes <- sizes / as.numeric(units[args[2]])
} else {
  sizes <- sizes / 1024**2
}

summary(sizes)