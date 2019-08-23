# Create a table containing the top 1,000 GeneATLAS height
# associations from chromosome 1, and an additional 1,000 "null" SNPs
# (SNPs that are not strongly associated with height), and save this
# table in a CSV file.

# SCRIPT PARAMETERS
# -----------------
# Select this many strongly associated SNPs, and this many SNPs that
# are not strongly associated ("null" SNPs).
n <- 1000

# The file containing the GeneATLAS association results.
geneatlas.file <-
  file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/summary",
            "gene-atlas/height/imputed/not-normalized",
            "imputed.allWhites.50-0.0.chr1.csv.gz")

# The table of association results for the 2n selected SNPs is saved
# to this CSV file.
out.file <- "../data/geneatlas-height.csv"

# SET UP ENVIRONMENT
# ------------------
# Load the packages used in the analysis below.
library(readr)
suppressMessages(library(dplyr))

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD ASSOCIATION RESULTS
# ------------------------
# Load the summary statistics from the CSV file.
cat("Loading GeneATLAS association results.\n")
geneatlas <- suppressMessages(read_delim(geneatlas.file,delim = " ",
                                         progress = FALSE))
colnames(geneatlas) <- c("SNP","ALLELE","iscores","NBETA","NSE","PV")
class(geneatlas)    <- "data.frame"

# Filter out SNPs that are not included in the dbSNP database.
geneatlas <- subset(geneatlas,substr(SNP,1,2) == "rs")

# SELECT ASSOCIATION RESULTS
# --------------------------
# Select the most strongly associated SNPs.
cat("Selecting SNPs.\n")
geneatlas.top <- top_n(geneatlas,-n,PV)

# Randomly select "null" SNPs.
rows           <- which(geneatlas$PV > 1e-6)
rows           <- sample(rows,n)
geneatlas.null <- geneatlas[rows,]

# Combine these into a single table.
geneatlas <- rbind(geneatlas.top,geneatlas.null)

# WRITE RESULTS TO FILE
# ---------------------
cat("Writing association results for selected SNPs to file.")
write.csv(geneatlas,out.file,row.names = FALSE,
          quote = FALSE)
