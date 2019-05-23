# TO DO: Briefly explain here what this script does, and how to use it.
#
# This should take about 1-2 hours to run, and require about 20 GB of
# memory.

# SCRIPT PARAMETERS
# -----------------
input.file  <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                         "12-feb-2019","ukb26140.csv.gz")
# output.file <- "/gpfs/data/stephens-lab/finemap-uk-biobank/height.csv"
output.file <- "/scratch/pcarbone/height.csv"

# The columns selected for subsequent analyses are as follows:
#
#   sex (31)
#   height (50)
#   UK Biobank assessment centre (54)
#   self-reported ethnic background (21000)
#   age (21022)
#   genetic ethnic grouping (22006)
#   genetic sex (22001)
#   genotype measurement batch (22000)
#   missingness (22005)
#   heterogeneity, PCA corrected (22004)
#   genetic PCs (22009-0.1 - 22009-0.40)
#   genetic relatedness pairing (22011)
#
cols       <- c("eid","31-0.0","50-0.0","54-0.0","21000-0.0","21022-0.0",
                "22006-0.0","22001-0.0","22000-0.0","22005-0.0","22004-0.0",
                paste0("22009-0.",1:40), paste0("22011-0.",0:4))
col_names  <- c("id","sex","height","assessment_centre","ethnic_self","age",
                "ethnic_genetic","sex_genetic","genotype_measurement_batch",
                "missingness","heterozygosity_pca",paste0("pc_genetic",1:40),
                paste0("relatedness_genetic",0:4))

# SET UP ENVIRONMENT
# ------------------
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

# LOAD DATA
# ---------
cat("Reading data from CSV file.\n")
out <- system.time(
  dat <- fread(input.file,sep = ",",header = TRUE,verbose = FALSE,
               showProgress = FALSE,colClasses = "character"))
class(dat) <- "data.frame"
cat(sprintf("Table imported with %d rows.\n",nrow(dat)))
    
# PREPARE DATA
# ------------
# Select the requested columns.
cat("Preparing data.\n")
dat        <- dat[,cols]
names(dat) <- col_names

# Convert all columns except the first one to numeric values, and set
# all empty strings to NA.
n <- length(dat)
for (i in 2:n) {
  x          <- dat[,i]
  x[x == ""] <- as.character(NA)
  dat[,i]    <- as.numeric(x)
}

# Remove all rows in which one or more of the values are missing
# (aside from the in the "relatedness_genetic" columns).
#
# When the "genetic ethnic grouping" column is included, this removes
# any samples that are not marked as being "White British".
cols <- which(!grepl("relatedness_genetic",names(dat)))
rows <- which(rowSums(is.na(dat[,cols])) == 0)
dat  <- dat[rows,]
cat(sprintf("After removing rows with NAs, %d rows remain.\n",nrow(dat)))

# SUMMARIZE DATA
# --------------
# Double-check that everything looks okay.
summary(dat)

# WRITE DATA TO FILE
# ------------------
cat("Writing prepared data to CSV file.\n")
write.csv(dat,output.file,row.names = FALSE,quote = FALSE)
