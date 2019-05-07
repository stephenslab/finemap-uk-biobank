# This should take about 1-2 hours to run, and require about 20 GB of
# memory.
library(data.table)

# SCRIPT PARAMETERS
# -----------------
pheno.file <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                        "12-feb-2019","ukb26140.csv.gz")
out.file   <- "/gpfs/data/stephens-lab/finemap-uk-biobank/height.csv"
cols       <- c("eid","31-0.0","50-0.0","54-0.0","21000-0.0","21022-0.0",
                "22006-0.0","22001-0.0","22000-0.0","22005-0.0",
                paste0("22009-0.",1:40), paste0("22011-0.",0:4))
col_names  <- c("id","sex","height","assessment_centre","ethnic_self","age",
                "ethnic_genetic","sex_genetic","genotype_measurement_batch",
                "missingness",paste0("pc_genetic",1:40),
                paste0("relatedness_genetic",0:4))

# LOAD DATA
# ---------
out <- system.time(
  dat <- fread(pheno.file,sep = ",",header = TRUE,verbose = TRUE,
               showProgress = TRUE,colClasses = "character",nrows = 100))
class(dat) <- "data.frame"

# PREPARE DATA
# ------------
# Select the reequested columns.
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
# (except for the NAs in the "relatedness_genetic" columns).
#
# When the "genetic ethnic grouping" column is included, this removes
# any samples that are not marked as being "White British".
cols <- which(!grepl("relatedness_genetic",names(dat)))
rows <- which(rowSums(is.na(dat[,cols])) == 0)
dat  <- dat[rows,]

# SUMMARIZE DATA
# --------------
# This is to double-check that everything looks okay.
summary(dat)

# WRITE DATA TO FILE
# ------------------
write.csv(dat,out.file,row.names = FALSE,quote = FALSE)
