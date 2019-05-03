# This should take about 1-2 hours to run, and require about 20 GB of
# memory.
library(data.table)

# SCRIPT PARAMETERS
# -----------------
pheno.file <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                        "12-feb-2019","ukb26140.csv.gz")
out.file   <- "height.csv"
cols       <- c("eid","31-0.0","50-0.0")
col_names  <- c("id","sex","height")

# LOAD DATA
# ---------
out <- system.time(
  dat <- fread("ukb26140.csv.gz",sep = ",",header = TRUE,verbose = TRUE,
               showProgress = TRUE,colClasses = "character",nrows = 10000))
class(dat) <- "data.frame"

# PREPARE DATA
# ------------
# Select the reequested columns.
dat        <- dat[cols]
names(dat) <- col_names

# Convert all columns except the first one to numeric values, and set
# all empty strings to NA.
n <- length(dat)
for (i in 2:n) {
  x          <- dat[,i]
  x[x == ""] <- as.character(NA)
  dat[,i]    <- as.numeric(x)
}

# Remove all rows in which one or more of the values are missing.
rows <- which(rowSums(is.na(dat)) == 0)
dat  <- dat[rows,]

# SUMMARIZE DATA
# --------------
# This is to double-check that everything looks okay.
summary(dat)

# WRITE DATA TO FILE
# ------------------
write.csv(dat,out.file,row.names = FALSE,quote = FALSE)
