# This should take about 1-2 hours to run, and require about 20 GB of
# memory.
library(data.table)

# SCRIPT PARAMETERS
# -----------------
pheno.file <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                        "12-feb-2019","ukb26140.csv.gz")
out.file   <- "relatedness.csv"
cols       <- c("eid",paste0("22011-0.",0:4))
col_names  <- c("id",paste0("relatedness_genetic",0:4))

# LOAD DATA
# ---------
out <- system.time(
  dat <- fread("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes/12-feb-2019/ukb26140.csv.gz",sep = ",",header = TRUE,verbose = TRUE,
               showProgress = TRUE,colClasses = "character"))
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
