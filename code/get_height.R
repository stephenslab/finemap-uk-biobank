# This should take about 1-2 hours to run, and require about 20 GB of
# memory.
library(data.table)

# SCRIPT PARAMETERS
# -----------------
cols      <- c("eid","31-0.0","50-0.0")
col_names <- c("id","sex","height")

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

# WRITE DATA TO FILE
# ------------------
write.csv(dat,"",row.names = FALSE,col.names = FALSE,quote = FALSE)
