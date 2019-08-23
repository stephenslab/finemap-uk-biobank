# Script to prepare the UK Biobank phenotype data for analyzing height
# associations. There are multiple input files because we made
# multiple requests for phenotype data. The output is a CSV file
# containing data for the selected phenotypes and samples. This should
# take at most 1-2 hours to run, and may require a good amount of
# memory (works well with about 24 GB, although may not require this
# much).

# SCRIPT PARAMETERS
# -----------------
input.file1 <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                         "12-feb-2019","ukb26140.csv.gz")
input.file2 <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                         "11-jun-2019","ukb32141.csv.gz")
output.file <- file.path("/gpfs/data/stephens-lab/finemap-uk-biobank",
                         "data/raw/height.csv")

# The columns selected for subsequent analyses are as follows:
#
#   sex (31)
#   height (50)
#   UK Biobank assessment centre (54)
#   age (21022)
#   genetic ethnic grouping (22006)
#   genetic sex (22001)
#   genotype measurement batch (22000)
#   missingness (22005)
#   genetic PCs (22009-0.1 - 22009-0.40)
#   genetic relatedness pairing (22011)
#   genetic kinship (22021)
#   outliers (22027)
#
cols      <- c("eid","31-0.0","50-0.0","54-0.0","21022-0.0","22006-0.0",
               "22001-0.0","22000-0.0","22005-0.0",paste0("22009-0.",1:40),
               paste0("22011-0.",0:4),"22021-0.0", "22027-0.0")
col_names <- c("id","sex","height","assessment_centre","age","ethnic_genetic",
               "sex_genetic","genotype_measurement_batch","missingness",
               paste0("pc_genetic",1:40),paste0("relatedness_genetic",0:4),
               "kinship_genetic", "outliers")

# SET UP ENVIRONMENT
# ------------------
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

# LOAD DATA
# ---------
# After loading the two tables from the CSV files, the tables are
# merged by the "eid" column. The second table has slightly fewer rows
# than the first table, and the sample ids in the second table are a
# subset of the ids in the first table.
cat("Reading data from the CSV files.\n")
out <- system.time({
  dat1 <- fread(input.file1,sep = ",",header = TRUE,verbose = FALSE,
                showProgress = FALSE,colClasses = "character");
  dat2 <- fread(input.file2,sep = ",",header = TRUE,verbose = FALSE,
                showProgress = FALSE,colClasses = "character")
})
class(dat1) <- "data.frame"
class(dat2) <- "data.frame"
cat(sprintf("Data loading step took %d seconds.\n",round(out["elapsed"])))
dat <- inner_join(dat1,dat2,by = "eid")
rm(dat1,dat2)
cat(sprintf("Merged table contains %d rows.\n",nrow(dat)))

# PREPARE DATA
# ------------
# Select the requested columns.
cat("Preparing data.\n")
dat        <- dat[,cols]
names(dat) <- col_names

# Convert all columns except the first one (the first column contains
# the sample ids) to numeric values, and set all empty strings to NA.
n <- length(dat)
for (i in 2:n) {
  x          <- dat[,i]
  x[x == ""] <- as.character(NA)
  dat[,i]    <- as.numeric(x)
}

# Remove all rows in which one or more of the values are missing,
# aside from the in the "outlier" and "relatedness_genetic" columns.
#
# When the "genetic ethnic grouping" column is included, this removes
# any samples that are not marked as being "White British". The
# "outliers" have value 1 when it is an outlier, NA otherwise.
cols <- !(names(dat) == "outliers" | grepl("relatedness_genetic",names(dat)))
rows <- which(rowSums(is.na(dat[,cols])) == 0)
dat  <- dat[rows,]
cat(sprintf("After removing rows with NAs, %d rows remain.\n",nrow(dat)))

# Remove rows with mismatches between self-reported and genetic sex
# This step should filter out 310 rows.
dat <- dat %>% filter(sex == sex_genetic)
cat(sprintf("After removing sex mismatches, %d rows remain.\n",nrow(dat)))

# Remove "missingness" and "heterozygosity" outliers as defined by UK
# Biobank. This step should filter out 723 rows. Note that this step
# will remove any samples in which the "missingness" column is greater
# than 5%.
dat <- dat %>% filter(is.na(outliers))
cat(sprintf("After removing outliers, %d rows remain.\n",nrow(dat)))

# Remove any individuals have at leat one relative based on the
# kinship calculations. This step should filter out 131,805 rows.
dat <- dat %>% filter(kinship_genetic == 0)
cat(sprintf(paste("After removing relatedness individuals based on kinship,",
                  "%d rows remain.\n"),nrow(dat)))

# Remove any individuals that have close relatives identified from the
# "relatendess" calculations that weren't already identified using the
# kinship calculations. This step should filter out 394 rows.
dat <- dat %>% filter(is.na(relatedness_genetic0))
cat(sprintf("After removing relatedness individuals, %d rows remain.\n",
            nrow(dat)))

# Remove individuals with "abnormal" height. This step should filter
# out 976 rows.
fit  <- lm(height ~ sex,dat)
y    <- resid(fit)
rows <- which(y >= median(y) - 3*sd(y) & y <= median(y) + 3*sd(y))
dat  <- dat[rows,]
cat(sprintf("After removing samples with abnormal height, %d rows remain.\n",
            nrow(dat)))

# Finally, remove the columns that are no longer needed for subsequent
# analyses.
cols.to.remove <- c("sex_genetic","ethnic_genetic",
                    paste0("relatedness_genetic",0:4),
                    "kinship_genetic","outliers")
cols <- which(!is.element(names(dat),cols.to.remove))
dat  <- dat[,cols]

# SUMMARIZE DATA
# --------------
# Double-check that everything looks okay.
summary(dat)

# WRITE DATA TO FILE
# ------------------
cat("Writing prepared data to CSV file.\n")
write.csv(dat,output.file,row.names = FALSE,quote = FALSE)
