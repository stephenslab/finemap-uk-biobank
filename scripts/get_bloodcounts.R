library(data.table)
library(dplyr)
input.file1 <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                         "12-feb-2019","ukb26140.csv.gz")
input.file2 <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                         "11-jun-2019","ukb32141.csv.gz")
input.file3 <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
                         "16-oct-2020","ukb44231.csv.gz")
output.file <- file.path("/gpfs/data/stephens-lab/finemap-uk-biobank",
                         "data/raw/bloodcounts.csv")

cols      <- c("eid","31-0.0","54-0.0","21022-0.0","22006-0.0",
               "22001-0.0","22000-0.0","22005-0.0",paste0("22009-0.",1:10),
               "22021-0.0", "22027-0.0",
               "30000-0.0", "30010-0.0", "30020-0.0", "30030-0.0", "30040-0.0",
               "30050-0.0", "30060-0.0", "30070-0.0", "30080-0.0", "30090-0.0",
               "30100-0.0", "30110-0.0", "30120-0.0", "30130-0.0", "30140-0.0",
               "30150-0.0", "30160-0.0", "30180-0.0", "30190-0.0", "30200-0.0",
               "30210-0.0", "30220-0.0", "30240-0.0", "30250-0.0", "30260-0.0",
               "30270-0.0", "30280-0.0", "30290-0.0", "30300-0.0")
col_names <- c("id","sex","assessment_centre","age","ethnic_genetic",
               "sex_genetic","genotype_measurement_batch","missingness",
               paste0("pc_genetic",1:10),
               "kinship_genetic", "outliers",
               "WBC#", "RBC#", "Haemoglobin", "Haematocrit", "MCV",
               "MCH", "MCHC", "RDW", "Platelet#", "Plateletcrit",
               "MPV", "PDW", "Lymphocyte#", "Monocyte#", "Neutrophill#",
               "Eosinophill#", "Basophill#", "Lymphocyte%", "Monocyte%", "Neutrophill%",
               "Eosinophill%", "Basophill%", "Reticulocyte%", "Reticulocyte#", "MRV",
               "MSCV", "IRF", "HLR%", "HLR#")

cat("Reading data from the CSV files.\n")
out <- system.time({
  dat1 <- fread(input.file1,sep = ",",header = TRUE,verbose = FALSE,
                showProgress = FALSE,colClasses = "character");
  dat2 <- fread(input.file2,sep = ",",header = TRUE,verbose = FALSE,
                showProgress = FALSE,colClasses = "character");
  dat3 <- fread(input.file3,sep = ",",header = TRUE,verbose = FALSE,
                showProgress = FALSE,colClasses = "character")
})
class(dat1) <- "data.frame"
class(dat2) <- "data.frame"
class(dat3) <- "data.frame"
cat(sprintf("Data loading step took %d seconds.\n",round(out["elapsed"])))
dat12 <- inner_join(dat1,dat2,by = "eid")
dat <- inner_join(dat12, dat3, by='eid')
rm(dat1,dat2,dat3,dat12)
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
cols <- !(names(dat) == "outliers")
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

# Remove individuals with "abnormal" measurements.
# ....

# Finally, remove the columns that are no longer needed for subsequent
# analyses.
cols.to.remove <- c("sex_genetic","ethnic_genetic",
                    "missingness",
                    "kinship_genetic","outliers")
cols <- which(!is.element(names(dat),cols.to.remove))
dat  <- dat[,cols]

cat("Writing prepared data to CSV file.\n")
write.csv(dat,output.file,row.names = FALSE,quote = FALSE)


