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
               "30000-0.0", "30010-0.0", "30020-0.0", "30040-0.0",
               "30070-0.0", "30080-0.0", "30090-0.0",
               "30110-0.0", "30180-0.0", "30190-0.0", "30200-0.0",
               "30210-0.0", "30220-0.0", "30240-0.0",
               "30270-0.0", "30290-0.0",
               paste0("41202-0.", 0:379), "3140-0.0")
col_names <- c("id","sex","assessment_centre","age","ethnic_genetic",
               "sex_genetic","genotype_measurement_batch","missingness",
               paste0("pc_genetic",1:10),
               "kinship_genetic", "outliers",
               "WBC_count", "RBC_count", "Haemoglobin", "MCV",
               "RDW", "Platelet_count", "Plateletcrit",
               "PDW", "Lymphocyte_perc", "Monocyte_perc", "Neutrophill_perc",
               "Eosinophill_perc", "Basophill_perc", "Reticulocyte_perc",
               "MSCV", "HLR_perc",
               paste0('ICD10.',0:379), "pregnancy")

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
## There are 502492 samples

# PREPARE DATA
# ------------
# Select the requested columns.
cat("Preparing data.\n")
dat        <- dat[,cols]
names(dat) <- col_names

# Convert numerical columns except the first one (the first column contains
# the sample ids) to numeric values, and set all empty strings to NA.
n <- length(dat)
for (i in 2:n) {
  x          <- dat[,i]
  x[x == ""] <- as.character(NA)
  if(grepl('ICD10', colnames(dat)[i])){
    dat[,i]    <- x
  }else{
    dat[,i]    <- as.numeric(x)
  }
}

# Remove any samples that are not marked as being "White British".
# This step should filter out 92887 rows.
dat = dat %>% filter(ethnic_genetic == 1)
cat(sprintf("After removing non White British, %d rows remain.\n",nrow(dat)))

# Remove all rows in which one or more of the values are missing,
# aside from the in the "outlier", "ICD10", "pregnancy" columns.
# The "outliers" have value 1 when it is an outlier, NA otherwise.
# The "pregnancy" have value NA for males.
# This step should filter out 18578 rows.
cols <- !(grepl(paste(c('ICD10', "outliers", "pregnancy"),collapse = '|'), names(dat)))
rows <- which(rowSums(is.na(dat[,cols])) == 0)
dat  <- dat[rows,]
cat(sprintf("After removing rows with NAs, %d rows remain.\n",nrow(dat)))

# Remove rows with mismatches between self-reported and genetic sex
# This step should filter out 287 rows.
dat <- dat %>% filter(sex == sex_genetic)
cat(sprintf("After removing sex mismatches, %d rows remain.\n",nrow(dat)))

# Remove "missingness" and "heterozygosity" outliers as defined by UK
# Biobank. This step should filter out 665 rows. Note that this step
# will remove any samples in which the "missingness" column is greater
# than 5%.
dat <- dat %>% filter(is.na(outliers))
cat(sprintf("After removing outliers, %d rows remain.\n",nrow(dat)))

# Remove any individuals have at leat one relative based on the
# kinship calculations. This step should filter out 126,236 rows.
dat <- dat %>% filter(kinship_genetic == 0)
cat(sprintf(paste("After removing relatedness individuals based on kinship,",
                  "%d rows remain.\n"),nrow(dat)))

# Remove any pregnant individuals
# This step should filter out 164 rows.
dat <- dat %>% filter(!(pregnancy %in% c(1,2)))
cat(sprintf(paste("After removing pregnant individuals,",
                  "%d rows remain.\n"),nrow(dat)))

# Remove any individuals with blood related diseases
# This step should filter out 6070 rows
icd10 = c('C94', 'C95', 'Z856', "C901", "C914", "C82", "C83", 'C84', "C85", "Z948",
  "Z511", "Z512", "Z542", "D46", paste0("D", 55:64), paste0("B", 20:24),
  "N180", "Z992", "Z491", "Z492", "K74", "C88", "C900", "C902", "C91", "C92",
  "D45", "D47", "E831")
daticd10 = dat %>% select(which(grepl('ICD10', colnames(dat)))) %>% as.matrix
icd_status = matrix(grepl(paste(icd10, collapse='|'), daticd10), nrow(daticd10), ncol(daticd10))
dat = dat %>% filter(rowSums(icd_status) == 0)
cat(sprintf(paste("After removing individuals with blood diseases,",
                  "%d rows remain.\n"),nrow(dat)))

# Remove individuals with "abnormal" measurements.
# This step should filter out 8625 rows
pheno_names = c("WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW", "Platelet_count",
                "Plateletcrit", "PDW", "Lymphocyte_perc", "Monocyte_perc",
                "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
                "Reticulocyte_perc", "MSCV", "HLR_perc")
## RINT
for(name in pheno_names){
  id = which(colnames(dat) == name)
  dat[,id] = qnorm((rank(dat[,id],na.last="keep")-0.5)/sum(!is.na(dat[,id])))
}
## compute empirical covariance matrix
covy = dat %>% select(pheno_names) %>% cov
D2 = stats::mahalanobis(dat %>% select(pheno_names), center=0, cov=covy) ## mahalanobis distance
dat = dat[D2 < qchisq(0.01, df=16, lower.tail = F),]
cat(sprintf(paste("After removing individuals with abnormal measurements,",
                  "%d rows remain.\n"),nrow(dat)))

# Finally, remove the columns that are no longer needed for subsequent
# analyses.
cols.to.remove <- c("sex_genetic","ethnic_genetic",
                    "missingness",
                    "kinship_genetic","outliers","pregnancy",paste0("ICD10.", 0:379))
cols <- which(!is.element(names(dat),cols.to.remove))
dat  <- dat[,cols]

# SUMMARIZE DATA
# --------------
# Double-check that everything looks okay.
summary(dat)

cat("Writing prepared data to CSV file.\n")
write.csv(dat,output.file,row.names = FALSE,quote = FALSE, na='NA')


