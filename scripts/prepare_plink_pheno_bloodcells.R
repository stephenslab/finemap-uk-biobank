# We create phenotype file and covariates files for plink GWAS.


# SCRIPT PARAMETERS
# -----------------
# CSV containing the phenotype and covariate data
pheno.file     <- file.path("/gpfs/data/stephens-lab/finemap-uk-biobank",
                            "data/raw/BloodCells/bloodcells.csv")
# Output phenotype file and covariants file names
out.pheno.file <- file.path("/gpfs/data/stephens-lab/finemap-uk-biobank",
                            "data/raw/BloodCells/bloodcells.pheno.txt")
out.covar.file <- file.path("/gpfs/data/stephens-lab/finemap-uk-biobank",
                            "data/raw/BloodCells/bloodcells.covar.txt")

# SET UP ENVIRONMENT
# ------------------
library(readr)
library(car)
library(dplyr)

# LOAD PHENOTYPE and COVARIATES DATA
# -------------------
# Load the previously prepared phenotype data from the CSV file.
cat("Reading phenotype data.\n")
dat         <- suppressMessages(read_csv(pheno.file))
class(dat) <- "data.frame"
dat$sex = factor(dat$sex)
dat$assessment_centre = factor(dat$assessment_centre)
dat$genotype_measurement_batch = factor(dat$genotype_measurement_batch)
dat$age2 = dat$age^2

pheno_names = c("WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW", "Platelet_count",
                "Plateletcrit", "PDW", "Lymphocyte_perc", "Monocyte_perc",
                "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
                "Reticulocyte_perc", "MSCV", "HLR_perc")
# Save phenotype file
cat("Writing phenotype data in the format accepted by PLINK.\n")
pheno.table <- dat %>% mutate(FID = id, IID = id) %>% select(FID, IID, pheno_names)
write.table(pheno.table, file = out.pheno.file, quote = FALSE,
            row.names = FALSE)

# Create covariates matrix.
Z = model.matrix(~ sex + age + age2 + assessment_centre + genotype_measurement_batch +
                   pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
                   pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10, data = dat)

# delete intercept
Z = Z[,-1]

# standardize quantitative columns
cols = which(colnames(Z) %in% c("age",paste0("pc_genetic", 1:10)))
Z[,cols] = scale(Z[,cols])
Z[,'age2'] = Z[,'age']^2

# Save covariates file
IID = dat$id
FID = dat$id
cov.table = cbind(FID, IID, Z)
write.table(cov.table, file= out.covar.file, quote=FALSE, row.names = FALSE)

# PLINK command
# plink2 --linear hide-covar no-x-sex omit-ref \
# --pfile height.chr3.140800000.141800000.0.01 \
# --covar covar.height.txt --pheno pheno.height.txt \
# --vif 100 \
# --out height.chr3.140800000.141800000.0.01.plink2

# PLINK computes vif for each column of Z, which treat factor levels
# separately. Some vif for one level of factor is very high.  Using R,
# we can compute one vif for all factor levels. In this case, vif is
# close to 1.
vif.Z = vif(lm(dat$WBC_count ~ sex1 + age + age2 + dat$assessment_centre + dat$genotype_measurement_batch +
                 pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
                 pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10, data = as.data.frame(Z)))
