# TO DO: Provide a more detailed description of what this script is
# for, and how to use it.

# Create phenotype and covariates files

# SCRIPT PARAMETERS
# -----------------
# CSV containing the phenotype and covariate data
pheno.file     <- file.path("/gpfs/data/stephens-lab/finemap-uk-biobank",
                            "data/raw/height.csv.gz")
# Output phenotype file and covariants file names
out.pheno.file <- "height.pheno.txt"
out.covar.file <- "height.covar.txt"

# SET UP ENVIRONMENT
# ------------------
library(readr)
# library(car)
# library(dplyr)

# LOAD PHENOTYPE and COVARIATES DATA
# -------------------
# Load the previously prepared phenotype data from the CSV file.
cat("Reading phenotype data.\n")
pheno        <- suppressMessages(read_csv(pheno.file))
class(pheno) <- "data.frame"
pheno$sex = factor(pheno$sex)
pheno$assessment_centre = factor(pheno$assessment_centre)
pheno$genotype_measurement_batch = factor(pheno$genotype_measurement_batch)
pheno$age2 = pheno$age^2
# match individual order with genotype file
ind = fread('height.chr3.140800000.141800000.0.01.psam')
match.idx = match(ind$IID, pheno$id)
pheno = pheno[match.idx,]

# Save phenotype file
cat("Writing phenotype data in the format accepted by PLINK.\n")
pheno.table <- data.frame(FID    = pheno$id,
                          IID    = pheno$id,
                          height = pheno$height)
write.table(pheno.table, file = out.pheno.file, quote = FALSE,
            row.names = FALSE) 

# Create covariates matrix.
Z = model.matrix(~ sex + age + age2 + assessment_centre + genotype_measurement_batch +
                   pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
                   pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10 +
                   pc_genetic11 + pc_genetic12 + pc_genetic13 + pc_genetic14 + pc_genetic15 +
                   pc_genetic16 + pc_genetic17 + pc_genetic18 + pc_genetic19 + pc_genetic20, data = pheno)

# delete intercept
Z = Z[,-1]

# standardize quantitative columns
cols = which(colnames(Z) %in% c("age","pc_genetic1","pc_genetic2","pc_genetic3","pc_genetic4",
                                "pc_genetic5","pc_genetic6","pc_genetic7","pc_genetic8","pc_genetic9", 
                                "pc_genetic10","pc_genetic11","pc_genetic12","pc_genetic13","pc_genetic14",
                                "pc_genetic15","pc_genetic16","pc_genetic17","pc_genetic18","pc_genetic19","pc_genetic20"))
Z[,cols] = scale(Z[,cols])
Z[,'age2'] = Z[,'age']^2

# Save covariates file
IID = pheno$id
FID = pheno$id
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
vif.Z = vif(lm(pheno$height~ sex + age + age2 + pheno$assessment_centre + pheno$genotype_measurement_batch +
                 pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
                 pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10 +
                 pc_genetic11 + pc_genetic12 + pc_genetic13 + pc_genetic14 + pc_genetic15 +
                 pc_genetic16 + pc_genetic17 + pc_genetic18 + pc_genetic19 + pc_genetic20, data = as.data.frame(Z)))
