#!/usr/bin/env Rscript
arg = commandArgs(trailingOnly = TRUE)
region.name = as.character(arg)

# Load pacakges
library(data.table)
library(readr)
library(Matrix)

# phenotpe file and genotype file names
pheno.file <- "/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/height.csv.gz"
geno.file = paste0('height.', region.name, '.0.01.raw.gz')

# Read genotype data
cat("Reading genotype data.\n")
geno <- fread(geno.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(geno) <- "data.frame"
# Extract the genotypes.
X <- as(as.matrix(geno[-(1:6)]),'dgCMatrix')

# Read phenotype data
cat("Reading phenotype data.\n")
pheno        <- suppressMessages(read_csv(pheno.file))
class(pheno) <- "data.frame"
pheno$sex = factor(pheno$sex)
pheno$assessment_centre = factor(pheno$assessment_centre)
pheno$genotype_measurement_batch = factor(pheno$genotype_measurement_batch)
pheno$age2 = pheno$age^2
## match individual order with genotype file
ind = fread(paste0('height.', region.name, '.0.01.psam'))
match.idx = match(ind$IID, pheno$id)
pheno = pheno[match.idx,]

# Covariates
Z = model.matrix(~ sex + age + age2 + assessment_centre + genotype_measurement_batch +
                   pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
                   pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10 +
                   pc_genetic11 + pc_genetic12 + pc_genetic13 + pc_genetic14 + pc_genetic15 +
                   pc_genetic16 + pc_genetic17 + pc_genetic18 + pc_genetic19 + pc_genetic20, data = pheno)

## Remove intercept
Z = Z[,-1]
## standardize quantitative columns
cols = which(colnames(Z) %in% c("age","pc_genetic1","pc_genetic2","pc_genetic3","pc_genetic4",
                                "pc_genetic5","pc_genetic6","pc_genetic7","pc_genetic8","pc_genetic9", 
                                "pc_genetic10","pc_genetic11","pc_genetic12","pc_genetic13","pc_genetic14",
                                "pc_genetic15","pc_genetic16","pc_genetic17","pc_genetic18","pc_genetic19","pc_genetic20"))
Z[,cols] = scale(Z[,cols])
Z[,'age2'] = Z[,'age']^2

# Compute XtX and Xty
y = pheno$height
names(y) = pheno$id

# Center y
y = y - mean(y)
# Center scale X
cm = Matrix::colMeans(X, na.rm = TRUE)
csd = susieR:::compute_colSds(X)
csd[csd == 0] = 1
X = as.matrix(t((t(X) - cm) / csd))

A   <- crossprod(Z) # Z'Z
# chol decomposition for (Z'Z)^(-1)
R = chol(solve(A)) # R'R = (Z'Z)^(-1)
W = R %*% crossprod(Z, X) # RZ'X
S = R %*% crossprod(Z, y) # RZ'y

# Load LD matrix from raw genotype
ld.matrix = as.matrix(fread(paste0('height.', region.name, '.0.01.matrix')))
# X'X
XtX = ld.matrix*(nrow(X)-1) - crossprod(W) # W'W = X'ZR'RZ'X = X'Z(Z'Z)^{-1}Z'X
rownames(XtX) = colnames(XtX) = colnames(X)
# X'y
Xty = as.vector(y %*% X)
Xty = Xty - crossprod(W, S) # W'S = X'ZR'RZ'y = X'Z(Z'Z)^{-1}Z'y

## SNP info
maf <- read.delim(paste0('height.', region.name, '.0.01.afreq'))
pos <- fread(paste0('height.', region.name, '.0.01.pvar'))
pos$maf = pmin(maf$ALT_FREQS, 1-maf$ALT_FREQS)

## Save results
saveRDS(list(XtX = XtX, Xty = Xty, yty = sum(y^2) - crossprod(S), n = length(y), pos=pos), 
        paste0('height.', region.name, '.0.01.XtX.Xty.rds'))

