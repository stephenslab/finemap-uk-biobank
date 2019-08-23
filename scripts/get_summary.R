# Before running this script, you will need to generate the PLINK
# files containing the genotypes for the selected SNP(s). This can be
# done on the gardner cluster with the following plink command:
#
# zcat /gpfs/data/stephens-lab/finemap-uk-biobank/height.csv.gz | awk 'BEGIN{FS=","} {print $1 OFS $1}' | tail -n +2 > height.id.txt
#
# plink2 --sample /gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample \
# --bgen /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr1_v3.bgen \
# --extract imputed.allWhites.chr1.2000.txt \
# --keep height.id.txt \
# --maf 0.05 \
# --make-pgen --threads 8 --memory 38000 \
# --out height.chr1.imputed.2000 > plink2.imputed.log
#
# plink2 --pfile height.chr1.imputed.2000 --threads 8 --export A-transpose --out height.chr1.imputed.2000
#
# gzip height.chr1.imputed.2000.traw

library(readr)
library(broom)
library(data.table)

# SCRIPT PARAMETERS
# -----------------
pheno.file <- "/gpfs/data/stephens-lab/finemap-uk-biobank/height.csv.gz"
geno.file  <- "height.chr1.imputed.2000.traw.gz"

# LOAD PHENOTYPE DATA
# -------------------
# Read the phenotype data from the CSV file.
cat("Reading phenotype data.\n")
pheno        <- suppressMessages(read_csv(pheno.file))
class(pheno) <- "data.frame"
rownames(pheno) <- pheno$id
pheno$sex = factor(pheno$sex)
pheno$assessment_centre = factor(pheno$assessment_centre)
pheno$genotype_measurement_batch = factor(pheno$genotype_measurement_batch)
pheno = pheno[,-1]

# LOAD GENOTYPE DATA
# ------------------
# Read the genotype data stored in the binary PLINK ("bed") file.
# Load the data from the .traw file.
cat("Reading genotype data.\n")
geno <- fread(geno.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(geno) <- "data.frame"

# Extract the SNP information.
map <- geno[1:6]

# Extract the genotypes.
geno <- geno[-(1:6)]
geno <- t(as.matrix(geno))
rownames(geno) = gsub("_.*", "", rownames(geno))

# ALIGN PHENOTYPES AND GENOTYPES
# ------------------------------
# Get the set of samples that are common to both the phenotype and
# genotype data.
cat("Aligning phenotypes and genotypes.\n")
ids   <- intersect(rownames(pheno),rownames(geno))
ids   <- sort(as.integer(ids))
ids   <- as.character(ids)
rows1 <- match(ids,rownames(pheno))
rows2 <- match(ids,rownames(geno))
pheno <- pheno[rows1,]
geno  <- geno[rows2,]

# ESTIMATE EFFECT OF GENOTYPE ON HEIGHT
# -------------------------------------
cat("Linear relationship between genotype and phenotype:\n")
beta = se = pv = numeric(ncol(geno))
for(i in 1:ncol(geno)){
  m = lm(height ~ geno[,i] + sex + age + I(age^2) + assessment_centre + genotype_measurement_batch +
           pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
           pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10 +
           pc_genetic11 + pc_genetic12 + pc_genetic13 + pc_genetic14 + pc_genetic15 +
           pc_genetic16 + pc_genetic17 + pc_genetic18 + pc_genetic19 + pc_genetic20, data = pheno)
  beta[i] = tidy(m)$estimate[2]
  se[i] = tidy(m)$std.error[2]
  pv[i] = tidy(m)$p.value[2]
  if(i %% 10 == 0){
    cat(paste0('Running ', i, '\n'))
  }
}

cat('Saving result')
result = cbind(beta, se, pv)
rownames(result) = map$SNP
saveRDS(result, 'height.chr1.imputed.2000.result.rds')
