# TO DO: Explain here what this script does, and how to use it.

# SCRIPT PARAMETERS
# -----------------
# This region on chromosome 3 contains one of the strongest
# associations for standing height; see http://big.stats.ox.ac.uk.
bim.file  <- file.path("/gpfs/data/pierce-lab/uk-biobank-genotypes",
                       "ukb_imp_chr3_v3.bim.gz")
out.file  <- "../data/region-variants-ZBTB38.txt"
start.pos <- 140800000
stop.pos  <- 141800000

# SET UP ENVIRONMENT
# ------------------
library(data.table)

# LOAD SNP DATA
# -------------
cat("Loading SNP data from PLINK file.\n")
map <- fread(bim.file,sep = "\t",header = FALSE,verbose = FALSE,
             showProgress = FALSE)
class(map) <- "data.frame"
names(map) <- c("chr","id","cM","pos","A1","A2")

# GET SNP IDs
# -----------
# Get the set of SNPs in which their base-pair positions are within
# the specified range.
cat("Getting SNPs within chosen region.\n")
map <- subset(map,pos > start.pos & pos < stop.pos)
map <- subset(map,!duplicated(id))

# WRITE SNP IDs to FILE
# ---------------------
cat("Writing SNP ids to file.\n")
write.table(map["id"],out.file,quote = FALSE,col.names = FALSE,
            row.names = FALSE)
