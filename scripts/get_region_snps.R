# Script to generate a list of all the ids of the genotyped genetic
# variants within a specified base-pair region on a chromosome.
#
# TO DO: Update this description.
#

# SCRIPT PARAMETERS
# -----------------
# This region on chromosome 3 contains one of the strongest
# associations for standing height; see http://big.stats.ox.ac.uk.
chr       <- 3
bim.file  <- file.path("/gpfs/data/pierce-lab/uk-biobank-genotypes",
                       "ukb_imp_chr3_v3.bim.gz")
out.file  <- "../data/region-variants-ZBTB38.csv"
start.pos <- 140.8e6
stop.pos  <- 141.8e6

# These files contain the GeneATLAS and Neale lab association results.
geneatlas.file <-
  file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/summary",
            "gene-atlas/height/imputed/not-normalized",
            "imputed.allWhites.50-0.0.chr3.csv.gz")
neale.file <-
  file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/summary/neale",
            "50_irnt.gwas.imputed_v3.both_sexes.tsv.gz")

# SET UP ENVIRONMENT
# ------------------
library(readr)
library(ggplot2)
library(cowplot)

# LOAD SNP DATA
# -------------
cat("Loading SNP data from PLINK file.\n")
map <- suppressMessages(read_delim(bim.file,delim = "\t",col_names = FALSE,
                                   progress = FALSE))
names(map) <- c("chr","id","cM","pos","A1","A2")
class(map) <- "data.frame"

# LOAD GENEATLAS ASSOCIATION RESULTS
# ----------------------------------
# Load the summary statistics from the space-delimited text file.
cat("Loading GeneATLAS association results.\n")
geneatlas <- suppressMessages(read_delim(geneatlas.file,delim = " ",
                                         progress = FALSE))
names(geneatlas) <- c("SNP","ALLELE","iscores","NBETA","NSE","PV")
class(geneatlas) <- "data.frame"

# LOAD NEALE LAB ASSOCIATION RESULTS
# ----------------------------------
# Load the summary statistics from the tab-delimited text file.
cat("Loading Neale lab association results.\n")
neale <- suppressMessages(read_delim(neale.file,delim = "\t",progress = FALSE))
class(neale) <- "data.frame"
neale        <- neale[c("variant","minor_AF","beta","se","tstat","pval")]

# Retain the SNPs on the selected chromosome.
out       <- strsplit(neale$variant,":",fixed = TRUE)
neale$chr <- factor(sapply(out,function (x) x[[1]]))
neale$pos <- as.numeric(sapply(out,function (x) x[[2]]))
neale     <- neale[neale$chr == chr,]

# SELECT SNPs
# -----------
# Get the set of SNPs in which their base-pair positions are within
# the specified range.
cat("Getting SNPs within chosen region.\n")
map <- subset(map,pos > start.pos & pos < stop.pos)
map <- subset(map,!duplicated(id))
map <- map[c("chr","id","pos")]
rownames(map) <- NULL
cat(sprintf("%d SNPs are selected.\n",nrow(map)))

# Incorporate the GeneATLAS and Neale lab association statistics into
# the table. Note that only SNPs with at least an association
# statistic computed by at least one group (GeneATLAS or Neale lab)
# are retained.
geneatlas    <- geneatlas[c("SNP","ALLELE","NBETA","NSE","PV")]
neale        <- neale[c("pos","beta","se","tstat","pval")]
names(geneatlas)[1] <- "id"
names(neale) <- c("pos","neale_beta","neale_se","neale_tstat","neale_pval")
map          <- merge(map,geneatlas,by = "id")
map          <- merge(map,neale,by = "pos",all.x = TRUE)
map          <- map[c("chr","pos","id","ALLELE","NBETA","NSE","PV",
                      "neale_beta","neale_se","neale_tstat","neale_pval")]
cat(sprintf("%d SNPs with association statistics are retained.\n",nrow(map)))

# TO DO: Create a scatterplot comparing the t-statistics and the p-values.

# WRITE SNP DATA to FILE
# ----------------------
cat("Writing SNP data to file.\n")
write.table(map,out.file,quote = FALSE,row.names = FALSE)
