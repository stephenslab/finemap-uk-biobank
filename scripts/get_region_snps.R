#!/usr/bin/env Rscript
arg = commandArgs(trailingOnly = TRUE)
chr = as.numeric(arg[1])
start.pos = as.numeric(arg[2])
stop.pos = as.numeric(arg[3])
region.name = as.character(arg[4])

print(arg)

print(chr)
# Script to generate a list of all the ids of the genotyped genetic
# variants within a specified base-pair region on a chromosome, and
# combine this list with other pertinent info, such as base-pair
# position and GeneATLAS association statistics.

# SCRIPT PARAMETERS
# -----------------
# This region on chromosome 3 contains one of the strongest
# associations for standing height; see http://big.stats.ox.ac.uk.
bim.file  <- file.path("/gpfs/data/pierce-lab/uk-biobank-genotypes",
                       paste0("ukb_imp_chr",chr,"_v3.bim.gz"))
out.file  <- paste0("../data/region-variants-",region.name,".csv")

# These files contain the GeneATLAS and Neale lab association results.
geneatlas.file <-
  file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/summary",
            "gene-atlas/height/imputed/not-normalized",
            paste0("imputed.allWhites.50-0.0.chr", chr,".csv.gz"))
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

# PLOT GENEATLAS VS. NEALE ASSOCIATIONS
# -------------------------------------
# Create a scatterplot comparing the t-statistics.
p1 <- ggplot(map,aes(x = NBETA/NSE,y = neale_tstat)) +
  geom_point(shape = 20,size = 2,na.rm = TRUE) +
  geom_abline(intercept = 0,slope = 1,color = "dodgerblue",
              linetype = "dotted") +
  xlim(c(-25,50)) +
  ylim(c(-25,50)) +
  theme_cowplot(font_size = 12) +
  labs(x = "GeneATLAS",y = "Neale lab",title = "t-statistics")

# Create a scatterplot comparing the p-values (on the logarithmic scale).
p2 <- ggplot(map,aes(x = PV + 1e-175,y = neale_pval + 1e-175)) +
  geom_point(shape = 20,size = 2,na.rm = TRUE) +
  geom_abline(intercept = 0,slope = 1,color = "dodgerblue",
              linetype = "dotted") +
  scale_x_continuous(limits = c(10^(-175),1),trans = "log10") +
  scale_y_continuous(limits = c(10^(-175),1),trans = "log10") +
  theme_cowplot(font_size = 12) +
  labs(x = "GeneATLAS",y = "Neale lab",title = "p-values")

# WRITE SNP DATA to FILE
# ----------------------
cat("Writing SNP data to file.\n")
write.csv(map,out.file,quote = FALSE,row.names = FALSE)
