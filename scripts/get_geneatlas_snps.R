# TO DO: Update this description.
#
# Create a table containing the top 1,000 GeneATLAS height
# associations from chromosome 1, and an additional 1,000 "null" SNPs
# (SNPs that are not strongly associated with height), and save this
# table in a CSV file.

# SCRIPT PARAMETERS
# -----------------
# Select this many strongly associated SNPs, and this many SNPs that
# are not strongly associated ("null" SNPs).
n <- 1000

# Select SNPs from this chromosome.
chr <- 1

# PLINK ".bim" file containing the information about the SNPs on the
# selected chromosome.
bim.file <- file.path("/gpfs/data/pierce-lab/uk-biobank-genotypes",
                      "ukb_imp_chr1_v3.bim.gz")

# These files contain the GeneATLAS and Neale lab association results.
geneatlas.file <-
  file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/summary",
            "gene-atlas/height/imputed/not-normalized",
            "imputed.allWhites.50-0.0.chr1.csv.gz")
neale.file <-
  file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/summary/neale",
            "50_irnt.gwas.imputed_v3.both_sexes.tsv.gz")

# The table of association results for the 2n selected SNPs is saved
# to this CSV file.
out.file <- "../data/geneatlas-neale-height.csv"

# SET UP ENVIRONMENT
# ------------------
# Load the packages used in the analysis below.
library(readr)
library(ggplot2)
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD SNP DATA
# -------------
cat("Loading SNP data from PLINK file.\n")
map <- suppressMessages(read_delim(bim.file,delim = "\t",col_names = FALSE,
                                   progress = FALSE))
class(map) <- "data.frame"
names(map) <- c("chr","SNP","cM","pos","A1","A2")
map        <- map[c("chr","pos","SNP")]

# Remove SNPs with multiple entries.
rows <- which(duplicated(map$pos))
snps <- map[rows,"pos"]
map  <- subset(map,!is.element(pos,snps))

# LOAD GeneATLAS ASSOCIATION RESULTS
# ----------------------------------
# Load the summary statistics from the space-delimited text file.
cat("Loading GeneATLAS association results.\n")
geneatlas <- suppressMessages(read_delim(geneatlas.file,delim = " ",
                                         progress = FALSE))
names(geneatlas) <- c("SNP","ALLELE","iscores","NBETA","NSE","PV")
class(geneatlas) <- "data.frame"

# Add the SNP base-pair positions to this table.
geneatlas <- merge(map,geneatlas,by = "SNP",sort = FALSE)

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

# Add these statistics to the GeneATLAS table.
neale        <- neale[c("minor_AF","beta","se","tstat","pval","pos")]
names(neale) <- c("minor_AF","neale_beta","neale_se","neale_tstat",
                  "neale_pval","pos")
geneatlas    <- merge(geneatlas,neale,by = "pos",sort = FALSE)
geneatlas    <- geneatlas[c(3,1,2,4:13)]

# Filter out SNPs with low minor allele frequencies.
geneatlas <- subset(geneatlas,minor_AF >= 0.01)

# SELECT ASSOCIATION RESULTS
# --------------------------
# Select the most strongly associated SNPs.
cat("Selecting SNPs.\n")
geneatlas.top <- top_n(geneatlas,-n,PV)

# Randomly select "null" SNPs.
rows           <- which(geneatlas$PV > 1e-6)
rows           <- sample(rows,n)
geneatlas.null <- geneatlas[rows,]

# Combine these into a single table.
geneatlas <- rbind(geneatlas.top,geneatlas.null)

# WRITE RESULTS TO FILE
# ---------------------
cat("Writing association results for selected SNPs to file.\n")
write.csv(geneatlas,out.file,row.names = FALSE,quote = FALSE)

stop()

# -----
# TO DO: Draw a scatterplot
# -----
p1 <- ggplot(geneatlas,aes(x = NBETA/NSE,y = neale_tstat)) +
  geom_point(shape = 20,size = 2) +
  geom_abline(intercept = 0,slope = 1,color = "dodgerblue",
              linetype = "dotted") +
  theme_cowplot() +
  xlim(c(-23,31)) +
  ylim(c(-23,31)) +
  labs(x = "GeneATLAS",y = "Neale lab",title = "t-statistics")

plot_grid(p1,p2)
                           
