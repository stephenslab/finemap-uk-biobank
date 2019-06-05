# This short script illustrates how to read genotype data from a
# ".traw" file (for details on this format, see
# https://www.cog-genomics.org/plink/2.0/formats#traw).
library(data.table)
geno.file <- "nod2.traw.gz"

# Load the data from the .traw file.
geno <- fread(geno.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(geno) <- "data.frame"

# Extract the SNP information.
map <- geno[1:6]

# Extract the genotypes.
geno <- geno[-(1:6)]
geno <- as.matrix(geno)
