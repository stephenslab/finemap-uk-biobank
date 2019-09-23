library(data.table)

markers <- read.table("../data/region-variants-ZBTB38.txt.gz",
                      stringsAsFactors = FALSE)[[1]]

map <- fread(file.path("/gpfs/data/pierce-lab/uk-biobank-genotypes",
                       "ukb_imp_chr3_v3.bim.gz"),sep = "\t",header = FALSE)
class(map) <- "data.frame"
names(map) <- c("chr","id","cM","pos","A1","A2")

geneatlas <- fread(file.path("/gpfs/data/stephens-lab/finemap-uk-biobank",
                   "data/summary/gene-atlas/height/imputed/not-normalized",
                   "imputed.allWhites.50-0.0.chr3.csv.gz"),header = TRUE)
class(geneatlas) <- "data.frame"

geneatlas <- subset(geneatlas,is.element(SNP,markers))
rows      <- match(geneatlas$SNP,map$id)
geneatlas <- cbind(map[rows,],geneatlas)
rownames(geneatlas) <- NULL
geneatlas <- geneatlas[c("chr","id","pos","PV-50-0.0")]
names(geneatlas)[4] <- "pvalue"

write.csv(geneatlas,"sample_locuszoom_data.csv",row.names = FALSE,
          quote = FALSE)
