arg = commandArgs(trailingOnly = TRUE)

chr = as.numeric(arg[1])
start = as.numeric(arg[2])
end = as.numeric(arg[3])

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))

gwas_dir = '/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/gwas_maf001_info6/'
genotype_dir = '/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/regions_genotype_maf001_info6/'
output_dir = '/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/regions_zscores_maf001_info6/'

pheno_names = c("WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW", "Platelet_count",
                "Plateletcrit", "PDW", "Lymphocyte_perc", "Monocyte_perc",
                "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
                "Reticulocyte_perc", "MSCV", "HLR_perc")

zscores = c()
for(trait in pheno_names){
  gwas = fread(paste0(gwas_dir,'bloodcells_gwas_chr', chr, '.', trait, '.glm.linear'))
  colnames(gwas)[1] = 'CHR'
  dt = gwas %>% filter(CHR == chr, POS >= start, POS <= end) %>% select(CHR, POS, ID, REF, ALT, T_STAT) %>% 
    rename(!!trait := T_STAT)
  if(is.null(nrow(zscores))){
    zscores = dt
  }else{
    zscores = inner_join(zscores, dt, by=c('ID', 'CHR', 'POS', 'REF', 'ALT'))
  }   
}

pos = zscores %>% select(CHR, POS, ID, REF, ALT)
zscores = zscores %>% select(-CHR, -POS, -ID, -REF, -ALT)

geno.file = paste0(genotype_dir, 'bloodcells_chr', 
                   chr, '.', start, '.', end, '.raw.gz')

cat("Reading genotype data.\n")
geno <- suppressMessages(fread(geno.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE))
class(geno) <- "data.frame"
# Extract the genotypes.
X <- as(as.matrix(geno[-(1:6)]),'dgCMatrix')
# X <- as.matrix(geno[-(1:6)])

pheno.file <- "/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells.pheno.txt"
cat("Reading phenotype data.\n")
pheno        <- suppressMessages(fread(pheno.file))
class(pheno) <- "data.frame"

ind = fread(paste0(genotype_dir, 'bloodcells_chr', 
                   chr, '.', start, '.', end, '.psam'))
match.idx = match(ind$IID, pheno$IID)
pheno = pheno[match.idx,]

Y = pheno %>% select(-FID, -IID) %>% as.matrix

# centering
Y = t(t(Y) - colMeans(Y))
XtY = as.matrix(Matrix::crossprod(X, Y))

X.cm = colMeans(X)
xtxdiag = colSums(X^2) - nrow(X) * X.cm^2

maf <- suppressMessages(read.delim(paste0(genotype_dir, 'bloodcells_chr', 
                         chr, '.', start, '.', end, '.afreq')))
maf = maf %>% mutate(maf = pmin(ALT_FREQS, 1-ALT_FREQS)) %>% select(ID, maf)
pos = inner_join(pos, maf, by='ID')

saveRDS(list(LD = paste0('bloodcells_chr', chr, '.', start, '.', end,'.matrix'),
             XtY = XtY, n = nrow(Y), pos=pos, Z = zscores, XtXD = xtxdiag), 
        paste0(output_dir, 'bloodcells_chr',chr, '.', start, '.', end, '.z.rds'))


