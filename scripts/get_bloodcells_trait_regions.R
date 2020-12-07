library(dplyr)
library(data.table)
pheno_names = c("WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW", "Platelet_count",
                "Plateletcrit", "PDW", "Lymphocyte_perc", "Monocyte_perc",
                "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
                "Reticulocyte_perc", "MSCV", "HLR_perc")

for( trait in pheno_names){
  dat.file <-
    file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/gwas_maf001_info6/",
              paste0("bloodcells_gwas_", trait))
  out.file <- file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/regions_raw/",
                        paste0(trait, '_regions'))
  cat(paste0("Loading association results, ", trait, ".\n"))
  dat <- fread(dat.file)
  class(dat) <- "data.frame"
  dat        <- dat %>% select('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'BETA', 'SE', 'T_STAT', 'P')
  colnames(dat)[1] = 'CHR'
  dat$P = as.numeric(dat$P)

  # filter on chr: delete chr 6 MHC region
  dat.na = dat %>% filter(! (CHR == 6 & POS>=25000000 & POS<=36000000))
  dat.na$logp = -log10(dat.na$P)
  chrpos = dat.na %>% group_by(CHR) %>% summarise(start = min(POS), end = max(POS))

  res = matrix(NA, 0, 5)

  i = 0
  while(max(dat.na$logp, na.rm=T) > -log10(5e-8)){
    i = i + 1
    signal = dat.na %>% filter(logp == max(logp, na.rm=T)) %>% top_n(1, POS) %>% select(CHR, POS, logp) %>%
      mutate(start = max(POS - 250000, chrpos$start[which(chrpos$CHR == CHR)]),
             end = min(POS + 250000, chrpos$end[which(chrpos$CHR == CHR)]))
    res = rbind(res, signal)
    dat.na[which(dat.na$CHR == signal$CHR & dat.na$POS <= signal$end & dat.na$POS >= signal$start),] = NA
  }

  fwrite(res, out.file, quote = FALSE, col.names=T, row.names = FALSE, sep = "\t")
}

