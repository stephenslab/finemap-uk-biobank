library(data.table)
library(dplyr)
pheno_names = c("WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW", "Platelet_count",
                "Plateletcrit", "PDW", "Lymphocyte_perc", "Monocyte_perc",
                "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
                "Reticulocyte_perc", "MSCV", "HLR_perc")


trait_regions = list()
for(trait in pheno_names){
  region = fread(paste0('/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/regions_raw/', trait, '_regions'))
  region = region %>% arrange(desc(logp))
  region_r = c()
  for(i in 1:22){
    region.chr = region %>% filter(CHR == i) %>% arrange(start)
    if(nrow(region.chr) == 0){
      next
    }
    tmp = region.chr %>% group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>%
      summarise(start = first(start), end = max(end), logp = max(logp),.groups = 'drop') %>% 
      mutate(length = end - start) %>%
      mutate(CHR = i) %>% select(CHR, start, end, length, logp)
    region_r = rbind(region_r, tmp)
  }
  trait_regions[[trait]] = region_r
}

saveRDS(trait_regions, '/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/regions/trait_regions.rds')

tb = bind_rows(trait_regions, .id = "column_label")
res.final = c()
for(i in 1:22){
  tb.chr = tb %>% filter(CHR == i) %>% arrange(start)
  if(nrow(tb.chr) == 0){
    next
  }
  tmp = tb.chr %>% group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>%
    summarise(start = first(start), end = max(end), logp = max(logp), .groups = 'drop') %>% 
    mutate(length = end - start) %>%
    mutate(CHR = i) %>% select(CHR, start, end, length, logp)
  res.final = rbind(res.final, tmp)
}

gwas_PDW = fread('/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/gwas/bloodcells_gwas_PDW')
colnames(gwas_PDW)[1] = 'CHR'
gwas_PDW$P = as.numeric(gwas_PDW$P)
gwas_PDW = gwas_PDW %>% mutate(logp = -log10(P))

snpsnum = c()
for(i in 1:nrow(res.final)){
  snpsnum = c(snpsnum, gwas_PDW %>% filter(CHR == res.final$CHR[i],
                                           POS >= res.final$start[i], 
                                           POS <= res.final$end[i]) %>% nrow )
}
res.final$snpsnum = snpsnum

fwrite(res.final, '/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/regions/regions.csv', quote=FALSE,
sep='\t', row.names=FALSE)

