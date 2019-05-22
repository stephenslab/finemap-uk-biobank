## Filter data
library(dplyr)

dat.file = '/gpfs/data/stephens-lab/finemap-uk-biobank/height.csv'
out.file = '/gpfs/data/stephens-lab/finemap-uk-biobank/height.filter.csv'
### Read data
dat = read.csv(dat.file, header = T)

### Remove individual with sex mismatched between self-reported sex and genetic sex
### This removes 310 individuals
dat = dat %>% filter(sex == sex_genetic)

### Remove individual with high missing rate
### This removes 198 individuals
dat = dat %>% filter(missingness < 0.05)

### Remove outliers defined by UK BioBank when we have the data

### Remove related individuals
dat = dat %>% filter(is.na(relatedness_genetic0))

### Remove individual with abnormal height
### remove individuals with height departing 3 standard deviations from their gender median
dat.m = dat %>% filter(sex == 1) %>% filter(height < (median(height) + 3*sd(height)) & height > (median(height) - 3*sd(height)))
dat.f = dat %>% filter(sex == 0) %>% filter(height < (median(height) + 3*sd(height)) & height > (median(height) - 3*sd(height)))
dat = rbind(dat.m, dat.f)

# SUMMARIZE DATA
# --------------
# Double-check that everything looks okay.
summary(dat)

cols <- c(which(grepl("relatedness_genetic",names(dat))), which(grepl('sex_genetic', names(dat))), which(grepl('ethic', names(dat))))

dat = dat[,-cols]

# WRITE DATA TO FILE
# ------------------
write.csv(dat,out.file,row.names = FALSE,quote = FALSE)
