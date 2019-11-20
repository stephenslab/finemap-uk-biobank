#!/bin/bash

regionchr="$1"
regionfrom="$2"
regionto="$3"
regionname="$4"

module load gcc/6.2.0 R/3.5.0

zcat /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/height.csv.gz | awk 'BEGIN{FS=","} {print $1 OFS $1}' | tail -n +2 > height.id.txt

zcat /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_mfi_chr${regionchr}_v3.txt.gz | awk '($8 < 0.9){print $2}' > ukb_mfi_chr${regionchr}_v3_exclude_id.txt

plink2 --sample /gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample \
--bgen /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr${regionchr}_v3.bgen \
--chr ${regionchr} --from-bp ${regionfrom} --to-bp ${regionto} \
--keep height.id.txt \
--exclude ukb_mfi_chr${regionchr}_v3_exclude_id.txt \
--maf 0.01 minor \
--snps-only --max-alleles 2 --rm-dup exclude-all --make-pgen --threads 8 --memory 38000 \
--out height.${regionname} > height.${regionname}.plink2.log

plink2 --pfile height.${regionname} --threads 8 --export A --out height.${regionname}

gzip height.${regionname}.raw

plink2 --pfile height.${regionname} --freq --out height.${regionname}

plink2 --pfile height.${regionname} --recode bgen-1.2 --out height.${regionname}

ldstore --bgen height.${regionname}.bgen --bcor height.${regionname}.bcor --n-threads 10 --accuracy high --ld-thold 0
ldstore --bcor height.${regionname}.bcor --merge 10
ldstore --bcor height.${regionname}.bcor --matrix height.${regionname}.matrix
rm height.${regionname}.bcor_*

Rscript prepare.susieinput.R ${regionname}
