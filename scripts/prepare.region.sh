#!/bin/bash

regionchr="$1"
regionfrom="$2"
regionto="$3"
regionname="$4"

zcat /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/height.csv.gz | awk 'BEGIN{FS=","} {print $1 OFS $1}' | tail -n +2 > height.id.txt

plink2 --sample /gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample \
--bgen /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr3_v3.bgen \
--chr ${regionchr} --from-bp ${regionfrom} --to-bp ${regionto} \
--keep height.id.txt \
--maf 0.01 minor \
--snps-only --max-alleles 2 --rm-dup exclude-all --make-pgen --threads 8 --memory 38000 \
--out height.${regionname}.0.01 > height.${regionname}.0.01.plink2.log

plink2 --pfile height.${regionname}.0.01 --threads 8 --export A --out height.${regionname}.0.01

gzip height.${regionname}.0.01.raw

plink2 --pfile height.${regionname}.0.01 --freq --out height.${regionname}.0.01

plink2 --pfile height.${regionname}.0.01 --recode bgen-1.2 --out height.${regionname}.0.01

ldstore --bgen height.${regionname}.0.01.bgen --bcor height.${regionname}.0.01.bcor --n-threads 10 --accuracy high --ld-thold 0 --variant-window-size 1000000 --n-variants-chunk 4000
ldstore --bcor height.${regionname}.0.01.bcor --merge 10
ldstore --bcor height.${regionname}.0.01.bcor --matrix height.${regionname}.0.01.matrix
rm height.${regionname}.0.01.bcor_*

Rscript prepare.susieinput.R ${regionname}