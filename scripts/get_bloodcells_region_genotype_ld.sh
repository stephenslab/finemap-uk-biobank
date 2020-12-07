#!/bin/bash

regionchr="$1"
regionfrom="$2"
regionto="$3"

module load gcc/6.2.0 R/3.5.0

outputfile=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/regions_genotype_maf001_info6/bloodcells_chr${regionchr}.${regionfrom}.${regionto}

ldfile=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/regions_ld_maf001_info6/bloodcells_chr${regionchr}.${regionfrom}.${regionto}

plink2 --pfile /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/genotypes_maf001_info6/bloodcells_chr${regionchr} \
--chr ${regionchr} --from-bp ${regionfrom} --to-bp ${regionto} \
--make-pgen --threads 20 --memory 48000 --out ${outputfile} > ${outputfile}.plink2.log

plink2 --pfile ${outputfile} --threads 10 --memory 48000 --export A --out ${outputfile}

gzip ${outputfile}.raw

plink2 --pfile ${outputfile} --threads 10 --memory 48000 --freq --out ${outputfile}

plink2 --pfile ${outputfile} --threads 10 --memory 48000 --export bgen-1.2 --out ${outputfile}

ldstore --bgen ${outputfile}.bgen --bcor ${ldfile}.bcor --n-threads 10 --accuracy high --ld-thold 0 
ldstore --bcor ${ldfile}.bcor --merge 10
ldstore --bcor ${ldfile}.bcor --matrix ${ldfile}.matrix
rm ${ldfile}.bcor_*

