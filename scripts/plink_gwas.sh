#!/bin/bash

# We request a compute nodes with 10 CPUs and 3 GB of memory, and at
# most 10 minutes of runtime.

#PBS -N gwas_plink
#PBS -S /bin/bash
#PBS -l mem=3gb
#PBS -l walltime=05:00:00 
#PBS -l nodes=1:ppn=10

CHR=1
phenos=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells.pheno.txt
covar=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells.covar.txt

# Get variant IDs with information score
zcat /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_mfi_chr${CHR}_v3.txt.gz | awk '{if ($8 != "NA" && $8 > 0.9) {print $1}}' > /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/ukb_mfi_chr${CHR}_info0.9.txt 

plink2 --sample /gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample \ 
--bgen /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr${CHR}_v3.bgen \
--chr $CHR \ 
--maf 0.01 minor \
--extract /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/ukb_mfi_chr${CHR}_info0.9.txt
--linear hide-covar no-x-sex omit-ref \
--covar $covar --pheno $phenos \ 
--vif 100
--out /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells_gwas_chr${CHR}
