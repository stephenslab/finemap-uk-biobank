#!/bin/bash


#PBS -N gwas_plink
#PBS -S /bin/bash
#PBS -l mem=60gb
#PBS -l walltime=65:00:00 
#PBS -l nodes=1:ppn=20
#PBS -t 0-21

CHR=$((${PBS_ARRAYID}+1))
phenos=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells.pheno.txt
covar=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells.covar.txt
idfile=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells.id.txt
genofile=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/genotypes/bloodcells_chr${CHR}

# Get variant IDs with information score
zcat /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_mfi_chr${CHR}_v3.txt.gz | awk '{if ($8 != "NA" && $8 > 0.9) {print $2}}' | sort | uniq > /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/ukb_mfi_chr${CHR}_info0.9.txt 

# Get genotype 
plink2 --sample /gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample \
       --bgen /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr${CHR}_v3.bgen \
       --chr $CHR \
       --keep $idfile \
       --extract /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/ukb_mfi_chr${CHR}_info0.9.txt \
       --maf 0.01 minor \
       --snps-only --max-alleles 2 --rm-dup exclude-all \
       --make-pgen --out $genofile

# GWAS
plink2 --pfile $genofile --glm hide-covar no-x-sex omit-ref \
       --covar $covar --pheno $phenos --vif 100 \
       --out /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/gwas/bloodcells_gwas_chr${CHR}
