#!/bin/bash


#PBS -N gwas_plink
#PBS -S /bin/bash
#PBS -l mem=50gb
#PBS -l walltime=65:00:00 
#PBS -l nodes=1:ppn=20
#PBS -t 0-21

CHR=$((${PBS_ARRAYID}+1))
phenos=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells.pheno.txt
covar=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells.covar.txt
idfile=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/bloodcells.id.txt
genofile=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/genotypes_maf001_info6/bloodcells_chr${CHR}

# Get variant IDs with information score
# zcat /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_mfi_chr${CHR}_v3.txt.gz | awk '{if ($8 > 0.6) {print $2}}' | sort | uniq > /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/ukb_mfi_chr${CHR}_info0.6.txt 

# Get genotype 
#plink2 --sample /gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample \
#       --bgen /gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr${CHR}_v3.bgen \
#       --chr $CHR \
#       --keep $idfile \
#       --extract /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/ukb_mfi_chr${CHR}_info0.6.txt \
#       --maf 0.001 minor \
#       --snps-only --max-alleles 2 --rm-dup exclude-all \
#       --make-pgen --threads 20 --memory 48000 --out $genofile

# GWAS
plink2 --pfile $genofile --glm hide-covar no-x-sex omit-ref \
       --covar $covar --pheno $phenos --vif 100 \
       --threads 20 --memory 48000 \
       --out /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/gwas_maf001_info6/bloodcells_gwas_chr${CHR}
