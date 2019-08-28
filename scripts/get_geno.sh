#!/bin/bash

# TO DO: Explain here what this script does, and how to use it.
#
#   qsub -I -l walltime=012:00:00 -l nodes=1:ppn=8 -l mem=40gb
#
# TO DO: How long (roughly) does this script take to complete?
#

# SCRIPT PARAMETERS
# -----------------
# TO DO: Explain what each of these script parameters are for.
PHENO_FILE=/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/height.csv.gz
PHENO_SAMPLE_FILE=ids
PLINK2_EXEC=/gpfs/data/xhe-lab/uk-biobank/tools/plink2
LDSTORE_EXEC=ldstore
BGEN_FILE=/gpfs/data/pierce-lab/uk-biobank-genotypes/ukb_imp_chr3_v3.bgen
GENO_SAMPLE_FILE=/gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample
SNPID_FILE=../data/region-variants-ZBTB38.txt.gz
OUT_FILE=ukb-imp-chr3-height-ZBTB38
CHR=3

# EXTRACT SAMPLE IDs
# ------------------
# Extract the ids of the phenotype samples from the first column of
# the CSV file.
echo "Extracting sample ids from phenotype file."
zcat $PHENO_FILE | tail -n +2 | cut -d "," -f 1 > $PHENO_SAMPLE_FILE

# EXTRACT GENOTYPES
# -----------------
# TO DO: Explain here what this plink command does in more detail.
echo "To do: Explain what this plink command does."
$PLINK2_EXEC --sample $GENO_SAMPLE_FILE --bgen $BGEN_FILE --chr $CHR \
  --keep $PHENO_SAMPLE_FILE --extract $SNPID_FILE --maf 0.01 minor \
  --snps-only --make-pgen --threads 8 --memory 38000 \
  --out $OUT_FILE
 
# plink2 --pfile height.chr3.140800000.141800000.0.01 --threads 8 --export A \
#        --out height.chr3.140800000.141800000.0.01

# gzip height.chr3.140800000.141800000.0.01.raw

# # There are 3,106 variants in the extracted file.

# # COMPUTE MINOR ALLELE FREQUENCIES
# # --------------------------------
# $PLINK2_EXEC --pfile height.chr3.140800000.141800000.0.01 --freq \
#	--out height.chr3.140800000.141800000.0.01

# # Transform pgen file to bgen and compute LD matrix using LDstore:

# # Note that `bgen-1.3` is not supported by LD store.

# $PLINK2_EXEC --pfile height.chr3.140800000.141800000.0.01 --recode bgen-1.2 --out height.chr3.140800000.141800000.0.01

# $LDSTORE_EXEC --bgen height.chr3.140800000.141800000.0.01.bgen --bcor height.chr3.140800000.141800000.0.01.bcor --n-threads 10 --accuracy high --ld-thold 0 --variant-window-size 1000000 --n-variants-chunk 5792
# $LDSTORE_EXEC --bcor height.chr3.140800000.141800000.0.01.bcor --merge 10
# $LDSTORE_EXEC --bcor height.chr3.140800000.141800000.0.01.bcor --matrix height.chr3.140800000.141800000.0.01.matrix

