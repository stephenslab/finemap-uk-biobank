#!/bin/bash

# This will run a basic association analysis in PLINK. See the "data"
# folder for the input files. To interpret output, see:
# https://www.cog-genomics.org/plink/1.9/formats#assoc_linear
plink --linear --allow-no-sex --output-missing-phenotype NA \
   --bfile 1kg --covar 1kg.covar --out demo
