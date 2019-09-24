# finemap-uk-biobank

Code and scripts implementing genetic association and fine-mapping
analyses of UK Biobank data.

## Current fine-mapping analysis pipeline for UK Biobank data

All scripts implementing the data processing and analysis can be found
in the [scripts](scripts) directory. The pipeline is currently
illustrated for fine-mapping standing height in a region on chromosome
3 near gene *ZBTB38*, and can be adapted for other traits. The steps
are as follows:

1. **Prepare phenotype data.** Run R script
   [get_pheno.R](scripts/get_pheno.R) to prepare a CSV file containing
   the phenotype and covariate data from the UK Biobank source
   files. For height, this step creates a new CSV file, `height.csv`,
   containing the phenotype and covariate data.

2. **Prepare SNP data.** Run R script
   [get_geneatlas_snps.R](scripts/get_geneatlas_snps.R) to create a
   table containing independently computed summary statistics. These
   are used to validate our association results. For height, this
   generates a new CSV file, `geneatlas-neale-height.csv`, containing
   the association results. Alternatively, run
   [get_region_snps.R](scripts/get_region_snps.R) to generate a text
   file containing the ids of the genetic variants within the selected
   region, accompanied by independently computed summary statistics,
   when available. This produces a new CSV file,
   `region-variants-ZBTB38.csv`, containing information about the
   selected genetic variants, such as base-pair positions, SNP variant
   ids, and association statistics.

3. **Prepare genotype data.** Run bash script
   [get_geno.sh](scripts/get_geno.sh) ...
