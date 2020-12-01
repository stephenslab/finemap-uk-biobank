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

3. **Prepare genotype data and SuSiE sufficient statistics** Run bash script
   [prepare.region.sh](scripts/prepare.region.sh) to create an RDS file 
   containing sufficient statistics using the genetype data and height. The 
   script requires 4 input: chromosome number, start base-pair position, stop 
   base-pair position, region name. For example,
   ```
   scripts/prepare.region.sh 3 140.8e6 141.8e6 ZBTB38
   ```

## Current analysis pipeline for UK Biobank Blood Cells data

1. **Prepare phenotype data.** Run R script
   [get_bloodcells.R](scripts/get_bloodcells.R) to prepare a CSV file containing
   the phenotype and covariate data from the UK Biobank source
   files.

2. **Run GWAS.** Run [plink_gwas.sh](scripts/plink_gwas.sh) and [gwas_results.sh](scripts/gwas_results.sh) to get GWAS results.

3. **Get fine-mapping regions.** Run R script [get_bloodcells_trait_regions.R](scripts/get_bloodcells_trait_regions.R) to get regions for each trait. Run R script [get_bloodcells_regions.R](scripts/get_bloodcells_regions.R) to combine overlapping regions for each trait and across traits.

4. **Prepare genotype data, LD and z scores for each region.** Run [get_bloodcells_region_genotype_ld.sh](scripts/get_bloodcells_region_genotype_ld.sh) to get genotype data and LD for each region. Run R script [get_bloodcells_zscores.R](scripts/get_bloodcells_zscores.R) to get z scores and XtY for each region.


