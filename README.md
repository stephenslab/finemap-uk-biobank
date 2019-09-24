# finemap-uk-biobank

Code and scripts implementing genetic association and fine-mapping
analyses of UK Biobank data.

## Current fine-mapping analysis pipeline for UK Biobank data

All scripts implementing the data processing and analysis can be found
in the [scripts](scripts) directory. The pipeline is currently
illustrated for fine-mapping standing height in a region on chromosome
3 near gene *ZBTB38*, and can be adapted for other traits. The steps
are as follows:

1. Run R script [get_pheno.R](scripts/get_pheno.R) to prepare a CSV
   file containing the phenotype and covariate data from the UK
   Biobank source files. For height, this step creates a new file,
   `height.csv`, containing the phenotype and covariate data.

2. Run R script [get_geneatlas_snps.R](scripts/get_geneatlas_snps.R)
   to create a table containing independently computed summary
   statistics. These are used to validate our association results. For
   height, this generates a file, `geneatlas-neale-height.csv`,
   containing the association results.

3. Optionally, run [get_snps.R](scripts/get_snps.R) to generate a text
   file containing the ids of the genetic variants within the selected
   region. (TO DO: Add association statistics to the output file.)

5. Run bash script [get_geno.sh](scripts/get_geno.sh) ...
