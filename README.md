# finemap-uk-biobank

Code and scripts implementing genetic association and fine-mapping
analyses of UK Biobank data.

## Current fine-mapping analysis pipeline for UK Biobank data

All scripts implementing the data processing and analysis can be found
in the [scripts](scripts) directory. The pipeline is currently
illustrated for height, and can be adapted for other traits. The steps
are as follows:

1. Run [get_pheno.R](scripts/get_pheno.R) to prepare a CSV file
   containing the phenotype and covariate data from the UK Biobank
   source files.

2. Run [get_geneatlas_snps.R](scripts/get_geneatlas_snps.R) to create
   a table containing independently computed summary statistics. These
   are used to validate our association results.

3. *Add third step here.*
