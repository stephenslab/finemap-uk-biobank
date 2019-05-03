# Steps to extract height related phenotypes

The following phenotypes are extracted from UK-BioBank data downloaded on `12-feb-2019`.

The extracted files are stored at `/gpfs/data/stephens-lab/finemap-uk-biobank/`. The whole UK-BioBank phenotypes downloaded on 12-feb-2019 are located at `/gpfs/data/xhe-lab/uk-biobank/data/phenotypes/12-feb-2019`.

  + **height_phenotypes**: ID, sex (31), height (50), UK Biobank assessment centre (54), self-reported ethnic background (21000), age (21022)
  
    1. Get column numbers for the phenotypes:
    ```
    grep -n '31-' ukb26140_header.csv
    grep -n '50-' ukb26140_header.csv
    grep -n '54-' ukb26140_header.csv
    grep -n '21000-' ukb26140_header.csv
    grep -n '21022-' ukb26140_header.csv
    ```
  
    2. Extract the phenotypes:
    ```
    zcat ukb26140.csv.gz | cut -d ',' -f 1,5,37-39,47-49,1156-1158,1171 > /gpfs/data/stephens-lab/finemap-uk-biobank/height_phenotypes.csv
    ```
  
  + **genetic_background**: ID, Genotype measurement batch (22000), Genetic sex (22001), Missingness (22005), Genetic ethnic grouping (22006)
  
    ```
    grep -n '22000-' ukb26140_header.csv
    grep -n '22001-' ukb26140_header.csv
    grep -n '22005-' ukb26140_header.csv
    grep -n '22006-' ukb26140_header.csv
    zcat ukb26140.csv.gz | cut -d ',' -f 1,1172,1173,1176,1177 > /gpfs/data/stephens-lab/finemap-uk-biobank/genetic_background.csv
    ```
  
  + **genetic_relatedness**: ID, Genetic relatedness pairing (22011)
    ```
    grep -n '22011-' ukb26140_header.csv
    zcat ukb26140.csv.gz | cut -d ',' -f 1,1221-1225 > /gpfs/data/stephens-lab/finemap-uk-biobank/genetic_relatedness.csv
    ```
  
  + **genetic_pc**: ID, Genetic principal components (22009)
    ```
    grep -n '22009-' ukb26140_header.csv
    zcat ukb26140.csv.gz | cut -d ',' -f 1,1180-1219 > /gpfs/data/stephens-lab/finemap-uk-biobank/genetic_pc.csv
    ```