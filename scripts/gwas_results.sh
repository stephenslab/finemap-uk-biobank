#!/bin/bash

path='/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/gwas/'

phenotypes=("WBC_count" "RBC_count" "Haemoglobin" "MCV" "RDW" "Platelet_count" "Plateletcrit" "PDW" "Lymphocyte_perc" "Monocyte_perc" "Neutrophill_perc" "Eosinophill_perc" "Basophill_perc" "Reticulocyte_perc" "MSCV" "HLR_perc")

for trait in "${phenotypes[@]}"
do
    for chr in {1..22}
    do
        if [[ ${chr} -eq 1 ]]
        then
          cat ${path}bloodcells_gwas_chr${chr}.${trait}.glm.linear > ${path}bloodcells_gwas_${trait}
        else
          tail -n +2 ${path}bloodcells_gwas_chr${chr}.${trait}.glm.linear >> ${path}bloodcells_gwas_${trait}
        fi
    done
done


