#!/bin/bash

#PBS -N bloodcells
#PBS -S /bin/bash
#PBS -l mem=50gb
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=10
#PBS -t 0-971
module load gcc/6.2.0 R/3.5.0

cd /home/yuxin/GitHub/finemap-uk-biobank/scripts

linenum=$((${PBS_ARRAYID}+2))
echo ${linenum}
line=$(cat /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/regions/regions.csv | head -n ${linenum} | tail -n 1)
echo ${line}

regionchr=$(echo ${line} | awk '{print $1;}')
regionfrom=$(echo ${line} | awk '{print $2;}')
regionto=$(echo ${line} | awk '{print $3;}')

./get_bloodcells_region_genotype_ld.sh ${regionchr} ${regionfrom} ${regionto}

# Rscript get_bloodcells_zscores.R ${regionchr} ${regionfrom} ${regionto}

