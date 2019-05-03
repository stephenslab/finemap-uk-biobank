# Notes on data sets

Larger data files for this project cannot be committed to the git
repository, so they are being stored separately on the CRI Cluster
("gardner").

## Summary data

+ **Gene ATLAS data.** We downloaded SNP association statistics for
standing height from the [GeneATLAS][gene-atlas] website on May
1, 2019. For full details about the data, see the
[paper][gene-atlas-paper]. These data are currently stored in the
following directory on the CRI cluster:
`/gpfs/data/stephens-lab/finemap-uk-biobank/data/gene-atlas`. There
are four sets of association statistics depending on whether: they
were computed for all SNPs, or only SNPs included on the genotype
array; and whether the phenotype was quantile-normalized or not. These
statistics were computed with an LMM-based association analysis,
controlling for genotype measurement batch, assessment center, age,
and PCs 1--20.

[gene-atlas]:       http://geneatlas.roslin.ed.ac.uk
[gene-atlas-paper]: https://doi.org/10.1038/s41588-018-0248-z

## UK-BioBank Phenotypes

+ The [steps](height.md) to extract height related phenotypes.
