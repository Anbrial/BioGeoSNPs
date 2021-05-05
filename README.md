# BioGeoSNPs
This repository includes the R code that was used in Llanos-Garrido et al. (2021) Environmental association model with loci under divergent selection predicts the distribution range of a lizard. DOI: TBD.

The genomic data required to fully replicate the analysis of this paper is in another repository (PANGAEA): https://doi.org/10.1594/PANGAEA.908220. We included the environmental PCA scores and 21 outlier genotypes in the 00_input folder.

The 01_randomized_GEAM folder contains the code that is needed to fully replicate the genotype-environment association analysis (GEAM) that is described in the paper. See the R script for further explanation on this analysis.

The 02_genotype_per_cell (work in progress) will include the resulting model in an QGIS-readable shapefile that will store the inferred genotypes presence at each 1x1 km geographic cell within the inferred species' range.
