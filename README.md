
# WNV_QTL_Resource
## Description

 This repostitory contains the workflow associated with the manuscript titled'Adaptive immune response to West Nile virus infection in the Collaborative Cross mouse model: A database of cellular phenotypes and Quantitative Trait Loci'. There are four Jupyter notebooks with step by step instructions, along with R programs that contain the code in these notebooks. These processes were used to generate all of the analyses presented in this manuscript.

## SIG_QTL_Mapping_HaplotypeGrouping_Workflow
This is the program/notebook for West Nile virus QTL mapping, haplotype grouping and gene candidate workflow. Thise code only runs one phenotype at a time and this needs to be updated for each time a different phenotype is processed.

## WNV_ QTL_Phenotype_Clustering_Publication
This is the program/notebook for West Nile virus QTL phenotype clustering. This purpose of this step was to expand the criteria for selecting QTL of interest for the hot spot analysis.

## WNV_QTL_Hotspots_Publication
This is the program/notebook for West Nile virus QTL hotspot analysis. This program looks at areas of the genome where there are a larger number of overlapping QTL than expected. Gene candidates in QTL that are also in or near these hot spots are given higher priority for further analysis.

## WNV_QTL_RandomWalk_Prioritization_Publication
This is the program/notebook for West Nile virus QTL random walk gene prioritization. This program creates additional priortization information for gene candidates by evaluating their proximity in a protein protein network to known genes associated with QTL phenotype that contains the candidates.

## WNV_rix_qtl_mapping_functions_publication.R
This is an R program that contains all of the functions necessary to run the above notebooks. These function are documented and most have help functions.



 
