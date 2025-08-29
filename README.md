# Code and metadata for the analysis of the soil fungi and bacteria communities

This code is part of [INSERT PUBLICATION WHEN DONE] and was created for the analysis of the soil microbial community of different _Acacia mangium_ plantations, compared to reference states of _Imperata cylindrica_ grasslands and remnant forests, all located in the Philippines. Data collection methods are described in [PUBLICATION].

Analyses are organised in separate folders for fungi and bacteria communities, with the related metadata provided. The analyses are conducted using the sequences obtained and processed through DADA2, as described in the publication.
Raw data for both fungi and bacteria can be found in the NCBI Short Read Archive under the BioProject ID and accession number PRJNA1295595.

Below are reported the required packages and useful references for analysis, understanding and results interpretation.

## Required packages to run the analyses:
# For taxonomy
```
library('remotes')
library("dada2")
library('BiocManager')
library('phyloseq')
library('openssl')
library('theseus')
library('tidyverse')
library('Biostrings')
library("ggpubr")
library("dplyr")
library("dendextend")
library('microViz')
library("tidyr")
library("reshape2")
library("gridExtra")
library("scales")
library("parallel")
library("permute")
library("sunburstR")
library("htmltools")
library("lattice")
library("Rmisc")
library("knitr")
library("kableExtra")
library("data.table")
library("grid")
library("microbiome")
library("ape")
library("adespatial")
library("gtable")
library("Biostrings")
library("BAT")
library("car")
library("tibble")
library("cowplot")
library("RVAideMemoire")
library("ggtern")
library("bestNormalize")
library("readxl")
library("ggsci")
library("ggrepel")
library("ggplot2")
library("biohelper")
library("vegan")  
library("agricolae")  
library("eulerr")  
library("stats")  #  For Kruskal-Wallis test
library("viridis")
library("ggforce")
library("reshape")
library("mice")
```
# For microeco analysis
```
library('devtools')
library('GUniFrac') #For phylogenetic tree
library('MiscMetabar') #For phylogenetic tree
library('microeco')
library('file2meco')
library('ggalluvial') #For alluvional plots
library('magrittr') #For PCoA visualisation
library('aplot') #For PCoA visualisation
library('mecoturn')
library ('ggcor')
library('glmmTMB')
library('lmerTest')
library('WGCNA')
library('ggraph')
library("networkD3")
library ('chorddiag')
library('circlize')
library('NetCoMi')
library('SpiecEasi')
library("pheatmap")
library ('ggtree')
library('caret')
library('rfPermute')
library("Boruta")
library("parallel")
library("rsample") 
library("randomForest")  
library("gridExtra")
library("multiROC")
library("ranger")
```
# For other statistical analysis
```
library(metagMisc)
library(rgexf)
library(meconetcomp)
library(indicspecies)
```
