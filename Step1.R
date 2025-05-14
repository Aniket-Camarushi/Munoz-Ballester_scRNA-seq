#load packages----

library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(edgeR)
library(matrixStats)
library(cowplot)
library(DT)
library(gt)
library(plotly)
library(limma) # venerable package for differential gene expression using linear modeling
library(RColorBrewer) #need colors to make heatmaps
library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps
library(d3heatmap) # install from github with: devtools::install_github("talgalili/d3heatmap")
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(enrichplot) # great for making the standard GSEA enrichment plots
library(msigdbr)
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(dplyr)
#update.packages(ask = FALSE, checkBuilt = TRUE)


targets <- read_tsv("studydesign.txt")# read in your study design


#Step2_dataWrangling----
sampleLabels <- targets$sample
XXX <- factor(targets$XXX)              #XXX and YYY are in this case proteins with different levels (KO,WT,mutant)
YYY <- factor(targets$YYY)
treatment <- factor(targets$treatment)  #treatment is for example, radiation
group <- factor(targets$type)

#Saved DGEList objects can be easily shared and loaded into an R environment
load(file = "myDGEList")
#myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = -1, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

