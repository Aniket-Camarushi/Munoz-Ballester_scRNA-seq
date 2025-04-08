### This is a setup file to load all necessary packages for this project ###

setup <- function()
{
  pkgs_list <- c("tidyverse", "hdf5r", "ggpubr", "Seurat", "enrichR", "cowplot", "patchwork", "VennDiagram")
  install.packages(pkgs_list[! pkgs_list %in% rownames(installed.packages())])
  
  bioc_pkgs <- c("limma", "org.Hs.eg.db", "celldex", "SingleR", "SingleCellExperiment", "viridis")
  bioc2install <- bioc_pkgs[! bioc_pkgs %in% rownames(installed.packages())]
  
  if (length(bioc2install) > 0) 
  {
    BiocManager::install(bioc2install)
  }
  
  
  # Load packages
  library(Seurat)
  library(hdf5r)
  library(tidyverse)
  library(ggpubr)
  library(limma)
  library(celldex)
  library(SingleR)
  library(enrichR)
  library(cowplot)
  library(patchwork)
  library(VennDiagram)
  library(org.Hs.eg.db)
  library(SingleCellExperiment)
  library(celldex)
  library(viridis)
  library(glue)
  library(dplyr)
  library(rafalib)
  library(CellChat)
  
  rm(pkgs_list, bioc_pkgs, bioc2install)
}
