### Purpose: This script is used to perform CellChat analysis on the 10k Mouse Brain dataset.
### Please note that this script is a work in progress and is not yet complete. 
### This script is an extension of the previous script "Mouse Brain 10k - Auto Annotations.R".
### Please ensure you understand the previous script before running this script.
### The script is intended to be used as a template for future analysis and may require modifications to suit your specific needs.
### Please set the working directory to your project folder before running this script.

# The working directory should contain the following folders:
# - Data
# - Scripts
# - Results
setwd("")
source("Setup.R")
setup()
source("Mouse Brain 10k Functions.R")

### Downloading Required Data ###
if (!dir.exists("../Data/10k_Adult_Mouse_Brain/"))
{
  download.file("https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Mouse_Brain_CNIK_3p_nextgem_10k_Mouse_Brain_CNIK_3p_nextgem/10k_Mouse_Brain_CNIK_3p_nextgem_10k_Mouse_Brain_CNIK_3p_nextgem_count_sample_filtered_feature_bc_matrix.tar.gz",
                destfile = "../Data/10k_Adult_Mouse_Brain.tar.gz")
  
  untar("../Data/10k_Adult_Mouse_Brain.tar.gz",
        exdir = "../Data/10k_Adult_Mouse_Brain/")
}


raw.data <- Read10X("../Data/10k_Adult_Mouse_Brain/sample_filtered_feature_bc_matrix/")

raw.data <- CreateSeuratObject(raw.data)

raw.data$mt.percent <- PercentageFeatureSet(raw.data, pattern = "^mt-")

raw.data$ribo.percent <- PercentageFeatureSet(raw.data, pattern = "^Rp[sl]")

brain2 <- QC_filtering(raw.data)

brain2 <- brain2 %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.8)

brain2 <- RunUMAP(brain2, dims = 1:20)
num_clusters <- length(unique(brain2$seurat_clusters))

mouse_cells <- celldex::MouseRNAseqData()
scebrain2 <- as.SingleCellExperiment(brain2)

brain2_pred_broad <- SingleR(test = scebrain2,
                            ref = mouse_cells,
                            labels = mouse_cells$label.main)

plotScoreHeatmap(brain2_pred_broad, max.labels = 16,
                 clusters = brain2$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = F) 

# DimPlot(brain2, reduction = "umap", label = T)

brain2_pred_fine <- SingleR(test = scebrain2,
                           ref = mouse_cells,
                           label = mouse_cells$label.fine)

plotScoreHeatmap(brain2_pred_fine, max.labels = 30,
                 clusters = brain2$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = F)

mark_lst <- list(
  'Astrocytes' = c("Aqp4", "Slc1a2", "Gja1", "Aldh1l1")
)

VlnPlot(brain2, features = mark_lst$Astrocytes)

# Only annotating known type by looking at the heatmap.
# Multiple types of cells are sometime listed. This will be resolved later when 
## documentation and gene markers are chosen for testing.
# If the cluster is marked as "unknown", then there are too many cell types expressing within the cluster.

brain2_auto_cluster_ids <- c(
  "0" = "Oligodendrocytes",
  "1" = "Oligodendrocytes",
  "2" = "Microglia",
  "3" = "Neurons",
  "4" = "Neurons & Astrocytes",
  "5" = "Oligodendrocytes",
  "6" = "Neurons",
  "7" = "Neurons",
  "8" = "Neurons",
  "9" = "Astrocytes",
  "10" = "Neurons",
  "11" = "Neurons",
  "12" = "Astrocytes", # IDK
  "13" = "Neurons",
  "14" = "Astrocytes",
  "15" = "Oli & Neur & Astro",
  "16" = "Neurons",
  "17" = "Fibro & Astro",
  "18" = "Neurons",
  "19" = "Neurons",
  "20" = "Glia & Oligo",
  "21" = "Neurons",
  "22" = "Endothelial",
  "23" = "Unknown",
  "24" = "Epithelial",
  "25" = "Fibroblasts",
  "26" = "Neurons",
  "27" = "Neurons",
  "28" = "Neurons"
)

brain2 <- RenameIdents(brain2, brain2_auto_cluster_ids)
brain2$cell_type <- Idents(brain2)
DimPlot(brain2, reduction = "umap", label = T)


#### CellChat Analysis ####
cell_subset <- subset(brain2, idents = c("Endothelial", "Astrocytes"))

cellchat <- createCellChat(object = cell_subset, group.by = "cell_type")

db_cellchat <- CellChatDB.mouse

cellchat@DB <- db_cellchat

# Process & Normalize Data
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

df_net <- subsetCommunication(cellchat, sources.use = "Endothelial", targets.use = "Astrocytes")

# View potential ligand-receptor interactions
head(df_net)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Plot network between endothelial cells and astrocytes
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_network(cellchat, signaling = "VEGF")


