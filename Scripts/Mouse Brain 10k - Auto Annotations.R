#### Files Used to Run This Script ####
# Set the working directory to your project folder before running this script.
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

# Reading raw data
raw.data <- Read10X("../Data/10k_Adult_Mouse_Brain/sample_filtered_feature_bc_matrix/")
raw.data <- CreateSeuratObject(raw.data)

# QC
# For brain cells, use the pattern "^mt-" which is lowercase
raw.data$mt.percent <- PercentageFeatureSet(raw.data, pattern = "^mt-")

# For ribosomal genes, use the pattern "^Rp[sl]" which is uppercase
# This is because the ribosomal genes are not in the same format as the mitochondrial genes.
# The ribosomal genes are in the format "Rpl" or "Rps" which is uppercase.
raw.data$ribo.percent <- PercentageFeatureSet(raw.data, pattern = "^Rp[sl]")

# The QC filtering function will filter out cells that have a high percentage of mitochondrial genes and ribosomal genes.
# If you want manually filter the cells, you can use the following code below.
brain2 <- QC_filtering(raw.data)

### Manually filtering the cells:
# VlnPlot(raw.data, features = c("nFeature_RNA", "nCount_RNA", "mt.percent", "ribo.percent"), ncol = 4)
# 
# FeatureScatter(raw.data, "nCount_RNA", "nFeature_RNA") +
#   geom_abline(slope = 0, intercept = 6500, col = "blue") +
#   geom_abline(slope = 0, intercept = 500, col = "red") 
# 
# FeatureScatter(raw.data, "nCount_RNA", "mt.percent") +
#   geom_abline(slope = 0, intercept = 10, col = "blue")
# 
# FeatureScatter(raw.data, "nCount_RNA", "ribo.percent") +
#   geom_abline(slope = 0, intercept = 5.25, col = "blue")
# 
# brain <- subset(raw.data, nFeature_RNA > 500 &
#                   nFeature_RNA < 6500 &
#                   mt.percent < 10 &
#                   ribo.percent < 5.25)

# Normalize, find variable features, scale data & Run PCA (Principal Component Analysis)
brain <- brain %>% 
  NormalizeData() %>%
  FindVariableFeatures()


var_feat_top10 <- head(VariableFeatures(brain), 10)
var_feat <- VariableFeaturePlot(brain)
var_feat <- LabelPoints(plot = var_feat, points = var_feat_top10, repel = T)

brain <- brain %>%
  ScaleData() %>%
  RunPCA()

DimPlot(brain, reduction = "pca")

VizDimLoadings(brain, dims = 1:2)

# PCA Comp 1
FeaturePlot(brain, features = c("Rbfox3", "Nrg3"), reduction = "pca")

ElbowPlot(brain, ndims = 50) +
  geom_abline(slope = 0, intercept = 2.2, col = "red")

# Find Nearest Neighbors, Cluster & run UMAP
brain <- brain %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.8)

brain <- RunUMAP(brain, dims = 1:20)

DimPlot(brain, reduction = "umap", label = T)

num_clusters <- length(unique(brain$seurat_clusters))


#### Auto Annotations ####
mouse_cells <- celldex::MouseRNAseqData()
sceBrain <- as.SingleCellExperiment(brain)

Brain_pred_broad <- SingleR(test = sceBrain,
                      ref = mouse_cells,
                      labels = mouse_cells$label.main)

plotScoreHeatmap(Brain_pred_broad, max.labels = 16,
                 clusters = brain$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = F) 

# DimPlot(brain, reduction = "umap", label = T)

Brain_pred_fine <- SingleR(test = sceBrain,
                           ref = mouse_cells,
                           label = mouse_cells$label.fine)

plotScoreHeatmap(Brain_pred_fine, max.labels = 30,
                 clusters = brain$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = F)

mark_lst <- list(
  'Astrocytes' = c("Aqp4", "Slc1a2", "Gja1", "Aldh1l1")
)

VlnPlot(brain, features = mark_lst$Astrocytes)

# Only annotating known type by looking at the heatmap.
# Multiple types of cells are sometime listed. This will be resolved later when 
## documentation and gene markers are chosen for testing.
# If the cluster is marked as "unknown", then there are too many cell types expressing within the cluster.

auto_cluster_ids <- c(
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

brain <- RenameIdents(brain, auto_cluster_ids)
brain$cell_type <- Idents(brain)
DimPlot(brain, reduction = "umap", label = T)


#### All UMAP Analyses ####
Idents(brain) <- brain$seurat_clusters

##### Clusters #####
all_markers <- FindAllMarkers(brain, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)

all_clusters <- cluster_markers(all_markers, seq(0, num_clusters - 1), DEG_filter = F, expression_filter = T)
clus9_matrix <- all_clusters$Cluster_9$Marker_Matrix
clus9_genes <- all_clusters$Cluster_9$Marker_Genes
6
clus14_matrix <- all_clusters$Cluster_14$Marker_Matrix
clus14_genes <- all_clusters$Cluster_14$Marker_Genes


#### Cluster 9 & Cluster 14 Comparison ####
##### Difference of Cluster 9 (Base) to Cluster 14 (Cluster to check) #####

## This is a list of genes that were present in cluster 14 
### than were not present in cluster 9, which was also ordered by the avg_log2FC
### where higher number represent higher expression within the cluster.
comparing_res <- comparative_genes(clus9_genes, clus14_genes)
diff_genes <- comparing_res$Diff14
diff_genes2 <- comparing_res$Diff9
same_genes <- comparing_res$Same

VlnPlot(brain, features = head(diff_genes))
VlnPlot(brain, features = head(diff_genes2))
VlnPlot(brain, features = head(same_genes))


##### EnrichR #####
dbs <- listEnrichrDbs()
to_check <- c("Human_Gene_Atlas", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
              "KEGG_2019_Human", "MSigDB_Hallmark_2020")

enr_clus9 <- enrichr(clus9.genes, databases = to_check)

# Get the top hits from each and plot them
enrichr_analysis(enr_clus9, "9")


