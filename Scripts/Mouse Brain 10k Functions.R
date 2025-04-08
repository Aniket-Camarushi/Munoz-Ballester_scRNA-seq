##### This file encodes the functions used for "Mouse Brain 10k - Auto Annotations"

qc_plots <- function(data, feat_min = NULL, feat_max = NULL, mito_max = NULL, ribo_max = NULL)
{
  vln <- VlnPlot(raw.data, features = c("nFeature_RNA", "nCount_RNA", "mt.percent", "ribo.percent"), ncol = 4)
  
  if (is.null(feat_min) | is.null(feat_max) | is.null(mito_max) | is.null(ribo_count))
  {
    feature_count <- FeatureScatter(raw.data, "nCount_RNA", "nFeature_RNA") +
      theme(legend.position = "none") +  
      ggtitle("Count vs Feature") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    mt_count <- FeatureScatter(raw.data, "nCount_RNA", "mt.percent") +
      theme(legend.position = "none") +  
      ggtitle("Count vs Mito %") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    ribo_count <- FeatureScatter(raw.data, "nCount_RNA", "ribo.percent") +
      theme(legend.position = "none") +  
      ggtitle("Count vs Ribo %") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  }
  else
  {
    feature_count <- FeatureScatter(raw.data, "nCount_RNA", "nFeature_RNA") +
      geom_abline(slope = 0, intercept = feat_max, col = "blue") +
      geom_abline(slope = 0, intercept = feat_min, col = "red") +
      theme(legend.position = "none") +  
      ggtitle("Count vs Feature") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    mt_count <- FeatureScatter(raw.data, "nCount_RNA", "mt.percent") +
      geom_abline(slope = 0, intercept = mito_max, col = "blue") +
      theme(legend.position = "none") +  
      ggtitle("Count vs Mito %") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    ribo_count <- FeatureScatter(raw.data, "nCount_RNA", "ribo.percent") +
      geom_abline(slope = 0, intercept = ribo_max, col = "blue") +
      theme(legend.position = "none") +  
      ggtitle("Count vs Ribo %") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  }
  
  final_plot <- (vln / (feature_count | (mt_count | ribo_count))) + plot_layout(nrow = 2)
  
  print(final_plot)
}

QC_filtering <- function(data)
{
  qc_choice <- as.numeric(readline(
    prompt = "Would you like to visualize QC cutoffs?\n\t1. Yes\n\t2. No\nEnter Choice: "))
  
  if (!is.na(qc_choice) & qc_choice == 1)
  {
    qc_plots(data)
    
    subsetting <- as.numeric(readline(
      prompt = "\nWould you like to subset cutoff parameters?\n\t1. Yes\n\t2. No\nEnter Choice: "))
    
    if (!is.na(subsetting) & subsetting == 1)
    {
      subset_choice <- 0
      
      while (subset_choice != 1) 
      {
        # Get user input for filtering conditions
        feat_min <- as.numeric(readline(prompt = "\nEnter minimum nFeature_RNA: "))
        feat_max <- as.numeric(readline(prompt = "Enter maximum nFeature_RNA: "))
        
        mito_max <- as.numeric(readline(prompt = "Enter maximum mito.percent: "))
        
        ribo_max <- as.numeric(readline(prompt = "Enter maximum ribo.percent: "))
        
        qc_plots(data, feat_min, feat_max, mito_max, ribo_max)
        
        subset_choice <- as.numeric(readline(
          prompt = "\nAre you sure about these QC cutoffs?\n\t1. Yes\n\t2. No\nEnter Choice: "))
      }
      
      if (!is.na(subset_choice) & subset_choice == 1)
      {
        cat("Subsetting data...\nDone!")
        brain <- subset(raw.data, nFeature_RNA > feat_min &
                          nFeature_RNA < feat_max &
                          mt.percent < mito_max &
                          ribo.percent < ribo_max)
        return(brain)
      }
      else
      {
        glue("\nYour last parameters were:
    nFeature_RNA > {feat_min} & nFeature_RNA < {feat_max}
    mito.percent < {mito_max}
    mito.percent < {ribo_max}")
      }
    }
    else
    {
      cat("Manual Subsetting!!")
    }
  }
  else
  {
    cat("Manual QC!!")
  }
}

HighDEG_Filter <- function(cluster_markers)
{
  genes <- c()
  
  for (i in VariableFeatures(brain))
  {
    if (i %in% cluster_markers)
    {
      genes <- c(genes, i)
    }
  }
  
  genes <- cluster_markers[order(cluster_markers$avg_log2FC, decreasing = T), ]
  
  return(genes)
}

specific_cluster_markers <- function(all_markers, clus_num = "0", DEG_filter = F, expression_filter = T)
{
  cluster.markers <- all_markers %>%
    filter(p_val_adj < 0.05) %>%
    filter(cluster == clus_num)
  
  cluster.markers.genes <- pull(cluster.markers, gene)
  
  if (DEG_filter) cluster.markers <- HighDEG_Filter(cluster.markers) 
  
  if (expression_filter) cluster.markers.genes <- cluster.markers.genes[
    order(cluster.markers$avg_log2FC, decreasing = T)]
  
  return(list("Marker_Matrix" = cluster.markers, "Marker_Genes" = cluster.markers.genes))
}

cluster_markers <- function(all_markers, cluster_numbers, DEG_filter = F, expression_filter = T)
{
  cluster_numbers <- as.character(cluster_numbers)
  
  all_clusters <- list()
  
  for (clus in cluster_numbers) 
  {
    res <- specific_cluster_markers(all_markers, clus_num = clus, 
                                    DEG_filter = DEG_filter, expression_filter = expression_filter)
    
    all_clusters[[glue("Cluster_{clus}")]] <- res
  }
  
  return(all_clusters)
}

comparative_genes <- function(base_genes, genes_to_check) 
{
  same_genes <- intersect(base_genes, genes_to_check)
  diff_genes <- setdiff(genes_to_check, base_genes)
  diff_genes2 <- setdiff(base_genes, genes_to_check)
  
  # diff_genes <- diff_genes[order(genes_to_check$avg_log2FC, decreasing = T), ]
  # same_genes <- same_genes[order(base_genes$avg_log2FC, decreasing = T), ]
  # 
  # return(c(diff_genes, same_genes))
  
  # Filter the rows for diff_genes and same_genes
  diff_genes_data <- genes_to_check[genes_to_check %in% diff_genes]
  diff_genes2_data <- base_genes[base_genes %in% diff_genes2]
  same_genes_data <- base_genes[base_genes %in% same_genes]
  
  # Sort both diff_genes and same_genes by avg_log2FC in decreasing order
  # diff_genes_data <- diff_genes_data[order(diff_genes_data$avg_log2FC, decreasing = TRUE), ]
  # diff_genes_data2 <- diff_genes2_data[order(diff_genes2_data$avg_log2FC, decreasing = TRUE), ]
  # same_genes_data <- same_genes_data[order(same_genes_data$avg_log2FC, decreasing = TRUE), ]
  
  res <- list("Diff14" = diff_genes, "Diff9" = diff_genes2, "Same" = same_genes)
  
  return(res)
}

comparative_genes_matrix <- function(base_cluster, cluster_to_check) 
{
  same_genes <- intersect(base_cluster$gene, cluster_to_check$gene)
  diff_genes <- setdiff(cluster_to_check$gene, base_cluster$gene)
  diff_genes2 <- setdiff(base_cluster$gene, cluster_to_check$gene)
  
  # diff_genes <- diff_genes[order(cluster_to_check$avg_log2FC, decreasing = T), ]
  # same_genes <- same_genes[order(base_cluster$avg_log2FC, decreasing = T), ]
  # 
  # return(c(diff_genes, same_genes))
  
  # Filter the rows for diff_genes and same_genes
  diff_genes_data <- cluster_to_check[cluster_to_check$gene %in% diff_genes, ]
  diff_genes2_data <- base_cluster[base_cluster$gene %in% diff_genes2, ]
  same_genes_data <- base_cluster[base_cluster$gene %in% same_genes, ]
  
  # Sort both diff_genes and same_genes by avg_log2FC in decreasing order
  diff_genes_data <- diff_genes_data[order(diff_genes_data$avg_log2FC, decreasing = TRUE), ]
  diff_genes_data2 <- diff_genes2_data[order(diff_genes2_data$avg_log2FC, decreasing = TRUE), ]
  same_genes_data <- same_genes_data[order(same_genes_data$avg_log2FC, decreasing = TRUE), ]
  
  res <- list("Diff14" = diff_genes_data, "Diff9" = diff_genes2_data, "Same" = same_genes_data)
  
  return(res)
}

plot_eres <- function(eres_name, eres_list, n = 10) 
{
  eres <- eres_list[[eres_name]]
  eres %>%
    top_n(n = n, wt = -log10(Adjusted.P.value)) %>%
    arrange(-log10(Adjusted.P.value)) %>%
    mutate(Term = factor(Term, levels = Term)) %>%
    ggplot(mapping = aes(x = Term, y = -log10(Adjusted.P.value), fill = Combined.Score)) +
    geom_bar(stat = "identity") +
    ggpubr::rotate() +
    theme_bw(base_size = 16) +
    rremove("ylab") +
    labs(title = eres_name)
}

enrichr_analysis <- function(func_list, clus_num)
{
  plotList <- lapply(names(func_list), plot_eres, eres_list = func_list)
  cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                     label_size = 18, align = "vh") +
    patchwork::plot_annotation(title = glue("Cluster {clus_num} Enrichr Analysis"), 
                               theme = theme(title = element_text(size = 20)))
}


