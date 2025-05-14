# looking into individual GENES differential expression

#normalized expression matrix
expr_matrix <- v.DEGList.filtered.norm$E
head(rownames(expr_matrix))

#load the reactome
reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, gene_symbol)

# Extract the pathways
pathway1 <- reactome %>% filter(term == "name of reactome pathway1 from myGSEA.df")
pathway2 <- reactome %>% filter(term == "name of reactome pathway2 from myGSEA.df")

# Combine gene lists and remove duplicates
genes_of_interest <- unique(c(pathway1$gene, pathway2$gene))

#filter expression matrix to only those genes
expr_subset <- expr_matrix[rownames(expr_matrix) %in% genes_of_interest, ]

# Check dimensions
dim(expr_subset)

#### MAKE A HEATMAP OF EXPRESSION ACROSS SAMPLES
# Optional: Scale by row (gene)
expr_scaled <- t(scale(t(expr_subset)))  # z-score per gene

# Set annotation for sample groups
sample_groups <- data.frame(
  Condition = c(rep("X", 3), rep("Y", 3))  #CHANGE FOR CONDITION
)
rownames(sample_groups) <- colnames(expr_scaled)

# Load pheatmap
#install.packages("pheatmap")
library(pheatmap)

pheatmap(expr_scaled,
         annotation_col = sample_groups,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         main = "") #title to explain heatmap. 


#### BOXPLOT FOR SOME GENES

# Get the first 5 most differentially expressed gene names (by absolute fold change)
#myTopHits <- myTopHits %>%
  #rownames_to_column("GeneSymbol")

#genes_in_subset <- rownames(expr_subset)
#top5_genes <- myTopHits %>%
  #filter(GeneSymbol %in% genes_in_subset) %>%
  #arrange(desc(abs(logFC))) %>%
  #slice(1:5) %>%
  #pull(GeneSymbol)

#top5_genes
# results: "Tk1"   "Tyms"  "Dnph1" "Ampd3" "Nme1" 


# Let's get the 10 most diff expressed among the most statistically significant
myTopHits <- myTopHits %>% rownames_to_column("GeneSymbol")
genes_in_subset <- rownames(expr_subset)

top12_genes <- myTopHits %>%
  filter(GeneSymbol %in% genes_in_subset) %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%
  slice(1:12) %>%
  pull(GeneSymbol)

top12_genes # gives you the list of the top 12 genes (printed)




# Convert to long format
library(reshape2)

expr_long <- melt(expr_subset)
colnames(expr_long) <- c("Gene", "Sample", "Expression")
expr_long$Condition <- ifelse(grepl("X", expr_long$Sample), "X", "Y") #edit

# Plot for the selected genes
library(ggplot2)

genes_to_plot <- c("")  # Fill with actual genes from top 12

ggplot(expr_long %>% filter(Gene %in% genes_to_plot),
       aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "Expression of Selected Genes in XXX Pathways")


# EXPLORE GENE FUNCTION, MAKE A TABLE OF THE TOP 12 with functions

library(biomaRt)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") #get the right one from ensembl!
genes_of_interest <- c("") # Fill with actual genes from top 12

gene_info <- getBM(
  attributes = c("external_gene_name", "description", "gene_biotype", "entrezgene_id"),
  filters = "external_gene_name",
  values = genes_of_interest,
  mart = mart
)

library(kableExtra)

gene_info %>%
  kable("html", caption = "Gene Functions and Annotations") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
