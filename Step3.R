# Define samples of interest
samples_of_interest <- c("XXX1", "XXX2", "XXX3", "YYY1", "YYY2", "YYY3")

# Subset log2.cpm.filtered.norm to only include samples of interest
log2.cpm.filtered.norm.subset <- log2.cpm.filtered.norm[, samples_of_interest]

# Step3_multivariate ----
pca.res <- prcomp(t(log2.cpm.filtered.norm.subset), scale. = FALSE, retx = TRUE)
pc.var <- pca.res$sdev^2  # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var / sum(pc.var) * 100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.res.df$Sample <- samples_of_interest
pca.res.df$Group <- group[sampleLabels %in% samples_of_interest]

# PCA Plot
library(ggrepel)

pca.plot <- ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, label = Sample, color = Group, group = Group) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, box.padding = 0.5, max.overlaps = Inf) + # <-- this avoids overlap
  geom_text(vjust = -1) +  # Adds sample labels near points
  xlab(paste0("PC1 (", pc.per[1], "%", ")")) + 
  ylab(paste0("PC2 (", pc.per[2], "%", ")")) +
  labs(title = "PCA plot", caption = paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)


# Subset data and calculate averages + logFC
mydata.df <- log2.cpm.filtered.norm.df %>%
  mutate(
    XXX.AVG = rowMeans(across(all_of(samples_of_interest[1:3]))),
    YYY.AVG = rowMeans(across(all_of(samples_of_interest[4:6]))),
    LogFC = XXX.AVG - YYY.AVG
  ) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(geneID, all_of(samples_of_interest), XXX.AVG, YYY.AVG, LogFC)

# Render datatable
#datatable(mydata.df,
#          extensions = c('KeyTable', "FixedHeader"),
#          filter = 'top',
#          options = list(keys = TRUE,
#                         searchHighlight = TRUE,
#                         pageLength = 10,
#                         lengthMenu = c("10", "25", "50", "100")))

# Group and design matrix
group <- factor(targets$type[targets$sample %in% samples_of_interest])
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

contrast.matrix <- makeContrasts(states = XXX - YYY, #keeping it simple
                                 levels=design)

#######################################################################################

