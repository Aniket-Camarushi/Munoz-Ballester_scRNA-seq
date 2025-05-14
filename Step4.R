

# Assign sample labels to DGEList
colnames(myDGEList.filtered.norm) <- sampleLabels # may need to switch this between samples_of_interest and sampleLabels

# Subset the DGEList to only include the selected samples
myDGEList.filtered.norm <- myDGEList.filtered.norm[, samples_of_interest]

# Optionally, ensure the geneIDs are retained (assuming they’re in rownames)
myDGEList.filtered.norm$genes <- myDGEList.filtered.norm$genes[rownames(myDGEList.filtered.norm$counts), , drop = FALSE]

# Verify the changes
myDGEList.filtered.norm



v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)


fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC") #only the pairwise comparison
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")




