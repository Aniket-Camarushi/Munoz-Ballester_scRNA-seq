#step7_functionalEnrichment----

# gProfiler2 for GO enrichment

# this code chunk assumes you already have a contrast matrix and ebFit object from the Step5_diffGenes script
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
#myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC") #only the pairwise comparison
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
#gost.res <- gost(rownames(myTopHits), organism = "mmusculus", correction_method = "fdr")
# produce an interactive manhattan plot of enriched GO terms
#gostplot(gost.res, interactive = TRUE, capped = FALSE)

# clusterProfiler for GSE

#### 1. Install msigdbr if not already
#install.packages("msigdbr")
library(msigdbr)
library(dplyr)

msigdbr_collections()

# 2. Load and get full Reactome (check the one you need to use, can also manually download .gmt file from GSEA Broad Institute)

reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, gene_symbol)

# 3. Rename columns to match TERM2GENE format
colnames(reactome) <- c("term", "gene")




# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC) #only the pairwise comparison
# construct a named vector
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

mydata.gsea <- mydata.gsea[!is.na(names(mydata.gsea)) & names(mydata.gsea) != ""] #fixes empty gene names

# run GSEA using the 'GSEA' function from clusterProfiler
set.seed(123) #set a random seed so that we can reproducible ordering for our GSEA results below
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=reactome, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)



#install.packages("DT")
library(DT)


# view results as an interactive table
datatable(myGSEA.df,
          extensions = c('KeyTable', "FixedHeader"),
         #caption = 'Signatures enriched in irradiation',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=2)




#### create enrichment plots using the enrichplot package

gseaplot2(myGSEA.res,
          geneSetID = c(), #can choose multiple signatures to overlay in this plot, manually get the numbers of the pathway of interest from myGSEA.df)
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          #title = myGSEA.res$Description[]
) #can also turn off this title



# CAMERA for GSEA

# this chunk assumes you have a contrast matrix and design set-up (see step 5 script), as well as a signature database (I'll use my local one again)
# read in signature database
#reactome <- getGmt("D:/Oncology/RNAseq Per2 p53/m2.cp.reactome.v2024.1.Mm.symbols.gmt", geneIdType=SymbolIdentifier())
#reactome <- geneIds(reactome) #extract as a list
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)

#camera.res <- camera(v.DEGList.filtered.norm$E, reactome, design, contrast.matrix[,1])
#camera.df <- as_tibble(camera.res, rownames = "setName")
#camera.df

# filter based on FDR and display as interactive table
#camera.df <- dplyr::filter(camera.df, FDR<=0.01)

#View(camera.df)

#datatable(camera.df,
#          extensions = c('KeyTable', "FixedHeader"),
#          caption = 'Signatures enriched in radiatio',
#          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
#  formatRound(columns=c(2,4,5), digits=2)

#as before, add a variable that maps up/down regulated pathways with phenotype
#camera.df <- camera.df %>%
  #mutate(phenotype = case_when(
    #Direction == "Up" ~ "NIR_dKO",  #NEED TO CHANGE THESE!!! "Disease"
    #Direction == "Down" ~ "NIR_wildtype")) #NEED TO CHANGE THESE!!! "Healthy"

# graph camera results as bubble chart
#ggplot(camera.df[1:25,], aes(x=phenotype, y=setName)) +
  #geom_point(aes(size=NGenes, color = Direction, alpha=-log10(FDR))) +
  #theme_bw()


#GSVA

# this chunk assumes you have a contrast matrix and design set-up (see step 5 script)
# read in your signature database
#reactome <- read.gmt("D:/Git/RNAseq - MDL FBRI/m2.cp.reactome.v2024.1.Mm.symbols.gmt")

#msigdbr_collections()

# 1. Prepare Reactome gene sets
library(msigdbr)
library(dplyr)
library(GSVA)

reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, gene_symbol)

colnames(reactome) <- c("term", "gene")

# Clean up
reactome$term <- as.character(reactome$term)
reactome$gene <- as.character(reactome$gene)
reactome <- reactome[!is.na(reactome$gene) & reactome$gene != "", ]

# Convert to list format for GSVA
reactome_list <- split(reactome$gene, reactome$term)
# Run GSVA directly, all in one step
gsva.res.reactome <- gsva(expr = v.DEGList.filtered.norm$E,
                          gset.idx.list = reactome_list,
                          method = "gsva",  # "ssgsea", "zscore", etc. are other options
                          min.sz = 5,
                          max.sz = 500,
                          verbose = TRUE)


# Apply linear model to GSVA result
# now using Limma to find significantly enriched gene sets in the same way you did to find diffGenes
# this means you'll be using topTable, decideTests, etc
# note that you need to reference your design and contrast matrix here
fit.reactome <- lmFit(gsva.res.reactome, design)
ebFit.reactome <- eBayes(fit.reactome)

# use topTable and decideTests functions to identify the differentially enriched gene sets
topPaths.reactome <- topTable(ebFit.reactome, adjust ="BH", coef=1, number=50, sort.by="logFC") #only the pairwise comparison
res.reactome <- decideTests(ebFit.reactome, method="global", adjust.method="BH", p.value=0.05, lfc=0.5)
# the summary of the decideTests result shows how many sets were enriched in induced and repressed genes in all sample types
summary(res.reactome)

# pull out the gene sets that are differentially enriched between groups
diffSets.reactome <- gsva.res.reactome[res.reactome[,1] !=0,]

# make a heatmap of differentially enriched gene sets
hr.reactome <- hclust(as.dist(1-cor(t(diffSets.reactome), method="pearson")), method="complete") #cluster rows by pearson correlation
hc.reactome <- hclust(as.dist(1-cor(diffSets.reactome, method="spearman")), method="complete") #cluster columns by spearman correlation

# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
mycl.reactome <- cutree(hr.reactome, k=2)
mycolhc.reactome <- rainbow(length(unique(mycl.reactome)), start=0.1, end=0.9)
mycolhc.reactome <- mycolhc.reactome[as.vector(mycl.reactome)]

# assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]). Type demo.col(20) to see more color schemes.
myheatcol <- colorRampPalette(colors=c("yellow","white","blue"))(100)


# plot the hclust results as a heatmap
heatmap.2(diffSets.reactome,
          Rowv=as.dendrogram(hr.reactome),
          Colv=as.dendrogram(clustColumns),
          col=myheatcol, scale="row",
          density.info="none", trace="none",
          cexRow=0.9, cexCol=1, 
          margins=c(10,36)
          )

# Set the output file path
#output_file <- "diffSets_reactome.csv"

# Write the data frame to CSV
#write.csv(diffSets.reactome, file = output_file, row.names = TRUE)

#cat("diffSets.reactome has been written to", output_file)
