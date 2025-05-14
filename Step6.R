#Step6_modules----
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, 
          margins=c(8,20)
          )

#modulePick <- 2 
#myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
#hrsub_up <- hclust(as.dist(1-cor(t(myModule_up), method="pearson")), method="complete") 

#heatmap.2(myModule_up, 
#          Rowv=as.dendrogram(hrsub_up), 
#          Colv=as.dendrogram(clustColumns),
#          labRow = NA,
#          col=myheatcolors, scale="row", 
#          density.info="none", trace="none", 
#          RowSideColors=module.color[module.assign%in%modulePick], 
#          margins=c(8,20)
#          )

#modulePick <- 1 
#myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
#hrsub_down <- hclust(as.dist(1-cor(t(myModule_down), method="pearson")), method="complete") 

#heatmap.2(myModule_down, 
#          Rowv=as.dendrogram(hrsub_down), 
#          Colv=as.dendrogram(clustColumns),
#          labRow = NA,
#          col=myheatcolors, scale="row", 
#          density.info="none", trace="none", 
#          RowSideColors=module.color[module.assign%in%modulePick], 
#          margins=c(8,20)
#          )
