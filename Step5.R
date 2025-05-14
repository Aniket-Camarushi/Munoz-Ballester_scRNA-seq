#Create a volcano plot for effect of radiation (based on logFC)
vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) + #only the pairwise comparison
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       #subtitle = "Effect of Per2KO when p53YC(5B)", #CHANGE THIS TITLE!!!
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)



results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)

# Subset to only the samples of interest
v.DEGList.filtered.norm$E <- v.DEGList.filtered.norm$E[, samples_of_interest]

# Update column names with sampleLabels (assuming sampleLabels has the same length as the six samples)

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")


