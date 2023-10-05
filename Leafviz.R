library(leafviz)
library(ggplot2) 
library(dplyr)
library(rtracklayer)
library(clusterProfiler)

# Use leafviz to quickly visualize splicing, including sashimi plots
leafviz("/Users/bmunn99/Desktop/RData/P1_P2.RData")

# Merge clusters and introns
P2_intron_clusters <- merge(clusters, introns, by = "clusterID")

# Write to .txt to load with RNA-seq, where you'll create a CPM column
write.table(P2_intron_clusters, file = "/Users/bmunn99/Desktop/dPSI>.1 P1&P2 vs all ctrls/P1_P2_intron_clusters.txt", sep = "\t")
# To quickly graph your clusters
ggplot(P2_intron_clusters, aes(x = deltapsi, y = -log10(FDR))) +
  geom_point() +
  theme_classic() 
tmp <- P2_intron_clusters

### Only run after merging median cpms for lc file
P2_intron_clusters <- read.table(file = "/Users/bmunn99/Desktop/dPSI>.1 P2 vs all ctrls/P2_lc_cpms.txt", sep = "\t", header = TRUE)
# Filter lowly expressed genes
tmp <- P2_intron_clusters %>% filter(medians>7)
# If you have a list of genes you'd like to exclude
exclude_genes <- read.table(file = "/Users/bmunn99/Desktop/dPSI>.1 lc ddx vs ago ctrls/ctrls_exclude_genes.txt", sep = "\t", header = TRUE)
tmp <- anti_join(tmp, exclude_genes, by = "gene_name")

# If you want to filter only for Up regulated genes
tmp <- tmp %>%
  filter(deltapsi > 0)

# Volcano Plot
# Load ggplot
library(ggrepel)
# Convert dispersion plot to data frame
# Basic scatter plot: x is "logFC" and y is "PValue"
ggplot(data=tmp, aes(x=deltapsi, y=-log10(FDR))) + geom_point()
# Convert directly in the aes()
p <- ggplot(data=tmp, aes(x=deltapsi, y=-log10(FDR))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=tmp, aes(x=deltapsi, y=-log10(FDR))) + geom_point() + theme_minimal()
# Add vertical lines for logFC thresholds and one horizontal line for PValue threshold
p2 <- p + geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column of NAs
tmp$diffexpressed <- "NO"
# If logFC > 0.6 and PValue < 0.05, set as "UP"
tmp$diffexpressed[tmp$deltapsi > 0.1 & tmp$FDR < 0.05] <- "UP"
# If logFC < -0.6 and PValue < 0.05, set as "DOWN"
tmp$diffexpressed[tmp$deltapsi < -0.1 & tmp$FDR < 0.05] <- "DOWN"
# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=tmp, aes(x=deltapsi, y=-log10(FDR), col=diffexpressed)) + geom_point() + theme_minimal()
# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
## Change point color
# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("purple", "orange", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_color_manual(values = mycolors)
# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
tmp$delabel <- NA
tmp$delabel[tmp$diffexpressed != "NO"] <- tmp$gene_name[tmp$diffexpressed != "NO"]
ggplot(data=tmp, aes(x=deltapsi, y=-log10(FDR), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text()
# Finally, we can organize the labels nicely using the "ggrepel" package and geom_text_repel () function
# Load library
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
# Plot adding up all layers we have seen so far
plot2 <- ggplot(data=tmp, aes(x=deltapsi, y=-log10(FDR), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(max.overlaps = 10, min.segment.length = 0,) +
  scale_color_manual(values=c("black", "orange", "purple")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme(axis.text = element_text(size=20), legend.title = element_text(size = 20), legend.text = element_text(size = 20), axis.title = element_text(size = 20), legend.position = "none") +
  ggtitle(NULL) + theme(plot.title = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(min(tmp$deltapsi), max(1)))  # Adjust x-axis range
plot2
ggsave('/Users/bmunn99/Desktop/dPSI>.1 P2 vs all ctrls/P2 ago ddx logCPM>7 volcano plot.png', plot2, width = 10, height = 10, dpi = 300)
# To convert dispersion plot to csv for GO Ontology
write.csv(tmp, "/Users/bmunn99/Desktop/dPSI>.1 P1&P2 vs all Ctrls/P1_P2_lc_dispersion.csv")


