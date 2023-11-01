library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Read in a df of gene names from Excel
gene_df <- read.csv('~/P2 diff expressed genes.csv', header = FALSE)
# Convert to a list
gene_list <- as.list(gene_df$V1)

# Convert your gene names to ENTREZ IDs
gene_list_entrez <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_vector <- gene_list_entrez$ENTREZID

# Create the GO term
ego <- enrichGO(gene           = gene_vector,
                OrgDb           = org.Hs.eg.db,
                ont             = "BP", 
                pAdjustMethod   = "BH", 
                qvalueCutoff    = 0.2, 
                readable        = TRUE)

# Extract the number of genes in each pathway
df$GeneCount <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[[", 1))

# Calculate -log10 of the FDR
df$negLog10FDR <- -log10(df$p.adjust)

# Calculate fold enrichment
df$FoldEnrichment <- df$GeneRatio / df$BgRatio

# Extract the actual ratios from the 'x/y' format
df$GeneRatioNum <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[[", 1)) / 
  as.numeric(sapply(strsplit(df$GeneRatio, "/"), "[[", 2))

df$BgRatioNum <- as.numeric(sapply(strsplit(df$BgRatio, "/"), "[[", 1)) / 
  as.numeric(sapply(strsplit(df$BgRatio, "/"), "[[", 2))

# Calculate Fold Enrichment
df$FoldEnrichment <- df$GeneRatioNum / df$BgRatioNum

# Lollipop plot using ggplot2
lollipop <- ggplot(df, aes(x=reorder(Description, GeneCount), y=GeneCount)) + 
  geom_segment(aes(xend=Description, yend=0), color="black") +
  geom_point(aes(size=FoldEnrichment, color=negLog10FDR), alpha=1) +
  coord_flip() +
  theme_minimal() +
  labs(x="Biological Process", y="Gene Count") +
  theme(
    axis.title = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 18)
  ) +
  scale_color_gradient(name="-log10(FDR)", low='blue', high='red') +
  theme(legend.position="right", legend.box="vertical", 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
lollipop
ggsave('~/Lollipop from R.png', lollipop, dpi = 300)

write.csv(df, '~/GO.csv', row.names = TRUE)
