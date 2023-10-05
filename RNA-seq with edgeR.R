library(edgeR)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(clusterProfiler)
library(tidyr)
library(pheatmap)

setwd("~/Desktop")

#import read counts generated from feature counts of all samples, ago, ddx23, and controls
BioCC_input <- read.table(file = "P2_three_ctrls_transcript.txt", sep = "\t", header = T)
rownames(BioCC_input) <- BioCC_input$Geneid
BioCC_input[,1:6] <- NULL

### Skip to line 66
### Only if you want to extract logCPMs and particular genes
keep_genes <- rowSums(cpm(BioCC_input) > 1) >= 2
counts_matrix_filtered <- BioCC_input[keep_genes, ]
dge <- DGEList(counts=counts_matrix_filtered)
dge <- calcNormFactors(dge, method = 'TMM')
logcpm_values <- cpm(dge, log=TRUE, normalized.lib.sizes = TRUE)
logcpm_values <- as.data.frame(logcpm_values)
logcpm_values$geneid <- rownames(logcpm_values)
logcpm_values <- logcpm_values[, c("geneid", colnames(logcpm_values)[1:ncol(logcpm_values)-1])]
# Filter the dataframe for the specified genes
filtered_df <- logcpm_values[logcpm_values$geneid %in% c("ENSG00000092847.13", "ENSG00000123908.12", "ENSG00000126070.20", "ENSG00000134698.11"),]
filtered_df <- filtered_df %>%
  mutate(geneid = recode(geneid,
                         'ENSG00000092847.13' = 'AGO1',
                         'ENSG00000123908.12' = 'AGO2',
                         'ENSG00000126070.20' = 'AGO3',
                         'ENSG00000134698.11' = 'AGO4'))
colnames(filtered_df) <- c('geneid', 'P1.1', 'P1.2', 'P2.1', 'P2.2', 'P2.3', 'Female.1', 'Female.2', 'Female.3', 'Female.4', 'Female.5', 'Female.6', 'Female.7',
                           'Male.1', 'Male.2', 'Male.3', 'Female.8', 'Female.9')
# Reshape data to long format
df_long <- gather(filtered_df, key = "sample", value = "logCPM", -geneid)
# Classify samples into groups
df_long$group <- case_when(
  df_long$sample %in% c("P1.1", "P1.2") ~ "P1",
  df_long$sample %in% c("P2.1", "P2.2", 'P2.3') ~ "P2",
  df_long$sample %in% c('Female.1', 'Female.2', 'Female.3', 'Female.4', 'Female.5', 'Female.6', 'Female.7', 'Female.8', 'Female.9') ~ "Female Ctrls",
  df_long$sample %in% c('Male.1', 'Male.2', 'Male.3') ~ "Male Ctrls",
  TRUE                                        ~ "Other"
)
# Calculate average logCPM for each group-gene combination
average_data <- df_long %>%
  group_by(geneid, group) %>%
  summarize(avg_logCPM = mean(logCPM, na.rm = TRUE))
# Plot using ggplot2 with average logCPM values over bars
ago_genes <- ggplot(df_long, aes(x = geneid, y = logCPM, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  # Add average logCPM values on top of the bars using the average_data
  geom_text(data = average_data, aes(y = avg_logCPM, label = sprintf("%.2f", avg_logCPM), group = group), 
            position = position_dodge(width = 0.9),  # Adjust width to match bar width
            vjust = 1.5,                            # Adjust vertical position
            size = 3.5,                              # Adjust text size
            color = "black") +                       # Text color
  
  labs(title = "LogCPM for AGO Genes", x = "Gene ID", y = "LogCPM Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ago_genes
ggsave("~/AGO_logCPM.png", ago_genes, dpi = 300, width = 12, height = 12)

### Resume edgeR analysis
# Adjust this filter for how many samples are being processed
BioCC_input_zero_filt <- BioCC_input[rowSums(BioCC_input == 0) <= 2,] #if the count sum is = 0 in less than OR equal to X, it will get filtered out

#median filtering
Medians = numeric()
for (i in 1:nrow(BioCC_input_zero_filt)){
  Medians[i] <- median(as.numeric(BioCC_input_zero_filt[i,]))
}
BioCC_final <- BioCC_input_zero_filt[Medians>=1,]
rpk <- BioCC_final
dim(rpk)
head(rpk)

#set groups
Groups <- as.factor(c(rep("NA", 1), rep("ctrls", 5), rep("ago", 3), rep("ctrls", 3)))

rpk.norm.g <- DGEList(counts = rpk, group = Groups)

rpk.norm.g <- calcNormFactors(rpk.norm.g, method = "TMM")

norm.counts.rpk.g <- cpm(rpk.norm.g, normalized.lib.sizes = TRUE, log = TRUE)

plotMDS(rpk.norm.g, labels = NULL, top = nrow(norm.counts.rpk.g), gene.selection = "common")

design.mat <- model.matrix(~ 0 + Groups)

rpk.norm.g <- estimateDisp(rpk.norm.g,design.mat)

names(rpk.norm.g)

mean(rpk.norm.g$tagwise.dispersion)

plotBCV(rpk.norm.g)

fit <- glmFit(rpk.norm.g, design.mat)

plotQLDisp(fit)

#contrast statements for glmQLFTesting 
P1vsWT <- makeContrasts(Groupsago-Groupsctrls,levels = design.mat)

test <- glmLRT(fit, contrast = P1vsWT)

ago3_P1_WT <- topTags(test, n=Inf)
topTags(ago3_P1_WT)
summary(dt_ago3_P1_WT_lrt<-decideTestsDGE(test,p.value = 0.05))

write.table(ago3_P1_WT, file = "~/P2_txt", sep = "\t")

tmp <- as.data.frame(ago3_P1_WT)
tmp$gene_id <- rownames(tmp)

#generate biexponential rocket plot
#import gencode annotation from local gtf file, subset it to only gene annotations and prep gene_id column
gtf <- readGFF("~/gencode.v43.annotation.gtf") # https://www.gencodegenes.org/human/release_43.html
genes <- gtf %>%
  filter(type == "gene") %>%
  select(gene_name, gene_id)
genes$gene_id <- gsub("\\.\\d+", "", genes$gene_id)
tmp <- as.data.frame(ago3_P1_WT)
tmp$gene_id <- rownames(tmp)
tmp$gene_id <- gsub("\\.\\d+", "", tmp$gene_id)
rownames(tmp) <- NULL
tmp <- merge(tmp, genes, by="gene_id")
### If you wish to write to csv to use for bootstrapping the miRNAs targeting genes
write.csv(tmp, '/Users/bmunn99/Desktop/P2 Data for Bootstrap/P2.csv', row.names = FALSE)
# Filter out lowly expressed genes
tmp <- tmp %>% filter(logCPM>7)
# Write to csv
write.table(tmp, file = "~/P2_genenames.csv", sep = ",")

### Exclude genes shared between P1 and ctrls_exclude_genes.txt
exclude_genes <- read.table(file = "~/ctrls_exclude_genes.txt", sep = "\t", header = TRUE)
tmp <- anti_join(tmp, exclude_genes, by = "gene_name")
### After filtering by logCPM for leafcutter, write out genes to exclude
exclude_genes <- tmp %>% filter(-1>logFC | logFC>1) %>% count(gene_name) %>%
  write.table(file = "~/ctrls_exclude_genes.txt", sep = "\t")



# Volcano Plot
# Load ggplot
library(ggplot2)
library(ggrepel)
library(dplyr)
# Convert dispersion plot to data frame
# Basic scatter plot: x is "logFC" and y is "PValue"
ggplot(data=tmp, aes(x=logFC, y=PValue)) + geom_point()
# Convert directly in the aes()
p <- ggplot(data=tmp, aes(x=logFC, y=-log10(PValue))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=tmp, aes(x=logFC, y=-log10(PValue))) + geom_point() + theme_minimal()
# Add vertical lines for logFC thresholds and one horizontal line for PValue threshold
p2 <- p + geom_vline(xintercept=c(-1.0, 1.0), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column of NAs
tmp$diffexpressed <- "NO"
# If logFC > 1 and PValue < 0.05, set as "UP"
tmp$diffexpressed[tmp$logFC > 1.0 & tmp$PValue < 0.05] <- "UP"
# If logFC < -1 and PValue < 0.05, set as "DOWN"
tmp$diffexpressed[tmp$logFC < -1.0 & tmp$PValue < 0.05] <- "DOWN"
# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=tmp, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) + geom_point() + theme_minimal()
# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-1.0, 1.0), col="red") +
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
ggplot(data=tmp, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text()
# Finally, we can organize the labels nicely using the "ggrepel" package and geom_text_repel () function
# Load library
library(ggrepel)
options(ggrepel.max.overlaps = Inf)

### To highlight particular genes in the volcano plot
# Add a vector of the genes of interest
cell_and_neuron_projection_genes <- c('ITGA3','SEMA3A','PLXNA2','ITGA6','LGALS1','PRRX1','PTPRG','SLIT2','SPOCK1','NTNG1','CFL1','TRAK1','MAGI2','LPAR1')
cell_projection_genes <- c('RALA', 'ARHGAP24', 'TENM2', 'EPS8', 'AUTS2', 'PFN1', 'ITGA2')
# Create a column in tmp telling which genes to color
cell_and_neuron_projection_genes_df <- tmp[tmp$gene_name %in% cell_and_neuron_projection_genes, ]
cell_projection_genes_df <- tmp[tmp$gene_name %in% cell_projection_genes, ]

# Plot adding up all layers we have seen so far
plot_P2_three_ctrls <- ggplot(data=tmp, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(max.overlaps = 10, min.segment.length = 0) +
  scale_color_manual(values=c("purple", "black", "orange")) +
  geom_vline(xintercept=c(-1.0, 1.0), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  theme(axis.text = element_text(size=20), legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = "none", axis.title = element_text(size = 20)) +
  ggtitle(NULL) + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size = 24)) + xlab("log2FC")
plot_P2_three_ctrls <- plot_P2_three_ctrls +
  geom_point(data = cell_and_neuron_projection_genes_df, aes(x=logFC, y=-log10(PValue)), color = 'red', shape = 18, size = 3) +
  geom_point(data = cell_projection_genes_df, aes(x=logFC, y=-log10(PValue)), color = 'green', shape = 18, size = 3)
plot_P2_three_ctrls
ggsave('~/P2 three ctrls logCPM>7 volcano plot.png', plot_P2_three_ctrls, width = 10, height = 10, dpi = 300)
# To convert dispersion plot to csv for GO Ontology
write.csv(tmp, "~/P2_dispersion.csv")


### To merge logCPMs for each gene with leafcutter introns_clusters.txt
tmp = mat
tmp$gene_id <- rownames(tmp)
tmp$gene_id <- gsub("\\.\\d+", "", tmp$gene_id)
rownames(tmp) <- NULL

#generate normalized count medians and standard deviations
gene_counts_long <- tmp %>% 
  tidyr::pivot_longer(cols = -gene_id, names_to = "Sample", values_to = "Count")
# Rename bam files to have the same header
gene_counts_long$Sample <- gsub("P2_S1.bam", "P2", gene_counts_long$Sample)
gene_counts_long$Sample <- gsub("P2_S2.bam", "P2.1", gene_counts_long$Sample)
gene_counts_long$Sample <- gsub("P2_S3.bam", "P2.2", gene_counts_long$Sample)
gene_counts_long <- gene_counts_long[grepl("P2", gene_counts_long$Sample), ]
medians <- gene_counts_long %>%
  group_by(gene_id) %>%
  summarise(medians = median(Count), sd = sd(Count))
test <- merge(medians, genes, by = "gene_id")

# Reassign <name>_intron_clusters to tmp and filter by deltapsi > 0
tmp <- read.table(file = "~/P2_intron_clusters.txt", sep = "\t", header = TRUE)
tmp <- tmp %>%
  filter(deltapsi > 0)
# Rename gene.y to gene_name and merge tmp with genes, then with test to generate csv with median CPMs
tmp <- tmp %>% rename(gene_name=gene.y)
merged_tmp <- merge(genes, tmp, by="gene_name")
merged_lc <- merge(merged_tmp, test, by="gene_name")
write.table(merged_lc, file="~/P2_lc_cpms.txt", sep = "\t")
