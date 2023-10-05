library(dplyr)
library(multiMiR)
library(ggplot2)
library(tidyr)
library(ggpubr)

# Read in P2 DGE and lc data from "RNA-seq with edgeR" script, created w/o logCPM filter
P2 <- read.csv("/Users/bmunn99/Desktop/P2 Data for Bootstrap/P2.csv", sep = ',', header = TRUE) 
# Read in miRTarBase data
miRNA_data <- read.csv("/Users/bmunn99/Desktop/hsa_MTI.csv", sep = ',', header = TRUE)
# Read in ECM_mirnas.csv
ECM_mirnas <- read.csv("/Users/bmunn99/Desktop/P2 Pathways/Cell Projection Development/unique_mirnas_unfiltered.csv", sep = ',', header = TRUE)
# Rename "mature_mirna_id" column in dataframes to "miRNA"
colnames(ECM_mirnas)[colnames(ECM_mirnas) == "mature_mirna_id"] <- "miRNA"
# Rename "target_symbol" column in dataframes to "Target Gene"
colnames(ECM_mirnas)[colnames(ECM_mirnas) == "target_symbol"] <- "Target Gene"

# Use the same logCPM filter for P2 as for ECMs
P2 <- P2 %>%
  filter(logCPM>=8)

### Create a separate function that extracts miRNAs from P2 and then returns the top 10 occurring miRNAs, then computes average for each list
#Pull random miRNAs from miRNA_data
get_mirnas <- function(gene) {
  miRNAs <- miRNA_data$miRNA[miRNA_data$Target.Gene == gene]
  # Make sure miRNAs is always length 10
  length(miRNAs) <- 10
  return(miRNAs)
}

# Bootstrap function for random miRNAs
random_bootstrap_top10_avg <- function(data, n) {
  replicate(n, {
    # Step 1: Sample 15 genes.
    genes <- sample(unique(data$gene_name), 14) # Sample the amount of genes you have
    
    # Step 2: Add miRNA column.
    df <- data.frame(gene_name = genes) %>%
      rowwise() %>%
      mutate(miRNA = list(get_mirnas(gene_name)))
    
    # Step 3: Count occurrences.
    df <- df %>%
      unnest(miRNA) %>%
      count(miRNA) %>%
      arrange(-n)  # Arrange in descending order
    
    # Step 4: Extract top 10 most occurring miRNAs
    top10 <- head(df, 10)
    
    # Step 5: Compute average occurrence of top 10
    avg_occurrence <- mean(top10$n)
    
    return(data.frame(miRNA = "average_top10", total_occurrences = avg_occurrence))
  }, simplify = FALSE) %>%
    bind_rows()  # Bind the dataframes into one
}

# Bootstrap each data frame 1000 times.
random_bootstrap_P2_top10 <- random_bootstrap_top10_avg(P2, 1000)
random_bootstrap_P2_top10 <- random_bootstrap_P2_top10 %>%
  mutate(Source = 'Random_P2')

# Bootstrap function for known miRNAs
filtered_bootstrap_top10_avg <- function(data, n) {
  # List of specific miRNAs to keep
  selected_miRNAs <- c('hsa-miR-129-2-3p','hsa-miR-1-3p','hsa-miR-124-3p','hsa-let-7b-5p','hsa-miR-107',
                       'hsa-miR-155-5p','hsa-miR-101-3p','hsa-miR-182-5p','hsa-miR-103a-3p','hsa-miR-16-5p')
  
  replicate(n, {
    # Step 1: Sample 15 genes.
    genes <- sample(unique(data$gene_name), 14) # Sample the amount of genes you have
    
    # Step 2: Add miRNA column.
    df <- data.frame(gene_name = genes) %>%
      rowwise() %>%
      mutate(miRNA = list(get_mirnas(gene_name)))
    
    # Step 3: Count occurrences.
    df <- df %>%
      unnest(miRNA) %>%
      count(miRNA) %>%
      filter(miRNA %in% selected_miRNAs)  # Step 4: Filter by selected miRNAs
    
    # Calculate average occurrences for this bootstrap sample
    avg_occurrence <- mean(df$n, na.rm = TRUE)
    
    return(data.frame(miRNA = "average", total_occurrences = avg_occurrence))
  }, simplify = FALSE) %>%
    bind_rows()  # Bind the dataframes into one
}

# Bootstrap each data frame 1000 times
filtered_bootstrap_P2_top10 <- filtered_bootstrap_top10_avg(P2, 1000) %>%
  mutate(Source = 'Filtered_P2')

# Prepare the ECM_mirnas data
selected_miRNAs <- c('hsa-miR-129-2-3p','hsa-miR-1-3p','hsa-miR-124-3p','hsa-let-7b-5p','hsa-miR-107',
                     'hsa-miR-155-5p','hsa-miR-101-3p','hsa-miR-182-5p','hsa-miR-103a-3p','hsa-miR-16-5p')
ECM_mirnas <- ECM_mirnas %>%
  filter(miRNA %in% selected_miRNAs)
ECM_mirnas <- ECM_mirnas %>%
  dplyr::select(miRNA, total_occurrences) %>%
  mutate(Source = 'ECM_mirnas')
# Remove duplicate rows so each row represents a single miRNA with the total_occurrences column intact
ECM_mirnas <- ECM_mirnas %>%
  distinct()
ECM_mirnas <- ECM_mirnas %>%
  mutate(miRNA = 'average')

# Bootstrap function to obtain top 10 miRNAs from a 1000 lists and bootstrap 1000 times
new_filtered_bootstrap_top10_avg <- function(data, n) {
  # Step 1: Sample 14 genes for the initial bootstrap.
  genes_initial <- sample(unique(data$gene_name), 14) # Sample the amount of genes you have
  
  # Step 2: Get miRNAs and identify top 10 for the initial bootstrap.
  df_initial <- data.frame(gene_name = genes_initial) %>%
    rowwise() %>%
    mutate(miRNA = list(get_mirnas(gene_name))) %>%
    unnest(miRNA) %>%
    count(miRNA)
  
  # Extract top 10 miRNAs from the initial bootstrap
  top10_miRNAs <- df_initial %>% arrange(desc(n)) %>% slice_head(n = 10) %>% pull(miRNA)
  
  replicate(n, {
    # Step 1: Sample 15 genes for subsequent iterations.
    genes <- sample(unique(data$gene_name), 14) # Sample the amount of genes you have
    
    # Step 2: Add miRNA column.
    df <- data.frame(gene_name = genes) %>%
      rowwise() %>%
      mutate(miRNA = list(get_mirnas(gene_name)))
    
    # Step 3: Count occurrences.
    df <- df %>%
      unnest(miRNA) %>%
      count(miRNA) %>%
      filter(miRNA %in% top10_miRNAs)  # Use the top10 miRNAs identified in the initial bootstrap
    
    # Calculate average occurrences for this bootstrap sample
    avg_occurrence <- mean(df$n, na.rm = TRUE)
    
    return(data.frame(miRNA = "average", total_occurrences = avg_occurrence))
  }, simplify = FALSE) %>%
    bind_rows()  # Bind the dataframes into one
}

# Bootstrap each data frame 1000 times
new_filtered_bootstrap_P2_top10 <- new_filtered_bootstrap_top10_avg(P2, 1000)
new_filtered_bootstrap_P2_top10 <- new_filtered_bootstrap_P2_top10 %>%
  mutate(Source = 'new_Filtered_P2')
  
### To compare top 10 ECM mirnas vs random top 10 mirnas from bootstrap
# P2
# Bind the bootstrap dataframes together
all_data <- rbind(random_bootstrap_P2_top10, ECM_mirnas)

# Calculate n values for each group
n_P2 <- all_data %>% 
  group_by(Source) %>%
  summarise(n = n())

# Plot
p2 <- ggplot(all_data, aes(x = Source, y = total_occurrences)) + 
  geom_jitter(data = subset(all_data, Source == "Random_P2"), width = 0.3, alpha = 0.3, aes(color = Source)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, color = "red", fill = "red") +
  stat_summary(fun.data = function(y) data.frame(y = mean(y), size = 6, label = round(mean(y), 2)), geom = "text", vjust = -1.5) +
  theme_minimal() +
  labs(title = NULL, y = 'Average Occurrences', x = NULL) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20)) +
  geom_signif(comparisons = list(c('Random_P2', 'ECM_mirnas')),
              map_signif_level = TRUE,
              textsize = 5,
              vjust = -0) +
  geom_signif(comparisons = list(c("Random_P2", "ECM_mirnas")), 
              map_signif_level = FALSE, 
              tip_length = 0, 
              textsize = 5, 
              vjust = -1.5) +
  scale_x_discrete(labels = c('ECM_mirnas' = 'Cell Projection miRNAs', 'Random_P2' = 'Random P2')) +
  scale_y_continuous(limits = c(2.5, 15)) +
  annotate("text", x = n_P2$Source, y = 14, label = paste("n =", n_P2$n), color = "black", size = 5)
p2
ggsave("/Users/bmunn99/Desktop/P2 Pathways/Cell Projection Development/Comparison of ECM miRNAs vs. Random P2.png", p2, width = 10, height = 10, dpi = 300)



### To compare top 10 ECM mirnas in P2 vs top 10 mirnas from a sampling of a bootstrap in P2
# P2
# Bind the bootstrap dataframes together
all_data <- rbind(filtered_bootstrap_P2_top10, new_filtered_bootstrap_P2_top10)

# Calculate n values for each group
n_P2 <- all_data %>% 
  group_by(Source) %>%
  summarise(n = n())

# Compute the IQR and the lower and upper thresholds for outliers
Q1 <- quantile(all_data$total_occurrences, 0.25)
Q3 <- quantile(all_data$total_occurrences, 0.75)
IQR <- Q3 - Q1

lower_threshold <- Q1 - 1.5 * IQR
upper_threshold <- Q3 + 1.5 * IQR

# Filter the data
all_data <- all_data %>%
  filter(total_occurrences >= lower_threshold & total_occurrences <= upper_threshold)

# Plot
p2 <- ggplot(all_data, aes(x = Source, y = total_occurrences)) + 
  geom_jitter(data = subset(all_data, Source == "new_Filtered_P2"), width = 0.3, alpha = 0.3, aes(color = Source), ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, color = "red", fill = "red") +
  stat_summary(fun.data = function(y) data.frame(y = mean(y), size = 6, label = round(mean(y), 2)), geom = "text", vjust = -1.5) +
  theme_minimal() +
  labs(title = NULL, y = 'Average Occurrences', x = NULL) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20)) +
  geom_signif(comparisons = list(c('new_Filtered_P2', 'Filtered_P2')),
              map_signif_level = TRUE,
              textsize = 5,
              vjust = 0) +
  geom_signif(comparisons = list(c("new_Filtered_P2", "Filtered_P2")), 
              map_signif_level = FALSE, 
              tip_length = 0, 
              textsize = 5, 
              vjust = -1.5) +
  scale_x_discrete(labels = c('Filtered_P2' = 'Cell Projection miRNAs in P2', 'new_Filtered_P2' = 'Random miRNAs in P2')) +
  scale_y_continuous(limits = c(0, 18)) +
  annotate("text", x = n_P2$Source, y = 16, label = paste("n =", n_P2$n), color = "black", size = 5) 
p2
ggsave("/Users/bmunn99/Desktop/P2 Pathways/Cell Projection Development/Top ECM miRNAs in P1 vs Top Random miRNAs in P2.png", p2, width = 10, height = 10, dpi = 300)

