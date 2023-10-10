library(zoo)
library(dplyr)

# You must first create dot plots of your SHAPE data using rnastructure's "Fold" and "ct2dot" commands https://rna.urmc.rochester.edu/Text/index.html
# Load the dot plots
AGO3_WT_dp <- read.table("/Users/bmunn99/Desktop/SHAPE/Nextseq/rnastructure/AGO3_WT.dp", header=TRUE, skip=1) # adjust skip as necessary
AGO3_E638A_dp <- read.table("/Users/bmunn99/Desktop/SHAPE/Nextseq/rnastructure/AGO3_E638A.dp", header=TRUE, skip=1) # adjust skip as necessary
AGO3_P2_dp <- read.table("/Users/bmunn99/Desktop/SHAPE/Nextseq/rnastructure/AGO3_P2.dp", header=TRUE, skip=1) # adjust skip as necessary
AGO3_P1_dp <- read.table("/Users/bmunn99/Desktop/SHAPE/Nextseq/rnastructure/AGO3_P1.dp", header=TRUE, skip=1) # adjust skip as necessary

# Rename the log10 column to Probability
AGO3_WT_dp <- AGO3_WT_dp %>%
  rename(Probability = X.log10.Probability.)
AGO3_E638A_dp <- AGO3_E638A_dp %>%
  rename(Probability = X.log10.Probability.)
AGO3_P2_dp <- AGO3_P2_dp %>%
  rename(Probability = X.log10.Probability.)
AGO3_P1_dp <- AGO3_P1_dp %>%
  rename(Probability = X.log10.Probability.)
# Transform from -log10
AGO3_WT_dp$Probability <- 10^(-AGO3_WT_dp$Probability)
AGO3_E638A_dp$Probability <- 10^(-AGO3_E638A_dp$Probability)
AGO3_P2_dp$Probability <- 10^(-AGO3_P2_dp$Probability)
AGO3_P1_dp$Probability <- 10^(-AGO3_P1_dp$Probability)

### Pearson correlation
# Step 1 (only if you wish to focus on a particular region of RNAs)
# Filter the data frames
region_WT <- AGO3_WT_dp[AGO3_WT_dp$i >= 750 & AGO3_WT_dp$i <= 790, ]
region_P2 <- AGO3_P2_dp[AGO3_P2_dp$i >= 750 & AGO3_P2_dp$i <= 790, ]

# Step 2: Merge the two filtered data frames and/or RNAs to compare based on 'i' and 'j'
merged_region <- merge(region_WT, region_P2, by = c("i", "j"), suffixes = c("_WT", "_P2"))
AGO3_WT_E638A <- merge(AGO3_WT_dp, AGO3_E638A_dp, by = c('i', 'j'), suffixes = c('_WT', "_E638A"))
# Write to csv
write.csv(AGO3_WT_E638A, '/Users/bmunn99/Desktop/SHAPE/Nextseq/rnastructure/AGO3_WT_E638A.csv', row.names = FALSE)
AGO3_WT_P2 <- merge(AGO3_WT_dp, AGO3_P2_dp, by = c('i', 'j'), suffixes = c('_WT', "_P2"))
# Write to csv
write.csv(AGO3_WT_P2, '/Users/bmunn99/Desktop/SHAPE/Nextseq/rnastructure/AGO3_WT_P2.csv', row.names = FALSE)
AGO3_WT_P1 <- merge(AGO3_WT_dp, AGO3_P1_dp, by = c('i', 'j'), suffixes = c('_WT', "_P1"))
# Write to csv
write.csv(AGO3_WT_P1, '/Users/bmunn99/Desktop/SHAPE/Nextseq/rnastructure/AGO3_WT_P1.csv', row.names = FALSE)

# Step 3: Compute the correlation between the 'Probability' columns
merged_region_correlation_value <- cor(merged_region$Probability_WT, merged_region$Probability_P2, method = "pearson", use = "complete.obs")
WT_E638A_correlation_value <- cor(AGO3_WT_E638A$Probability_WT, AGO3_WT_E638A$Probability_E638A, method = "pearson", use = "complete.obs")
WT_P2_correlation_value <- cor(AGO3_WT_P2$Probability_WT, AGO3_WT_P2$Probability_P2, method = "pearson", use = "complete.obs")
WT_P1_correlation_value <- cor(AGO3_WT_P1$Probability_WT, AGO3_WT_P1$Probability_P1, method = "pearson", use = "complete.obs")

# Print the results
print(merged_region_correlation_value)
print(WT_E638A_correlation_value)
print(WT_P2_correlation_value)
print(WT_P1_correlation_value)

# To index individual nucleotides for further analysis
example <- AGO3_WT_P2[AGO3_WT_P2$i >= 755 & AGO3_WT_P2$i <= 787, ] %>%
  filter(Probability_P2 > 0.4)
