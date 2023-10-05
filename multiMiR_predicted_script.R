library(edgeR)
library(multiMiR)
library(tidyverse)
library(writexl)
library(biomaRt)

# Plug target genes into multiMiR and get validated miRNAs
predicted_mirnas <- get_multimir(org     = 'hsa',
                                 target  = c('ITGA3','SEMA3A','PLXNA2','ITGA6','LGALS1','PRRX1','PTPRG','SLIT2','SPOCK1','NTNG1','CFL1','TRAK1','MAGI2','LPAR1'),
                                 table   = 'validated',
                                 summary = TRUE,
                                 predicted.cutoff = 35,
                                 predicted.cutoff.type = 'p',
                                 predicted.site = 'all')
names(predicted_mirnas)
table(predicted_mirnas@data$type)
predicted_mirnas <- as.data.frame(predicted_mirnas@data)

# Filter rows where 'target_symbol' appears in multiple rows of 'mature_mirna_id'
# Remove duplicate rows based on "target_symbol"
predicted_mirnas <- predicted_mirnas %>%
  distinct(mature_mirna_id, target_symbol, .keep_all = TRUE)
# Group by "mature_mirna_id" and summarize the occurrences
summed_mirna_occurrences <- predicted_mirnas %>%
  group_by(mature_mirna_id) %>%
  summarise(total_occurrences = n())
# Merge the data frames summed_mirna_occurrences and predicted_mirnas
merged_df <- merge(predicted_mirnas, summed_mirna_occurrences, by = "mature_mirna_id")
# Keep only the necessary columns
merged_df <- merged_df[, c("mature_mirna_id", "target_symbol", "total_occurrences")]
# Write to csv the unique mirnas that are unfiltered
write.csv(merged_df, "/Users/bmunn99/Desktop/P2 Pathways/Neuron Projection Development/unique_mirnas_unfiltered.csv", row.names = FALSE)
# Filter by mirnas with occurrences >= 5
merged_df <- merged_df %>%
  filter(total_occurrences >= 5)
unique_mirnas <- merged_df %>%
  select(mature_mirna_id, total_occurrences) %>%
  distinct(mature_mirna_id, total_occurrences)
# Write to csv the unique mirnas filtered by >=5 occurrences
write.csv(unique_mirnas, "/Users/bmunn99/Desktop/P2 Pathways/Neuron Projection Development/unique_mirnas.csv", row.names = FALSE)
# Write to csv the mirnas targeting target genes and their occurrences
write.csv(merged_df, file = "/Users/bmunn99/Desktop/P2 Pathways/Neuron Projection Development/ECM_mirnas.csv", row.names = FALSE)
