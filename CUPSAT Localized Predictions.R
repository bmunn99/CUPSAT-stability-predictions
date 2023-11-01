library(ggplot2)
library(tidyverse)
### Amino acid change
# Read in data from CUPSAT predictions: https://cupsat.brenda-enzymes.org/
df <- read.csv("/Users/bmunn99/OneDrive/Clemson/Thesis Figures/Additional Files/CUPSAT predictions.csv", sep = ',', header = TRUE)
# Remove "C" from the chain column
df <- df %>%
  filter(Chain != "C")
# If the amino acid changes are categorical and there are a lot of them, it might be useful to factor it
# This way you ensure that ggplot treats it as a categorical variable
df$Amino.acid <- as.factor(df$Amino.acid)
df$Residue.ID <- as.numeric(df$Residue.ID)

### Normalize the data for amino acid changes at each residue of interest
# Compute the mean and standard deviation of delta_G values
df <- df %>%
  group_by(Residue.ID) %>%
  mutate(
    mean_delta_G = mean(Predicted.DDG..kcal.mol.),
    sd_delta_G = sd(Predicted.DDG..kcal.mol.),
    normalized_delta_G = (Predicted.DDG..kcal.mol. - mean_delta_G) / sd_delta_G
  ) %>%
  ungroup()

# Plot the whole data as a scatterplot
plot <- ggplot(df, aes(x=Residue.ID, y=normalized_delta_G, color=Amino.acid)) +
  geom_boxplot(width = 0.3, height = 0, alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), 
    axis.title = element_text(size = 14), 
    axis.text.y = element_text(size = 12)) +
  labs(x="Residue ID", y="Localized Z-Score", color="Amino Acid Change")  
plot
# Save the plot
ggsave("/Users/bmunn99/Desktop/Alignments/CUPSAT Predictions Localized AGO3.png", plot, width = 10, height = 10, dpi = 300)


### P2
# Filter to include only the residues of interest
df_subset <- df %>% filter(between(Residue.ID, 224, 226))

# Create a new column for color
df_subset$color <- ifelse(df_subset$Amino.acid == "LEU" & df_subset$Residue.ID == 225, "red", "green")
# Create a new column for point size
df_subset$point_size <- ifelse(df_subset$Amino.acid == "LEU" & df_subset$Residue.ID == 225, 1.5, 1)

# Ensure that variables are in the correct format
df_subset$Residue.ID <- as.factor(df_subset$Residue.ID)
df_subset$Amino.acid <- as.factor(df_subset$Amino.acid)

# Plot just the region of interest +- 10 in a scatterplot
plot2 <- ggplot(df_subset, aes(x=Residue.ID, y=normalized_delta_G, color=color, size=point_size)) +
  geom_boxplot(aes(color = color), width = 0.3, height = 0, alpha = 0.8) +
  annotate("rect", xmin = as.numeric(2) - 0.5, xmax = as.numeric(2) + 0.5, ymin = as.numeric(0) - 2, ymax = as.numeric(0) + 2, alpha = 0.1, fill = "blue") +
  scale_color_identity() +
  scale_size_identity() +
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25), 
    axis.title.x = element_text(size = 30), 
    axis.title.y = element_blank(), # Remove y-axis title if you are putting graphs on same axis
    axis.text.y = element_text(size = 25)) +
  labs(x="Residue", y="Localized Z-Score", color="Amino Acid Change") +
  scale_y_continuous(limits = c(-2, 2))
plot2
# Save the plot
ggsave("/Users/bmunn99/Desktop/Alignments/CUPSAT Localized 225.png", plot2, width = 10, height = 10, dpi = 300)


### P1
# Filter to include only the residues of interest
df_subset <- df %>% filter(between(Residue.ID, 506, 508))

# Create a new column for color
df_subset$color <- ifelse(df_subset$Amino.acid == "TRP" & df_subset$Residue.ID == 507, "red", "green")
# Create a new column for point size
df_subset$point_size <- ifelse(df_subset$Amino.acid == "TRP" & df_subset$Residue.ID == 507, 1.5, 1)

# Ensure that variables are in the correct format
df_subset$Residue.ID <- as.factor(df_subset$Residue.ID)
df_subset$Amino.acid <- as.factor(df_subset$Amino.acid)

# Plot just the region of interest +- 10 in a scatterplot
plot3 <- ggplot(df_subset, aes(x=Residue.ID, y=normalized_delta_G, color=color, size=point_size)) +
  geom_boxplot(aes(color = color), width = 0.3, height = 0, alpha = 0.8) +
  annotate("rect", xmin = as.numeric(2) - 0.5, xmax = as.numeric(2) + 0.5, ymin = as.numeric(0) - 2, ymax = as.numeric(0) + 2, alpha = 0.1, fill = "blue") +
  scale_color_identity() +
  scale_size_identity() +
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25), 
    axis.title.x = element_text(size = 30), 
    axis.title.y = element_blank(), # Remove y-axis title if you are putting graphs on same axis
    axis.text.y = element_text(size = 25)) +
  labs(x="Residue", y="Localized Z-Score", color="Amino Acid Change") +
  scale_y_continuous(limits = c(-2, 2))
plot3
# Save the plot
ggsave("/Users/bmunn99/Desktop/Alignments/CUPSAT Localized 507.png", plot3, width = 10, height = 10, dpi = 300)



### E638A
# Filter to include only the residues of interest
df_subset <- df %>% filter(between(Residue.ID, 637, 639))

# Create a new column for color
df_subset$color <- ifelse(df_subset$Amino.acid == "ALA" & df_subset$Residue.ID == 638, "red", "green")
# Create a new column for point size
df_subset$point_size <- ifelse(df_subset$Amino.acid == "ALA" & df_subset$Residue.ID == 638, 1.5, 1)

# Ensure that variables are in the correct format
df_subset$Residue.ID <- as.factor(df_subset$Residue.ID)
df_subset$Amino.acid <- as.factor(df_subset$Amino.acid)

# Plot just the region of interest +- 10 in a scatterplot
plot4 <- ggplot(df_subset, aes(x=Residue.ID, y=normalized_delta_G, color=color, size=point_size)) +
  geom_boxplot(aes(color = color), width = 0.3, height = 0, alpha = 0.8) +
  annotate("rect", xmin = as.numeric(2) - 0.5, xmax = as.numeric(2) + 0.5, ymin = as.numeric(0) - 2, ymax = as.numeric(0) + 2, alpha = 0.1, fill = "blue") +
  scale_color_identity() +
  scale_size_identity() +
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25), 
    axis.title = element_text(size = 30), 
    axis.text.y = element_text(size = 25)) +
  labs(x="Residue", y="Localized Z-Score", color="Amino Acid Change") + 
  scale_y_continuous(limits = c(-2, 2))
plot4
# Save the plot
ggsave("/Users/bmunn99/Desktop/Alignments/CUPSAT Localized 638.png", plot4, width = 10, height = 10, dpi = 300)
