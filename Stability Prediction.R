library(ggplot2)
library(tidyverse)
# Read in data from CUPSAT predictions: https://cupsat.brenda-enzymes.org/
df <- read.csv("/Users/Desktop/Alignments/CUPSAT predictions.csv", sep = ',', header = TRUE)
# Remove "C" from the chain column
df <- df %>%
  filter(Chain != "C")
# If the amino acid changes are categorical and there are a lot of them, it might be useful to factor it
# This way you ensure that ggplot treats it as a categorical variable
df$Amino.acid <- as.factor(df$Amino.acid)
df$Residue.ID <- as.numeric(df$Residue.ID)

### Normalize the data for each residue of interest
# Compute the mean and standard deviation of delta_G values
mean_delta_G <- mean(df$Predicted.DDG..kcal.mol.)
std_dev_delta_G <- sd(df$Predicted.DDG..kcal.mol.)

# Normalize the delta_G values
df <- df %>%
  mutate(normalized_delta_G = (Predicted.DDG..kcal.mol. - mean_delta_G) / std_dev_delta_G)

# Plot the whole data as a scatterplot
plot <- ggplot(df, aes(x=Residue.ID, y=Predicted.DDG..kcal.mol., color=Amino.acid)) +
  geom_boxplot(width = 0.3, height = 0, alpha = 0.5) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12), 
    axis.title = element_text(size = 14), 
    axis.text.y = element_text(size = 12)) +
  ggtitle("CUPSAT Predictions for AGO3 Substitutions") + 
  labs(x="Residue ID", y="∆∆G Value", color="Amino Acid Change")  
plot
# Save the plot
ggsave("/Users/Desktop/Alignments/CUPSAT Predictions AGO3.png", plot, width = 10, height = 10, dpi = 300)


### P1
# Filter to include only the residues of interest
df_subset <- df %>% filter(between(Residue.ID, 220, 230))

# Create a new column for color
df_subset$color <- ifelse(df_subset$Amino.acid == "LEU" & df_subset$Residue.ID == 225, "red", "green")

# Ensure that variables are in the correct format
df_subset$Residue.ID <- as.factor(df_subset$Residue.ID)
df_subset$Amino.acid <- as.factor(df_subset$Amino.acid)

# Get the y limits for the rectangle
y_min <- min(df_subset$normalized_delta_G[df_subset$Residue.ID == 225] -1)
y_max <- max(df_subset$normalized_delta_G[df_subset$Residue.ID == 225] +1)

# Plot just the region of interest +- 10 in a scatterplot
plot2 <- ggplot(df_subset, aes(x=Residue.ID, y=normalized_delta_G, color=color)) +
  geom_boxplot(aes(color = color), width = 0.3, height = 0, alpha = 0.8) +
  annotate("rect", xmin = as.numeric(6) - 0.5, xmax = as.numeric(6) + 0.5, ymin = y_min, ymax = y_max, alpha = 0.1, fill = "red") +
  scale_color_identity() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20), 
    axis.title = element_text(size = 20), 
    axis.text.y = element_text(size = 20)) +
  labs(x="Residue", y="Normalized ∆∆G", color="Amino Acid Change")
plot2
# Save the plot
ggsave("/Users/Desktop/Alignments/CUPSAT 225.png", plot2, width = 10, height = 10, dpi = 300)


### P2
# Filter to include only the residues of interest
df_subset <- df %>% filter(between(Residue.ID, 502, 512))

# Create a new column for color
df_subset$color <- ifelse(df_subset$Amino.acid == "TRP" & df_subset$Residue.ID == 507, "red", "green")

# Ensure that variables are in the correct format
df_subset$Residue.ID <- as.factor(df_subset$Residue.ID)
df_subset$Amino.acid <- as.factor(df_subset$Amino.acid)

# Get the y limits for the rectangle
y_min <- min(df_subset$normalized_delta_G[df_subset$Residue.ID == 507] -1)
y_max <- max(df_subset$normalized_delta_G[df_subset$Residue.ID == 507] +1)

# Plot just the region of interest +- 10 in a scatterplot
plot3 <- ggplot(df_subset, aes(x=Residue.ID, y=normalized_delta_G, color=color)) +
  geom_boxplot(aes(color = color), width = 0.3, height = 0, alpha = 0.8) +
  annotate("rect", xmin = as.numeric(6) - 0.5, xmax = as.numeric(6) + 0.5, ymin = y_min, ymax = y_max, alpha = 0.1, fill = "red") +
  scale_color_identity() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20), 
    axis.title = element_text(size = 20), 
    axis.text.y = element_text(size = 20)) +
  labs(x="Residue", y="Normalized ∆∆G", color="Amino Acid Change")
plot3
# Save the plot
ggsave("/Users/Desktop/Alignments/CUPSAT 507.png", plot3, width = 10, height = 10, dpi = 300)



### E638A
# Filter to include only the residues of interest
df_subset <- df %>% filter(between(Residue.ID, 633, 643))

# Create a new column for color
df_subset$color <- ifelse(df_subset$Amino.acid == "ALA" & df_subset$Residue.ID == 638, "red", "green")

# Ensure that variables are in the correct format
df_subset$Residue.ID <- as.factor(df_subset$Residue.ID)
df_subset$Amino.acid <- as.factor(df_subset$Amino.acid)

# Get the y limits for the rectangle
y_min <- min(df_subset$normalized_delta_G[df_subset$Residue.ID == 638] -1)
y_max <- max(df_subset$normalized_delta_G[df_subset$Residue.ID == 638] +1)

# Plot just the region of interest +- 10 in a scatterplot
plot4 <- ggplot(df_subset, aes(x=Residue.ID, y=normalized_delta_G, color=color)) +
  geom_boxplot(aes(color = color), width = 0.3, height = 0, alpha = 0.8) +
  annotate("rect", xmin = as.numeric(6) - 0.5, xmax = as.numeric(6) + 0.5, ymin = y_min, ymax = y_max, alpha = 0.1, fill = "red") +
  scale_color_identity() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20), 
    axis.title = element_text(size = 20), 
    axis.text.y = element_text(size = 20)) +
  labs(x="Residue", y="Normalized ∆∆G", color="Amino Acid Change")
plot4
# Save the plot
ggsave("/Users/Desktop/Alignments/CUPSAT 638.png", plot4, width = 10, height = 10, dpi = 300)


