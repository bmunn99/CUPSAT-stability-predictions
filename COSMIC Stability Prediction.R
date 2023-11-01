library(ggplot2)
library(tidyverse)

# Load in data from COSMIC
cosmic_data <- read.csv("~/COSMIC AGO3 variants.csv", sep = ',', header = TRUE)
highlight_positions <- c(225, 507, 638)
# Graph the data by position and counts
p <- ggplot(cosmic_data, aes(x=Position, y=Count)) +
  geom_bar(stat='identity', width = 1, fill='black') +  # Plots the data points
  geom_vline(xintercept = highlight_positions, color="red", linetype="dashed", size=0.5) + # Highlight positions in red
  xlim(0, 860) +  # Set x-axis limits
  labs(x="Position", y="Count") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20))# Optional theme for cleaner look
print(p)  # Display the plot
ggsave("~/COSMIC_whole.png", height = 5, width = 10, p, dpi = 300)

### E638A
# Filter to include only the residues of interest
cosmic_subset <- cosmic_data %>% filter(between(Position, 618, 658))
highlight_E638A <- (638)
# Graph the data 
p2 <- ggplot(cosmic_subset, aes(x=Position, y=Count)) +
  geom_bar(stat='identity', width = 0.5, fill='black') +  # Plots the data points
  geom_vline(xintercept = highlight_E638A, color="red", linetype="dashed", size=1) + # Highlight positions in red
  labs(x="Position", y="Count") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20))# Optional theme for cleaner look
p2
ggsave("~/COSMIC_638.png", p2, dpi = 300)

### P1
# Filter to include only the residues of interest
cosmic_subset <- cosmic_data %>% filter(between(Position, 487, 527))
highlight_P1 <- (507)
# Graph the data 
p3 <- ggplot(cosmic_subset, aes(x=Position, y=Count)) +
  geom_bar(stat='identity', width = 0.5, fill='black') +  # Plots the data points
  geom_vline(xintercept = highlight_P1, color="red", linetype="dashed", size=1) + # Highlight positions in red
  labs(x="Position", y="Count") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20))# Optional theme for cleaner look
p3
ggsave("~/COSMIC_507.png", p3, height = 10, width = 10, dpi = 300)

### P2
# Filter to include only the residues of interest
cosmic_subset <- cosmic_data %>% filter(between(Position, 205, 245))
highlight_P2 <- (225)
# Graph the data 
p4 <- ggplot(cosmic_subset, aes(x=Position, y=Count)) +
  geom_bar(stat='identity', width = 0.5, fill='black') +  # Plots the data points
  geom_vline(xintercept = highlight_P2, color="red", linetype="dashed", size=1) + # Highlight positions in red
  labs(x="Position", y="Count") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20))# Optional theme for cleaner look
p4
ggsave("~/COSMIC_225.png", p4, height = 10, width = 10, dpi = 300)
