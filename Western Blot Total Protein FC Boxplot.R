library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# Load in csv containing densities of western blots
# Must contain three columns: 1) Source, 2) Sample, 3) Density
df <- read.csv('/Users/bmunn99/OneDrive/Clemson/Thesis Figures/Western Blots/Densities.csv', header = TRUE)
df$Sample <- factor(df$Sample, levels=c('WT', 'E638A', 'P1', 'P2'))
# Plot the data in a boxplot
p <- ggplot(df, aes(x=Sample, y=Density, fill=Sample)) + 
  theme_minimal() +
  theme(legend.position = 'none') +
  theme(
    axis.text = element_text(size=16), 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size=20)) +
  geom_boxplot(position=position_dodge(width=0.8), width=0.7, outlier.shape = NA) +  # Boxplot
  geom_text(data=data.frame(Sample=unique(df$Sample), Density=max(df$Density)),
            aes(label="n=6"), vjust=-2, position=position_dodge(width=0.8)) +  # Add "n=6" above each boxplot
  stat_summary(fun=mean, geom="point", shape=4, size=7, color="red", position=position_dodge(width=0.8)) +  # Mean as red cross
  stat_summary(fun=mean, geom="text", 
               vjust=3.5,    # Adjust this value to position the text above the mean point
               aes(label=sprintf("%.2f", ..y..)),  # Format to 2 decimal places
               position=position_dodge(width=0.8)) + # Display mean value
  scale_fill_brewer(palette = 'Set3')
p
ggsave('/Users/bmunn99/OneDrive/Clemson/Thesis Figures/Western Blots/Densities.png', p, height = 10, width = 10, dpi = 300)

