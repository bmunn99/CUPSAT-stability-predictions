# Load library and HPA data
library(ggplot2)
library(dplyr)
library(viridis)

### AGO2
hpa_AGO2 <- read.table("/Users/bmunn99/Downloads/normal_tissue AGO2.csv", sep = ',', header = TRUE) # https://www.proteinatlas.org/about/download
### AGO3
hpa_AGO3 <- read.table("/Users/bmunn99/Downloads/normal_tissue 2.tsv", sep = '\t', header = TRUE) # https://www.proteinatlas.org/about/download

# Filter rows with Gene.name == 'AGO3'
hpa_AGO3 <- hpa_AGO3 %>%
  filter(Gene.name %in% 'AGO3')

# Bind the two dfs together
hpa_both <- rbind(hpa_AGO2, hpa_AGO3)

# Create new columns 'Tissue.group'
hpa_both <- hpa_both %>%
  mutate(Tissue.group = case_when(
    grepl("cerebral cortex|cerebellum|hippocampus|caudate", Tissue, ignore.case = TRUE) ~ "Brain",
    grepl("thyroid gland|parathyroid gland|adrenal gland", Tissue, ignore.case = TRUE) ~ "Gland",
    grepl("nasopharynx|bronchus|lung", Tissue, ignore.case = TRUE) ~ "Respiratory",
    grepl("oral mucosa|salivary gland|esophagus", Tissue, ignore.case = TRUE) ~ "Mucosal",
    grepl("stomach|duodenum|small intestine|colon|rectum", Tissue, ignore.case = TRUE) ~ "Gastrointestinal",
    grepl("pancreas", Tissue, ignore.case = TRUE) ~ "Pancreas",
    grepl("liver|gallbladder", Tissue, ignore.case = TRUE) ~ "Liver",
    grepl("kidney|urinary bladder", Tissue, ignore.case = TRUE) ~ "Kidney",
    grepl("testis|epididymis|seminal vesicle|prostate", Tissue, ignore.case = TRUE) ~ "Male Reproductive",
    grepl("vagina|ovary|fallopian tube|endometrium|cervix|placenta|breast", Tissue, ignore.case = TRUE) ~ "Female Reproductive",
    grepl("heart muscle|smooth muscle|skeletal muscle", Tissue, ignore.case = TRUE) ~ "Muscle",
    grepl("soft tissue|adipose tissue", Tissue, ignore.case = TRUE) ~ "Fat",
    grepl("appendix|spleen|lymph node|tonsil|bone marrow", Tissue, ignore.case = TRUE) ~ "Immune",
    grepl("skin", Tissue, ignore.case = TRUE) ~ "Epithelial",
    TRUE ~ "Other"
  ))

# Filter and summarize the data
hpa_both <- hpa_both %>%
  group_by(Gene.name, Level) %>%
  summarize(count = n())

# Reorder levels of 'Level' factor
hpa_both$Level <- factor(hpa_both$Level, levels = c('High', 'Medium', 'Low', 'Not detected'))

# Plot
hpa_both_plot <- ggplot(hpa_both, aes(x = Gene.name, y = count, fill = Level)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_viridis(discrete = TRUE, option = "H", name = "Protein Score") +
  labs(x = "Protein", y = "Number of Tissues") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),  
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 24))
hpa_both_plot
ggsave('/Users/bmunn99/OneDrive/Clemson/Thesis Figures/HPA/Tissues Per Score.png', hpa_both_plot, dpi = 300)
