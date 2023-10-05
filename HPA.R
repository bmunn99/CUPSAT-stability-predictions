# Load library and HPA data
library(ggplot2)
library(dplyr)
library(viridis)
hpa <- read.table("~/normal_tissue.tsv", sep = '\t', header = TRUE) # https://www.proteinatlas.org/about/download

# Filter rows with Gene.name == 'AGO3'
hpa_filtered <- hpa[hpa$Gene.name == 'AGO3', ]
hpa_filtered <- subset(hpa_filtered, Tissue != "soft tissue 2")
hpa_filtered <- subset(hpa_filtered, Tissue != "stomach 1")
# Group by Tissue and filter
hpa_filtered <- hpa_filtered %>%
  group_by(Tissue) %>%
  filter(!(n_distinct(Level) > 1 & Level == "Not detected"))

# Create new columns 'Tissue.group'
hpa_filtered <- hpa_filtered %>%
  mutate(Tissue.group = case_when(
    grepl("cerebral cortex|cerebellum|hippocampus|caudate", Tissue, ignore.case = TRUE) ~ "Brain",
    grepl("thyroid gland|parathyroid gland|adrenal gland", Tissue, ignore.case = TRUE) ~ "Gland",
    grepl("nasopharynx|bronchus|lung", Tissue, ignore.case = TRUE) ~ "Respiratory",
    grepl("oral mucosa|salivary gland|esophagus", Tissue, ignore.case = TRUE) ~ "Mucosal",
    grepl("stomach|duodenum|small intestine|colon|rectum", Tissue, ignore.case = TRUE) ~ "Gastrointestinal",
    grepl("pancreas", Tissue, ignore.case = TRUE) ~ "Pancreas",
    grepl("liver|gallbladder", Tissue, ignore.case = TRUE) ~ "Liver",
    grepl("kidney|urinary bladder", Tissue, ignore.case = TRUE) ~ "Kidney",
    grepl("testis|epididymis|seminal vesicle|prostate", Tissue, ignore.case = TRUE) ~ "Male",
    grepl("vagina|ovary|fallopian tube|endometrium|cervix|placenta|breast", Tissue, ignore.case = TRUE) ~ "Female",
    grepl("heart muscle|smooth muscle|skeletal muscle", Tissue, ignore.case = TRUE) ~ "Muscle",
    grepl("soft tissue|adipose tissue", Tissue, ignore.case = TRUE) ~ "Fat",
    grepl("appendix|spleen|lymph node|tonsil|bone marrow", Tissue, ignore.case = TRUE) ~ "Immune",
    grepl("skin", Tissue, ignore.case = TRUE) ~ "Epithelial",
    TRUE ~ "Other"
  ))

# Reorder levels of 'Level' factor
hpa_filtered$Level <- factor(hpa_filtered$Level, levels = c("Not detected", "Low", "Medium", "High"))
                                 
# Remove duplicates of the same tissue
hpa_filtered <- distinct(hpa_filtered, Tissue, .keep_all = TRUE)

# Order by Tissue.group and Tissue
hpa_filtered <- hpa_filtered %>%
  arrange(Tissue.group, Tissue)

# Convert Tissue to a factor respecting the order in the dataframe
hpa_filtered$Tissue <- factor(hpa_filtered$Tissue, levels = unique(hpa_filtered$Tissue))

# Create the bar plot without legend
p <- ggplot(hpa_filtered, aes(x=Tissue, y=Level, fill=Tissue.group)) + 
  geom_bar(stat='identity') +
  scale_fill_viridis(discrete = TRUE, option = "H") +
  labs(x = "Tissue", y = "Protein Score", title = "Protein Expression for AGO3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate x-axis labels for better visibility
        axis.text.y = element_text(size = 12),
        legend.position = "none",  # Remove legend
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) 
p
ggsave("~/HPA AGO3.png", p, width = 10, height = 10, dpi = 300)
