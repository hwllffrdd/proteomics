library(dplyr)
library(ggplot2)

setwd("R:/...") #set your working directory

#load the data - this script presumes the use of ProteomeDiscoverer output
data_1 <- read.csv("data_1.csv", sep=";")
data_2 <- read.csv("data_2.csv", sep=";")

#getting rid of proteins identified with just one peptide (column 9)
#if needed, the expression data (col 29) can be reversed
data_1_filtered <- data_1 %>% 
  filter(.[[9]] >= 2) %>%
  mutate(.[[29]] = 1 / .[[29]])

#the same for data_2, here without the 1/x reversion
data_2_filtered <- data_2 %>%
  filter(.[[9]] >= 2)

#cleaning data, so that only statistically significant proteins (cols 44 and 18)
#are present and getting rid of proteins identified with just one peptide (col 9)
data_1_clean <- data_1 %>%
  filter(.[[44]] <= 0.05, .[[9]] >= 2)
data_2_clean <- data_2 %>%
  filter(.[[18]] <= 0.05, .[[9]] >= 2)

#full join of the initial data
combined_data <- full_join(data_1, data_2, by = "Accession")
write.csv(combined_data, "merge_full.csv")
#overlapping data
combined_data1 <- inner_join(data_1_filtered, data_2_filtered, by = "Accession")
#calculation of the correlation coefficient... replace Column_name
correlation_coefficient <- cor(combined_data1$Column_name, combined_data2$Colum_name, use = "pairwise.complete.obs")
correlation_coefficient
#dot plot visualizing the correlation of the two proteomic datasets
ggplot(combined_data1, aes(x = Column_name, y = Column_name)) + 
  geom_point(alpha = 0.5, color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a linear regression line without a confidence band
  scale_x_continuous(name = "Data1 fold change", trans = "log2", limits = c(0.25, 4), breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16)) +  # X-axis on log2 scale
  scale_y_continuous(name = "Data2 fold change", trans = "log2", limits = c(0.25, 4), breaks = c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16))  +
  labs(title = "Comparison of Protein Expression between ...") +
  theme_minimal(base_size = 14) +  # Use a minimal theme with a base font size
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none"  # Hide legend (adjust if you have a legend)
  ) +
  annotate("text", x = Inf, y = Inf, label = paste("R =", round(correlation_coefficient, 3)), 
           hjust = 2, vjust = 25, size = 5, color = "darkred")

#listing the overlapping statistically significant up- and downregulated proteins
#1.5x fold change applied (and reversed: 0.667)
combined_up <- combined_data %>%
  filter(.[[29]] >= 1.5, .[[162]] >= 1.5)
combined_down <- combined_data %>%
  filter(.[[29]] <= 0.667, .[[162]] <= 1.5)
write.csv(combined_up, "merge_up.csv")
write.csv(combined_down, "merge_down.csv")

#EnrichPathway pipeline:
#create uniprot lists
combined_up_list <- as.vector(combined_up$Accession)
combined_down_list <- as.vector(combined_down$Accession)
background <- as.vector(combined_data2$Accession)

#convert them to ensembl gene id
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
entrez_ids_up <- getBM(attributes = c('uniprotswissprot', 'entrezgene_id'),
                       filters = 'uniprotswissprot',
                       values = combined_up_list,
                       mart = ensembl)
entrez_ids_down <- getBM(attributes = c('uniprotswissprot', 'entrezgene_id'),
                         filters = 'uniprotswissprot',
                         values = combined_down_list,
                         mart = ensembl)
entrez_background <- getBM(attributes = c('uniprotswissprot', 'entrezgene_id'),
                         filters = 'uniprotswissprot',
                         values = background,
                         mart = ensembl)

#convert the numbers to characters
list_up <- as.character(entrez_ids_up[["entrezgene_id"]])
list_down <- as.character(entrez_ids_down[["entrezgene_id"]])
list_background <- as.character(entrez_background[["entrezgene_id"]])

library("ReactomePA")

#run the EnrichPathway and create the graphic output
enrich_up <- enrichPathway(list_up, pvalueCutoff=0.05, universe = list_background)
enrich_down <- enrichPathway(list_down, pvalueCutoff=0.05, universe = list_background)
dotplot_up <- dotplot(enrich_up, showCategory=20, font.size=8)
dotplot_down <- dotplot(enrich_down, showCategory=20, font.size=8)
dotplot_up
dotplot_down
#create the summaries
summary_up <- as.data.frame(summary(enrich_up))
summary_down <- as.data.frame(summary(enrich_down))
head(summary_up)
head(summary_down)
write.csv(summary_up, file="summary_common_up.csv")
write.csv(summary_down, file="summary_common_down.csv")
