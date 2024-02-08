library(ggplot2)
setwd("R:/Documents/r/czppgl/")
proteodata <- read.csv("czppgl_corr.csv", sep=";")
#clreate the basic scatterplot with transparency
basic_scatterplot <- ggplot(data=proteodata, aes(x=RNAseq, y=proteo)) +
  geom_point(aes(alpha=I(0.7), color = ifelse(RNAseq > 1 & proteo > 1, "Condition1",
  ifelse(RNAseq < -1 & proteo < -1, "Condition2", "Other"))))
#set the colors for significant & non-significant values
scatterplot <- basic_scatterplot + scale_color_manual(values = c("Condition1" = "red", 
  "Condition2" = "blue", "Other" = "grey"))
scatterplot
final <- scatterplot + ggtitle("CZPPGL RNAseq vs. Proteomic data") + 
  scale_color_manual(values = c("Condition1" = "blue", "Condition2" = "red",
  "Other" = "grey"), name = "Identified Genes", # Change the legend title
  labels = c("Condition1" = "Upregulated in RNAseq & Proteomic Data", #Change legend labels
  "Condition2" = "Downregulated in RNAseq & Proteomic Data", 
  "Other" = "No Consistent Difference")) + labs(x = "RNAseq log(2) fold change", #Change axis labels
  y = "Proteomics log(2) fold change") + theme_minimal() + guides(alpha = FALSE) + #remove alpha label
  scale_x_continuous(breaks = seq(-8, 8, by = 1)) + #set scales and grid steps
  scale_y_continuous(breaks = seq(-15, 15, by = 1))
final

