library(ggplot2)
setwd("R:/Documents/r/hpheo1/")
proteodata <- read.csv("tmt.csv", sep=";")
ncol(proteodata)
nrow(proteodata)
#create vectors containing the basic data
accessions <- proteodata[,1]
foldchange <- proteodata[,2]
pval <- proteodata[,3]
#log2 transform foldchange
foldchange2 <- log(foldchange,2)
#-log10 transform pval
pval2 <- -log10(pval)
#create the data frame
vulkan <- data.frame(accessions, foldchange2, pval2)
#basis for the plot
basic_vulkan <- ggplot(data=vulkan, aes(x=foldchange2, y=pval2)) +
  geom_point(aes(color = ifelse(foldchange2 > 1 & pval2 > 1.31, "Condition1",
  ifelse(foldchange2 < -1 & pval2 > 1.31, "Condition2", "Other"))))
#set colors
basic_vulkan <- basic_vulkan + scale_color_manual(values = c("Condition1" = "red", 
  "Condition2" = "blue", "Other" = "grey"))
#title, legend, labels appearance setting
basic_vulkan + ggtitle("hPheo1 parental vs. SDHB KO proteomic data") +
  scale_color_manual(values = c("Condition1" = "blue", "Condition2" = "red",
  "Other" = "grey"), name = "Quantified protein groups", # Change the legend title
  labels = c("Condition1" = "Significantly upregulated in SDHB KO", #Change legend labels
  "Condition2" = "Significantly downregulated in SDHB KO", 
  "Other" = "Not significant")) + labs(x = "log(2) fold change", y = "-log(10)p-value") +
  theme(plot.title = element_text(hjust=0.5)) + #plot title justification
  geom_text(data = subset(vulkan, accessions == "P78540"), aes(label = "Arg2"), #label a datapoint
  vjust = -1, hjust = 1, check_overlap = TRUE) + theme_minimal()