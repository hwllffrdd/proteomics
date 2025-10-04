library(mixOmics)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)

# Load & check your data
setwd("R:/Documents/r/ppgl_tmt")
# wanna read csv or tab-delimited txt?
data <- read.csv("ppgl_tmt.csv", sep=";")
data <- read.table("ppgl_tmt_wo_ctrl_perseus.txt", sep="\t", header=TRUE)

# Assuming 'data' is your dataframe and it has an 'Accession' column
head(data)
exprs <- as.matrix(data[,-1]) # Convert to matrix, excluding the first column
rownames(exprs) <- data$Accession # Assign row names
head(exprs)
write.csv(ppgldata, "ppgldata.csv", row.names = TRUE, quote = TRUE, sep = ";")

# log2 transform the matrix
exprs_log2 <- log2(exprs + 1)

# Create the design matrix for your experimental groups
group <- factor(c(rep("ctrl", 3), rep("nt", 3), rep("ret", 4), rep("sdhb", 4), rep("vhl", 4)))
design <- model.matrix(~0+group) # No intercept, treating groups as different levels
colnames(design) <- levels(group)

# Apply the lmFit and eBayes functions from limma to fit the linear model and compute statistics.
fit <- lmFit(exprs_log2, design)
fit <- eBayes(fit)

# Perform ANOVA-like comparison using limma to test for differences across the groups
# Change the number you are dividing by depending on the number of groups in parentheses
cont.matrix <- makeContrasts(
  CtrlVsOthers = ctrl - (nt + ret + sdhb + vhl)/4,
  NtVsOthers = nt - (ctrl + ret + sdhb + vhl)/4,
  RetVsOthers = ret - (ctrl + nt + sdhb + vhl)/4,
  SdhbVsOthers = sdhb - (ctrl + nt + ret + vhl)/4,
  VhlVsOthers = vhl - (ctrl + nt + ret + sdhb)/4,
  levels=design
)

# Apply contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
# Apply treat for additional stringency on fold change
fit2 <- treat(fit2, lfc=log2(1.5))
# Compute statistics with eBayes
fit2 <- eBayes(fit2)

# Filter out non-significant accessions based on the ANOVA results
results <- topTable(fit2, adjust="BH", number=nrow(exprs), genelist=rownames(exprs)) # Adjust method for P-values, and select all rows
sigResults <- results[results$adj.P.Val < 0.05,] # Adjust threshold as needed

write.csv(results, file="allresults.csv")

sigResults <- sigResults[complete.cases(sigResults), ]
rownames(sigResults) <- sigResults$ID


# Normalize the significant results using Z-score normalization
sigExprs <- exprs[rownames(sigResults),]
zScores <- t(apply(sigExprs, 1, scale))

num_srn_mx <- zScores
head(num_srn_mx)
write.csv(num_srn_mx, file="ppglmatrix.csv")

#create a factor for sPLS-DA's "Y" values:
ppglcat <- group

#transverse the data, so that the individual samples' values are in rows:
ppgldata <- t(num_srn_mx)
head(ppgldata)

#in case if something goes wrong with transversing the data:
#ppgldata0 <- read.csv("transversed.csv", sep=";", row.names=1)
#ppgldata <- data.matrix(ppgldata0)
#str(ppgldata)


#perform sPLS-DA
splsda <- splsda(
  ppgldata,
  group,
  ncomp = 2, keepX = c(50, 25)
)
summary(splsda)
#create the plot for initial sPLS-DA
plotIndiv(splsda, legend = TRUE, ellipse = TRUE, ind.names = FALSE, title = "PPGL TMT sPLS-DA")
#optimize sPLS-DA:
plotVar(splsda, cutoff = 0.8)
list.keepX <- c(5:10,  seq(20, 100, 10))
tune.splsda._ <- tune.splsda(ppgldata, ppglcat, ncomp = 2, 
                             validation = 'Mfold',
                             folds = 3, dist = 'max.dist', progressBar = FALSE,
                             measure = "overall", test.keepX = list.keepX, # option to be concidered: measure = "BER"
                             nrepeat = 50)
error <- tune.splsda._$error.rate
error
ncomp <- tune.splsda._$choice.ncomp$ncomp
ncomp
select.keepX <- tune.splsda._$choice.keepX[1:ncomp]
select.keepX

#create the final sPLS-DA plot with optimal parameters
final_splsda <- splsda(ppgldata, ppglcat, ncomp = ncomp, keepX = select.keepX)
background <- background.predict(final_splsda, comp.predicted=2, dist = "centroids.dist")


# Define the layout: 2 columns, with the second one being narrower for the legend
layout_matrix <- matrix(c(1, 2), ncol = 2, nrow = 1, byrow = TRUE)
layout(layout_matrix, widths = c(4, 1)) # Adjust 'widths' as needed to control space allocation

# Reset the graphics device
par(mar = c(5, 4, 4, 2) + 0.1)  # Reset to default margins

# Set up the plot area: fig=c(left, right, bottom, top)
# First plot uses 80% of width
par(fig = c(0, 0.8, 0, 1))

# Reset any existing graphics parameters
dev.off()  # Only if there's an active plot
par(mar = c(5, 4, 4, 2) + 0.1)  # Default margins

# Create a layout matrix: 1=plot area, 2=legend area
# This gives 80% of width to plot, 20% to legend
layout_matrix <- matrix(c(1, 2), ncol=2, byrow=TRUE)
layout(layout_matrix, widths=c(4, 1))  # 4:1 ratio for plot:legend

# First plot area - the main plot
par(mar = c(5, 4, 4, 0))  # Reduce right margin since legend will be there
plotIndiv(final_splsda, comp=c(1,2), 
          ellipse = FALSE, title="PPGL TMT sPLS-DA", 
          ind.names = FALSE, pch=16,
          legend = FALSE,
          col=c("navy", "darkmagenta", "darkorange", "forestgreen"),
          style="graphics")

# Second plot area - just the legend
par(mar = c(5, 0, 4, 2))  # Reduce left margin
plot.new()
legend("center", 
       legend = c("Not assigned", "RET", "SDHB", "VHL"),
       col = c("navy", "darkmagenta", "darkorange", "forestgreen"),
       pch = 16, 
       bty = "n")

# Reset the layout
layout(1)

# get the sPLS-DA data
loadings_comp1 <- final_splsda$loadings$X[,1] # For the first component
loadings_comp2 <- final_splsda$loadings$X[,2] # For the second component
loadings_all <- splsda$loadings$X
write.csv(loadings_comp1, file="splsda_loadings_comp1.csv")
write.csv(loadings_comp2, file="splsda_loadings_comp2.csv")
write.csv(loadings_all, file="splsda_loadings_all_wo_ctrls_perseus.csv")

legend2=list(legend = levels(ppglcat), # set of classes
            col = unique(color.mixo(ppglcat)), # set of colours
            title = "Sample", # legend title
            cex = 0.7) # legend size

cim <- cim(final_splsda, row.sideColors = color.mixo(ppglcat), 
           legend = legend2)

vip_scores <- vip(final_splsda)
write.csv(vip_scores, file="splsda_vip_wo_ctrls_perseus.csv")


library(dplyr)

summary_table <- data.frame(Variable = rownames(loadings_all),
                            Loading_Comp1 = loadings_all[,1],
                            Loading_Comp2 = loadings_all[,2],
                            VIP_Score = vip_scores) %>%
  arrange(desc(abs(Loading_Comp1))) # Sort by the magnitude of loading scores for component 1

write.csv(summary_table, file="slpsda_summary_wo_ctrls_perseus.csv")

# Calculate mean intensity across all samples
mean_intensity <- rowMeans(exprs, na.rm = TRUE)

# Keep proteins above a threshold (adjust based on your data distribution)
# Look at the distribution first:
hist(log10(mean_intensity), breaks = 50, main = "Distribution of Mean Intensities")

# Example: Keep top 50% or proteins above certain threshold
intensity_threshold <- quantile(mean_intensity, 0.5)  # Remove bottom 50%
high_quality_proteins <- mean_intensity > intensity_threshold

exprs_filtered <- exprs[high_quality_proteins, ]

library(preprocessCore)
exprs_quantile <- normalize.quantiles(as.matrix(exprs_filtered))
rownames(exprs_quantile) <- rownames(exprs_filtered)
colnames(exprs_quantile) <- colnames(exprs_filtered)

exprs_quantile_log2 <- log2(exprs_quantile + 1)
ppgldata_pca <- t(exprs_quantile_log2)

ppglpca <- pca(ppgldata_pca, 
                 ncomp = 2,
                 center = TRUE,
                 scale = FALSE,
                 max.iter = 500,
                 tol = 1e-09,
                 logratio = c("none", "CLR", "ILR"),
                 ilr.offset = 0.001,
                 V = NULL,
                 multilevel = NULL,
                 verbose.call = FALSE
)
#create the PCA plot:
layout(matrix(1))
plotIndiv(ppglpca, group = ppglcat,  legend = TRUE, title = 'PPGL TMT PCA',
          ellipse = FALSE, ind.names = FALSE, cex=1.2, point.lwd = 1.05,
          pch = 16, col = c("navy", "darkmagenta", "darkorange", "forestgreen"), style="graphics")

# Install/load pheatmap package
install.packages("pheatmap")
library(pheatmap)

# Z-score normalization by row (protein)
exprs_scaled <- t(apply(exprs_filtered_log2, 1, scale))
colnames(exprs_scaled) <- colnames(exprs_filtered_log2)

# Create annotation dataframe for pheatmap
annotation_col <- data.frame(Group = ppglcat)
rownames(annotation_col) <- colnames(exprs_scaled)

# Unsupervised hierarchical clustering heatmap
pheatmap(exprs_scaled, 
         annotation_col = annotation_col,
         show_rownames = FALSE,
         scale = "none",  # Already scaled above
         color = colorRampPalette(c("green", "black", "red"))(50),
         clustering_distance_rows = "minkowski",  # or "correlation"
         clustering_distance_cols = "minkowski",
         clustering_method = "ward.D2")  # or "average", "complete"
