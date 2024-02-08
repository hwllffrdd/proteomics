library(mixOmics)
setwd("R:/Documents/r/czppgl")
czppgldata <- read.csv("czppgl.csv", sep=";")
row.names(czppgldata) <- czppgldata$X
czppgldata$X <- NULL
#create a factor for sPLS-DA's "Y" values:
czppglcat <- factor(rep(c("pheo", "ctr"), times = c(13, 3))) 
#transverse the data, so that the individual samples' values are in rows:
czppgldata2 <- t(czppgldata)
czppgldata2_scaled <- scale(czppgldata2)
#perform sPLS-DA
czppglsplsda <- splsda(
  czppgldata2,
  czppglcat,
  ncomp = 2, keepX = c(50, 25)
)
summary(czppglsplsda)
#create the plot for initial sPLS-DA
plotIndiv(czppglsplsda, legend = TRUE, ellipse = TRUE, ind.names = FALSE, title = "CZPPGL sPLS-DA")
#optimize sPLS-DA:
plotVar(czppglsplsda, cutoff = 0.8)
list.keepX <- c(5:10,  seq(20, 100, 10))
tune.splsda._ <- tune.splsda(czppgldata2, czppglcat, ncomp = 3, 
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
plot(tune.splsda.srbct, col = color.jet(ncomp))
#create the final sPLS-DA plot with optimal parameters
final_splsda <- splsda(czppgldata2, czppglcat, ncomp = 2, keepX = select.keepX)
background <- background.predict(final_splsda, comp.predicted=2, dist = "centroids.dist")
background
czppglcat <- factor(czppglcat, levels = c("pheo", "ctr"))
plotIndiv(final_splsda, comp=c(1,2), legend=TRUE,
          ellipse = TRUE, title="CZPPGL sPLS-DA", ind.names = FALSE, pch=16:17,
          col=c("forestgreen", "red"), style="graphics")
plotLoadings(final_splsda, contrib = 'max', method = 'median', comp = 2, ndisplay = 20)
final_splsda$loadings

#using the same data, perform PCA:
czppglpca <- pca(czppgldata2, 
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
plotIndiv(czppglpca, group = czppglcat,  legend = TRUE, title = 'CZPPGL PCA',
          ellipse = TRUE, ind.names = FALSE, cex=1.2, point.lwd = 1.05,
          pch = 16:17, col = c("forestgreen", "red"), style="graphics")