# plotIndiv present in the package mixOmics
# methods(plotIndiv)
# plotIndiv.ipca*   plotIndiv.pca*    plotIndiv.pls*    plotIndiv.plsda*  plotIndiv.rcc*    plotIndiv.rgcca*  plotIndiv.sgcca*  plotIndiv.sipca*  plotIndiv.spls*   plotIndiv.splsda*

# New dependencies with ggplot2 and ellipse
# require(mixOmics); require(ellipse); require(ggplot2)
# source("../../mixOmics/R/plotIndiv.R"); 

#source("./R scripts/wrapper.sgcca.v6.R"); 
#source("./R scripts/wrappers.R"); source("./R scripts/helpers.R")


rm(list=ls())

require(mixOmics); require(ellipse); require(ggplot2)


sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# Florian
##sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R

# KA
sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R


#############
### IPCA ####
#############

if (plotIndiv.ipca){
data(liver.toxicity)
ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode="deflation")

# Style
plotIndiv(ipca.res)
plotIndiv(ipca.res, style = "lattice")
plotIndiv(ipca.res, style = "graphics")

# comp
plotIndiv(ipca.res, comp = c(1, 2));
plotIndiv(ipca.res, comp = c(1, 3));
plotIndiv(ipca.res, comp = c(2, 3))
plotIndiv(ipca.res, style = "lattice", comp = c(1, 2));
plotIndiv(ipca.res, style = "lattice", comp = c(1, 3));
plotIndiv(ipca.res, style = "lattice", comp = c(2, 3))
plotIndiv(ipca.res, style = "graphics", comp = c(1, 2))
plotIndiv(ipca.res, style = "graphics", comp = c(1, 3))
plotIndiv(ipca.res, style = "graphics", comp = c(2, 3))

# Ellipse
plotIndiv(ipca.res, plot.ellipse = TRUE)
plotIndiv(ipca.res, style = "lattice", plot.ellipse = TRUE)
plotIndiv(ipca.res, style = "graphics", plot.ellipse = TRUE)

# Ellipse & confidence
plotIndiv(ipca.res, ellipse.level = 0.5, plot.ellipse = TRUE)
plotIndiv(ipca.res, style = "lattice", ellipse.level = 0.5, plot.ellipse = TRUE)
plotIndiv(ipca.res, style = "graphics", ellipse.level = 0.5, plot.ellipse = TRUE)

# Legend
plotIndiv(ipca.res, add.legend = TRUE)
plotIndiv(ipca.res, style = "lattice", add.legend = TRUE)
plotIndiv(ipca.res, style = "graphics", add.legend = TRUE)

# Abline.line
plotIndiv(ipca.res, abline.line = TRUE)
plotIndiv(ipca.res, style = "lattice", abline.line = TRUE)
plotIndiv(ipca.res, style = "graphics", abline.line = TRUE)

# X.label
plotIndiv(ipca.res, X.label = "I am here")
plotIndiv(ipca.res, style = "lattice", X.label = "I am here")
plotIndiv(ipca.res, style = "graphics", X.label = "I am here")

# Y.label
plotIndiv(ipca.res, Y.label = "I am here")
plotIndiv(ipca.res, style = "lattice", Y.label = "I am here")
plotIndiv(ipca.res, style = "graphics", Y.label = "I am here")

# Color
plotIndiv(ipca.res, col = "purple")
plotIndiv(ipca.res, style = "lattice", col = "purple")
plotIndiv(ipca.res, style = "graphics", col = "purple")

# pch
plotIndiv(ipca.res, pch = 2, ind.names = FALSE)
plotIndiv(ipca.res, style = "lattice", pch = 2, ind.names = FALSE)
plotIndiv(ipca.res, style = "graphics", pch = 2, ind.names = FALSE)

plotIndiv(ipca.res, pch = 2, ind.names = FALSE, add.legend = TRUE)
plotIndiv(ipca.res, style = "lattice", pch = 2, ind.names = FALSE, add.legend = TRUE)
plotIndiv(ipca.res, style = "graphics", pch = 2, ind.names = FALSE, add.legend = TRUE)

# cex
plotIndiv(ipca.res, cex = 2, ind.names = FALSE)
plotIndiv(ipca.res, style = "lattice", cex = 2, ind.names = FALSE)
plotIndiv(ipca.res, style = "graphics", cex = 2, ind.names = FALSE)

plotIndiv(ipca.res, cex = 2, ind.names = FALSE, add.legend = TRUE)
plotIndiv(ipca.res, style = "lattice", cex = 2, ind.names = FALSE, add.legend = TRUE)
plotIndiv(ipca.res, style = "graphics", cex = 2, ind.names = FALSE, add.legend = TRUE)

# rep.space = "X-variate"
plotIndiv(ipca.res, rep.space = "X-variate")
plotIndiv(ipca.res, style = "lattice", rep.space = "X-variate")
plotIndiv(ipca.res, style = "graphics", rep.space = "X-variate")

plotIndiv(ipca.res, rep.space = "XY-variate")
plotIndiv(ipca.res, style = "lattice", rep.space = "XY-variate")
plotIndiv(ipca.res, style = "graphics", rep.space = "XY-variate")

plotIndiv(ipca.res, rep.space = "Y-variate")
plotIndiv(ipca.res, style = "lattice", rep.space = "Y-variate")
plotIndiv(ipca.res, style = "graphics", rep.space = "Y-variate")

# block
plotIndiv(ipca.res, blocks = "Y")
plotIndiv(ipca.res, style = "lattice", blocks = "Y")
plotIndiv(ipca.res, style = "graphics", blocks = "Y")

# Test group.training
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 5, group.training = liver.toxicity$treatment[, 4])
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 5, group.training = liver.toxicity$treatment[, 4], add.legend = TRUE)
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 5, group.training = liver.toxicity$treatment[, 4], add.legend = TRUE, col = c("pink", "blue", "green", "black"))
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 5, group = liver.toxicity$treatment[, 4], plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group.training = liver.toxicity$treatment[, 4], style = "lattice")
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group.training = liver.toxicity$treatment[, 4], style = "lattice", add.legend = TRUE)
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group.training = liver.toxicity$treatment[, 4], style = "lattice", add.legend = TRUE, col = c("pink", "blue", "green", "black"))
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = liver.toxicity$treatment[, 4], style = "lattice", plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group.training = liver.toxicity$treatment[, 4], style = "graphics")
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group.training = liver.toxicity$treatment[, 4], style = "graphics", add.legend = TRUE)
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group.training = liver.toxicity$treatment[, 4], style = "graphics", add.legend = TRUE, col = c("pink", "blue", "green", "black"))
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = liver.toxicity$treatment[, 4], style = "graphics", plot.ellipse = TRUE)

# Test ind.names
plotIndiv(ipca.res, ind.names = FALSE, cex = 5, group.training = liver.toxicity$treatment[, 4])
plotIndiv(ipca.res, ind.names = FALSE, cex = 5, group = liver.toxicity$treatment[, 4], plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, cex = 1, group.training = liver.toxicity$treatment[, 4], style = "lattice")
plotIndiv(ipca.res, ind.names = FALSE, cex = 1, group = liver.toxicity$treatment[, 4], style = "lattice", plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, cex = 1, group.training = liver.toxicity$treatment[, 4], style = "graphics")
plotIndiv(ipca.res, ind.names = FALSE, cex = 1, group = liver.toxicity$treatment[, 4], style = "graphics", plot.ellipse = TRUE)

# Test cex
plotIndiv(ipca.res, ind.names = FALSE, cex = c(3,8,13,18), group.training = liver.toxicity$treatment[, 4])
plotIndiv(ipca.res, ind.names = FALSE, cex = c(3,8,13,18), group = liver.toxicity$treatment[, 4], plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, cex = c(3,8,13,18), group = liver.toxicity$treatment[, 4], plot.ellipse = TRUE, add.legend = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, cex = c(1, 2, 3, 4), group.training = liver.toxicity$treatment[, 4], style = "lattice")
plotIndiv(ipca.res, ind.names = FALSE, cex = c(1, 2, 3, 4), group = liver.toxicity$treatment[, 4], style = "lattice", plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, cex = c(1, 2, 3, 4), group = liver.toxicity$treatment[, 4], style = "lattice", plot.ellipse = TRUE, add.legend = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, cex = c(1, 2, 3, 4), group.training = liver.toxicity$treatment[, 4], style = "graphics")
plotIndiv(ipca.res, ind.names = FALSE, cex = c(1, 2, 3, 4), group = liver.toxicity$treatment[, 4], style = "graphics", plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, cex = c(1, 2, 3, 4), group = liver.toxicity$treatment[, 4], style = "graphics", plot.ellipse = TRUE, add.legend = TRUE)

# Test pch
plotIndiv(ipca.res, ind.names = FALSE, pch = c(3,8,13,18), group.training = liver.toxicity$treatment[, 4])
plotIndiv(ipca.res, ind.names = FALSE, pch = c(3,8,13,18), group = liver.toxicity$treatment[, 4], plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, pch = c(3,8,13,18), group.training = liver.toxicity$treatment[, 4], style = "lattice")
plotIndiv(ipca.res, ind.names = FALSE, pch = c(3,8,13,18), group = liver.toxicity$treatment[, 4], style = "lattice", plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, pch = c(3,8,13,18), group.training = liver.toxicity$treatment[, 4], style = "graphics")
plotIndiv(ipca.res, ind.names = FALSE, pch = c(3,8,13,18), group = liver.toxicity$treatment[, 4], style = "graphics", plot.ellipse = TRUE)

# abline.line
plotIndiv(ipca.res, ind.names = FALSE, abline = TRUE, group.training = liver.toxicity$treatment[, 4])
plotIndiv(ipca.res, ind.names = FALSE, abline = TRUE, group = liver.toxicity$treatment[, 4], plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, abline = TRUE, group.training = liver.toxicity$treatment[, 4], style = "lattice")
plotIndiv(ipca.res, ind.names = FALSE, abline = TRUE, group = liver.toxicity$treatment[, 4], style = "lattice", plot.ellipse = TRUE)
plotIndiv(ipca.res, ind.names = FALSE, abline = TRUE, group.training = liver.toxicity$treatment[, 4], style = "graphics")
plotIndiv(ipca.res, ind.names = FALSE, abline = TRUE, group = liver.toxicity$treatment[, 4], style = "graphics", plot.ellipse = TRUE)
}

data(liver.toxicity)
sipca.res <- sipca(liver.toxicity$gene, ncomp = 3, mode="deflation", keepX=c(50,50,50))
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 5, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])))
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 5, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), plot.ellipse = TRUE)
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), style = "lattice")
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), style = "lattice", plot.ellipse = TRUE)
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), style = "graphics")
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), style = "graphics", plot.ellipse = TRUE)

############
### PCA ####
############

if (plotIndiv.pca){
data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)

data(liver.toxicity)
ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode="deflation")

# Style
plotIndiv(pca.res)
plotIndiv(pca.res, style = "lattice")
plotIndiv(pca.res, style = "graphics")

# comp
plotIndiv(pca.res, comp = c(1, 2));
plotIndiv(pca.res, comp = c(1, 3));
plotIndiv(pca.res, comp = c(2, 3))
plotIndiv(pca.res, style = "lattice", comp = c(1, 2));
plotIndiv(pca.res, style = "lattice", comp = c(1, 3));
plotIndiv(pca.res, style = "lattice", comp = c(2, 3))
plotIndiv(pca.res, style = "graphics", comp = c(1, 2))
plotIndiv(pca.res, style = "graphics", comp = c(1, 3))
plotIndiv(pca.res, style = "graphics", comp = c(2, 3))

# Ellipse
plotIndiv(pca.res, plot.ellipse = TRUE)
plotIndiv(pca.res, style = "lattice", plot.ellipse = TRUE)
plotIndiv(pca.res, style = "graphics", plot.ellipse = TRUE)

# Ellipse & confidence
plotIndiv(pca.res, ellipse.level = 0.5, plot.ellipse = TRUE)
plotIndiv(pca.res, style = "lattice", ellipse.level = 0.5, plot.ellipse = TRUE)
plotIndiv(pca.res, style = "graphics", ellipse.level = 0.5, plot.ellipse = TRUE)

# Legend
plotIndiv(pca.res, add.legend = TRUE)
plotIndiv(pca.res, style = "lattice", add.legend = TRUE)
plotIndiv(pca.res, style = "graphics", add.legend = TRUE)

# Abline.line
plotIndiv(pca.res, abline.line = TRUE)
plotIndiv(pca.res, style = "lattice", abline.line = TRUE)
plotIndiv(pca.res, style = "graphics", abline.line = TRUE)

# X.label
plotIndiv(pca.res, X.label = "I am here")
plotIndiv(pca.res, style = "lattice", X.label = "I am here")
plotIndiv(pca.res, style = "graphics", X.label = "I am here")

# Y.label
plotIndiv(pca.res, Y.label = "I am here")
plotIndiv(pca.res, style = "lattice", Y.label = "I am here")
plotIndiv(pca.res, style = "graphics", Y.label = "I am here")

# Color
plotIndiv(pca.res, col = "purple")
plotIndiv(pca.res, style = "lattice", col = "purple")
plotIndiv(pca.res, style = "graphics", col = "purple")

# pch
plotIndiv(pca.res, pch = 2, ind.names = FALSE)
plotIndiv(pca.res, style = "lattice", pch = 2, ind.names = FALSE)
plotIndiv(pca.res, style = "graphics", pch = 2, ind.names = FALSE)

plotIndiv(pca.res, pch = 2, ind.names = FALSE, add.legend = TRUE)
plotIndiv(pca.res, style = "lattice", pch = 2, ind.names = FALSE, add.legend = TRUE)
plotIndiv(pca.res, style = "graphics", pch = 2, ind.names = FALSE, add.legend = TRUE)

# cex
plotIndiv(pca.res, cex = 2, ind.names = FALSE)
plotIndiv(pca.res, style = "lattice", cex = 2, ind.names = FALSE)
plotIndiv(pca.res, style = "graphics", cex = 2, ind.names = FALSE)

plotIndiv(pca.res, cex = 2, ind.names = FALSE, add.legend = TRUE)
plotIndiv(pca.res, style = "lattice", cex = 2, ind.names = FALSE, add.legend = TRUE)
plotIndiv(pca.res, style = "graphics", cex = 2, ind.names = FALSE, add.legend = TRUE)

# rep.space = "X-variate"
plotIndiv(pca.res, rep.space = "X-variate")
plotIndiv(pca.res, style = "lattice", rep.space = "X-variate")
plotIndiv(pca.res, style = "graphics", rep.space = "X-variate")

plotIndiv(pca.res, rep.space = "XY-variate")
plotIndiv(pca.res, style = "lattice", rep.space = "XY-variate")
plotIndiv(pca.res, style = "graphics", rep.space = "XY-variate")

plotIndiv(pca.res, rep.space = "Y-variate")
plotIndiv(pca.res, style = "lattice", rep.space = "Y-variate")
plotIndiv(pca.res, style = "graphics", rep.space = "Y-variate")

# block
plotIndiv(pca.res, blocks = "Y")
plotIndiv(pca.res, style = "lattice", blocks = "Y")
plotIndiv(pca.res, style = "graphics", blocks = "Y")

# Test group.training
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 5, col =  as.numeric(as.factor(multidrug$cell.line$Class)))
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 5, group =  as.numeric(as.factor(multidrug$cell.line$Class)))
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 5, group =  as.numeric(as.factor(multidrug$cell.line$Class)), add.legend = TRUE)

plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 5, group.training = as.numeric(as.factor(multidrug$cell.line$Class)))
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 5, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), add.legend = TRUE)
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 5, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), add.legend = TRUE, col = 1:9)
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 5, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice")
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice", add.legend = TRUE)
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice", add.legend = TRUE, col = 1:9)
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice", plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics")
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics", add.legend = TRUE)
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics", add.legend = TRUE, col = 1:9)
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics", plot.ellipse = TRUE)

# Test ind.names
plotIndiv(pca.res, ind.names = FALSE, cex = 5, group.training = as.numeric(as.factor(multidrug$cell.line$Class)))
plotIndiv(pca.res, ind.names = FALSE, cex = 5, group = as.numeric(as.factor(multidrug$cell.line$Class)), plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = FALSE, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice")
plotIndiv(pca.res, ind.names = FALSE, cex = 1, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice", plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = FALSE, cex = 1, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics")
plotIndiv(pca.res, ind.names = FALSE, cex = 1, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics", plot.ellipse = TRUE)

# Test cex
plotIndiv(pca.res, ind.names = FALSE, cex = 1:9, group.training = as.numeric(as.factor(multidrug$cell.line$Class)))
plotIndiv(pca.res, ind.names = FALSE, cex = 1:9, group = as.numeric(as.factor(multidrug$cell.line$Class)), plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = FALSE, cex = 1:9, group = as.numeric(as.factor(multidrug$cell.line$Class)), plot.ellipse = TRUE, add.legend = TRUE)
plotIndiv(pca.res, ind.names = FALSE, cex = 1:9, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice")
plotIndiv(pca.res, ind.names = FALSE, cex = 1:9, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice", plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = FALSE, cex = 1:9, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice", plot.ellipse = TRUE, add.legend = TRUE)
plotIndiv(pca.res, ind.names = FALSE, cex = 1:9, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics")
plotIndiv(pca.res, ind.names = FALSE, cex = 1:9, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics", plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = FALSE, cex = 1:9, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics", plot.ellipse = TRUE, add.legend = TRUE)

# Test pch
plotIndiv(pca.res, ind.names = FALSE, pch = 10:18, group = as.numeric(as.factor(multidrug$cell.line$Class)))
plotIndiv(pca.res, ind.names = FALSE, pch = 10:18, group = as.numeric(as.factor(multidrug$cell.line$Class)), plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = FALSE, pch = 10:18, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice")
plotIndiv(pca.res, ind.names = FALSE, pch = 10:18, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice", plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = FALSE, pch = 10:18, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics")
plotIndiv(pca.res, ind.names = FALSE, pch = 10:18, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics", plot.ellipse = TRUE)

# abline.line
plotIndiv(pca.res, ind.names = FALSE, abline = TRUE, group.training = as.numeric(as.factor(multidrug$cell.line$Class)))
plotIndiv(pca.res, ind.names = FALSE, abline = TRUE, group = as.numeric(as.factor(multidrug$cell.line$Class)), plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = FALSE, abline = TRUE, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice")
plotIndiv(pca.res, ind.names = FALSE, abline = TRUE, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "lattice", plot.ellipse = TRUE)
plotIndiv(pca.res, ind.names = FALSE, abline = TRUE, group.training = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics")
plotIndiv(pca.res, ind.names = FALSE, abline = TRUE, group = as.numeric(as.factor(multidrug$cell.line$Class)), style = "graphics", plot.ellipse = TRUE)
}

data(liver.toxicity)
spca.rat <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))
plotIndiv(spca.rat, ind.names = liver.toxicity$treatment[, 3], cex = 5,  group = as.numeric(liver.toxicity$treatment[, 3]))
plotIndiv(spca.rat, ind.names = liver.toxicity$treatment[, 3], cex = 1,  group = as.numeric(liver.toxicity$treatment[, 3]), style = "lattice")
plotIndiv(spca.rat, ind.names = liver.toxicity$treatment[, 3], cex = 1,  group = as.numeric(liver.toxicity$treatment[, 3]), style = "graphics")


#############
### RGCCA ###
#############

if (plotIndiv.rgcca){
  
  data(nutrimouse)
  Y = unmap(nutrimouse$diet)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
  design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
  nutrimouse.rgcca <- wrapper.rgcca(blocks = data,
                                    design = design1,
                                    tau = "optimal",
                                    ncomp = c(2, 2, 2),
                                    scheme = "centroid",
                                    verbose = FALSE)
  
  nutrimouse.rgcca
  
  plotIndiv(nutrimouse.rgcca)
  plotIndiv(nutrimouse.rgcca, style = "lattice")
  plotIndiv(nutrimouse.rgcca, style = "graphics")
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet)
  plotIndiv(nutrimouse.rgcca, style = "lattice", group = nutrimouse$diet)
  plotIndiv(nutrimouse.rgcca, style = "graphics", group = nutrimouse$diet)
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, add.legend =TRUE)
  plotIndiv(nutrimouse.rgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE)
  plotIndiv(nutrimouse.rgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE)
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, add.legend =TRUE, plot.ellipse = TRUE)
  plotIndiv(nutrimouse.rgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE)
  plotIndiv(nutrimouse.rgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE)
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, add.legend =TRUE, plot.ellipse = TRUE, ellipse.level = 0.5)
  plotIndiv(nutrimouse.rgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5)
  plotIndiv(nutrimouse.rgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5)
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, add.legend =TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = 2)
  plotIndiv(nutrimouse.rgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = 2)
  plotIndiv(nutrimouse.rgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = 2)
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, add.legend =TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid")
  plotIndiv(nutrimouse.rgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid")
  plotIndiv(nutrimouse.rgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid")
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, cex = c(1,2,3,4,5))
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", cex = c(1,2,3,4,5))
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", cex = c(1,2,3,4,5))
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, cex = c(1,2,3,4,5), blocks = "gene")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", cex = c(1,2,3,4,5), blocks = "gene")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, tyle = "graphics", cex = c(1,2,3,4,5), blocks = "gene")
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, pch = c(9:13), ind.names = FALSE, add.legend = TRUE)
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", pch = c(9:13), ind.names = FALSE, add.legend = TRUE)
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", pch = c(9:13), ind.names = FALSE, add.legend = TRUE)
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, pch = c(9:13), blocks = "Y", ind.names = FALSE)
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", pch = c(9:13), blocks = "Y", ind.names = FALSE)
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", pch = c(1,2,3,4,5), blocks = "Y", ind.names = FALSE)
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, X.label = "I am here")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", X.label = "I am here")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", X.label = "I am here")
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, Y.label = "I am here")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", Y.label = "I am here")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", Y.label = "I am here")
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, main = "I am here")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", main = "I am here")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", main = "I am here")
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, rep.space = "X-variate")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", rep.space = "X-variate")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", rep.space = "X-variate")
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, rep.space = "Y-variate")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", rep.space = "Y-variate")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", rep.space = "Y-variate")
  
  # KA: from what I understand, this option is not available anymore:
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, rep.space = "XY-variate")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", rep.space = "XY-variate")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", rep.space = "XY-variate")
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, col = c(10:14), add.legend = TRUE)
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", col = c(10:14), add.legend = TRUE)
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "graphics", col = c(10:14), add.legend = TRUE)
  
  plotIndiv(nutrimouse.rgcca, col = "pink", add.legend = TRUE)
  plotIndiv(nutrimouse.rgcca, style = "lattice", col = "pink", add.legend = TRUE)
  plotIndiv(nutrimouse.rgcca, style = "graphics", col = "pink", , add.legend = TRUE)
  
  plotIndiv(nutrimouse.rgcca, col = "pink", add.legend = TRUE, abline.line = TRUE)
  plotIndiv(nutrimouse.rgcca, style = "lattice", col = "pink", add.legend = TRUE, abline.line = TRUE)
  plotIndiv(nutrimouse.rgcca, style = "graphics", col = "pink", , add.legend = TRUE, abline.line = TRUE)
  
  plotIndiv(nutrimouse.rgcca, comp = c(2, 1))
  plotIndiv(nutrimouse.rgcca, style = "lattice", comp = c(2, 1))
  plotIndiv(nutrimouse.rgcca, style = "graphics", comp = c(2, 1))
  
  nutrimouse.rgcca <- wrapper.rgcca(blocks = data,
                                    design = design1,
                                    tau = "optimal",
                                    ncomp = c(2, 2, 1),
                                    scheme = "centroid",
                                    verbose = FALSE)
  
  nutrimouse.rgcca
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, blocks = "gene")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", blocks = "gene")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, tyle = "graphics", blocks = "gene")
  
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, blocks = "lipid")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, style = "lattice", blocks = "lipid")
  plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, tyle = "graphics", blocks = "lipid")  
}


#############
### SGCCA ###
#############

if (plotIndiv.sgcca){
  
  data(nutrimouse)
  Y = unmap(nutrimouse$diet)
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
  design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
  nutrimouse.sgcca <- wrapper.sgcca(blocks = data,
                                    design = design1,
                                    penalty = c(0.3, 0.5, 1),
                                    ncomp = c(2, 2, 3),
                                    scheme = "centroid",
                                    verbose = FALSE, 
                                    bias = FALSE)
  
  nutrimouse.sgcca
  
  plotIndiv(nutrimouse.sgcca)
  plotIndiv(nutrimouse.sgcca, style = "lattice")
  plotIndiv(nutrimouse.sgcca, style = "graphics")
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet)
  plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet)
  plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet)
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE)
  plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE)
  plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE)
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE, plot.ellipse = TRUE)
  plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE)
  plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE)
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE, plot.ellipse = TRUE, ellipse.level = 0.5)
  plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5)
  plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5)
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = 2)
  plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = 2)
  plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = 2)
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid")
  plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid")
  plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid")
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, cex = c(1,2,3,4,5))
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", cex = c(1,2,3,4,5))
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", cex = c(1,2,3,4,5))
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, cex = c(1,2,3,4,5), blocks = "gene")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", cex = c(1,2,3,4,5), blocks = "gene")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, tyle = "graphics", cex = c(1,2,3,4,5), blocks = "gene")
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, pch = c(9:13), ind.names = FALSE, add.legend = TRUE)
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", pch = c(9:13), ind.names = FALSE, add.legend = TRUE)
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", pch = c(9:13), ind.names = FALSE, add.legend = TRUE)
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, pch = c(9:13), blocks = "Y", ind.names = FALSE)
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", pch = c(9:13), blocks = "Y", ind.names = FALSE)
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", pch = c(1,2,3,4,5), blocks = "Y", ind.names = FALSE)
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, X.label = "I am here")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", X.label = "I am here")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", X.label = "I am here")
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, Y.label = "I am here")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", Y.label = "I am here")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", Y.label = "I am here")
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, main = "I am here")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", main = "I am here")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", main = "I am here")
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, rep.space = "X-variate")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", rep.space = "X-variate")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", rep.space = "X-variate")
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, rep.space = "Y-variate")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", rep.space = "Y-variate")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", rep.space = "Y-variate")
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, rep.space = "XY-variate")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", rep.space = "XY-variate")
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", rep.space = "XY-variate")
  
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, col = c(10:14), add.legend = TRUE)
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "lattice", col = c(10:14), add.legend = TRUE)
  plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, style = "graphics", col = c(10:14), add.legend = TRUE)
  
  plotIndiv(nutrimouse.sgcca, col = "pink", add.legend = TRUE)
  plotIndiv(nutrimouse.sgcca, style = "lattice", col = "pink", add.legend = TRUE)
  plotIndiv(nutrimouse.sgcca, style = "graphics", col = "pink", , add.legend = TRUE)
  
  plotIndiv(nutrimouse.sgcca, col = "pink", add.legend = TRUE, abline.line = TRUE)
  plotIndiv(nutrimouse.sgcca, style = "lattice", col = "pink", add.legend = TRUE, abline.line = TRUE)
  plotIndiv(nutrimouse.sgcca, style = "graphics", col = "pink", , add.legend = TRUE, abline.line = TRUE)
  
  plotIndiv(nutrimouse.sgcca, comp = c(2, 1))
  plotIndiv(nutrimouse.sgcca, style = "lattice", comp = c(2, 1))
  plotIndiv(nutrimouse.sgcca, style = "graphics", comp = c(2, 1))

  plotIndiv(nutrimouse.sgcca, comp = c(1, 3), blocks = "Y")
  plotIndiv(nutrimouse.sgcca, style = "lattice", comp = c(1, 3), blocks = "Y")
  plotIndiv(nutrimouse.sgcca, style = "graphics", comp = c(1, 3), blocks = "Y")
  
  nutrimouse.sgcca <- wrapper.sgcca(blocks = data,
                                    design = design1,
                                    penalty = c(0.3, 0.5, 1),
                                    ncomp = c(2, 2, 1),
                                    scheme = "centroid",
                                    verbose = FALSE, 
                                    bias = FALSE)
  
  
  plotIndiv(nutrimouse.sgcca, add.legend = FALSE, blocks = "gene", ind.names = FALSE, group = nutrimouse$diet)
  plotIndiv(nutrimouse.sgcca, style = "lattice", add.legend = FALSE, blocks = "gene", ind.names = FALSE, group = nutrimouse$diet)
  plotIndiv(nutrimouse.sgcca, style = "graphics", add.legend = FALSE, blocks = "gene", ind.names = FALSE, group = nutrimouse$diet) 
  
  plotIndiv(nutrimouse.sgcca, add.legend = FALSE, blocks = "lipid", ind.names = FALSE, group = nutrimouse$diet)
  plotIndiv(nutrimouse.sgcca, style = "lattice", add.legend = FALSE, blocks = "lipid", ind.names = FALSE, group = nutrimouse$diet)
  plotIndiv(nutrimouse.sgcca, style = "graphics", add.legend = FALSE, blocks = "lipid", ind.names = FALSE, group = nutrimouse$diet) 
  
  s.match(nutrimouse.sgcca$variates$gene, nutrimouse.sgcca$variates$lipid, col = as.numeric(nutrimouse$diet))
}  
  
  

data(liver.toxicity)
sipca.res <- sipca(liver.toxicity$gene, ncomp = 3, mode="deflation", keepX=c(50,50,50))
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 5, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])))
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 5, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), plot.ellipse = TRUE)
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), style = "lattice")
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), style = "lattice", plot.ellipse = TRUE)
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), style = "graphics")
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 1, group = as.numeric(as.factor(liver.toxicity$treatment[, 4])), style = "graphics", plot.ellipse = TRUE)

data(multidrug)

## this data set contains missing values, therefore 
## the 'prcomp' function cannot be applied
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)

plotIndiv(pca.res, group.training = c(rep("group1", 20), rep("group2", 20), rep("group3", 20)), style = "ggplot2", add.legend = FALSE, main = "main")
plotIndiv(pca.res, pch = c(15, 16, 17), add.legend = FALSE, ind.names = FALSE, col = c("red", "green", "black"),
          style = "graphics", cex = c(1, 1, 1), group.training = c(rep("group1", 20), rep("group2", 20), rep("group3", 20)))


plotIndiv(pca.res, ind.names = FALSE, comp = c(1, 3), plot.ellipse = TRUE, add.legend = TRUE, group.training = c(rep("group1", 20), rep("group2", 20), rep("group3", 20)))
data(liver.toxicity)

# implement IPCA on a microarray dataset
ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode="deflation")
ipca.res

plotIndiv(ipca.res, ind.names = TRUE, comp = c(1, 3), plot.ellipse = TRUE, add.legend = TRUE, group.training = c(rep("group1", 14), rep("group2", 50)))#, group.training = c(rep("group1", 10), rep("group2", 50)))#, group.training = c(rep("group1", 10), rep("group2", 50)))#, cex = 7, pch = 2, group.training = c(rep("group1", 10), rep("group2", 50)))
plotIndiv(pca.res, ind.names = FALSE, cex = c(rep(1, 60),  group.training = c(rep("group1", 14), rep("group2", 50)))

group.training = factor(c(rep("group1", 10), rep("group2", 50)))
col = c(1, 2)

plotIndiv(pca.res, ind.names = TRUE, pch = c(1, 3), plot.ellipse = TRUE, group.training = c(rep("group1", 10), rep("group2", 50)))
  plotIndiv(pca.res, ind.names = FALSE, pch = 1, plot.ellipse = FALSE, cex = rep(1, 60))


plotIndiv(pca.res, plot.ellipse = FALSE, ind.names = FALSE, group.training = c(rep("group1", 10), rep("group2", 50)), cex = c(1, 5))


data(linnerud)
X1 <- linnerud$exercise
Y2 <- linnerud$physiological
linn.res <- rcc(X1, Y2, ncomp = 3)

plotIndiv(linn.res, ellipse.level = 0.5, rep.space = "X-variate")

plotIndiv(linn.res, plot.ellipse = FALSE, comp = c(2, 3))


data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
plotIndiv(nutri.res, ellipse.level = 0.9, group.training = c(rep("group1", 20), rep("group2", 20)), style = "graphics", add.legend = FALSE, blocks = "Y")


data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, gene expression and lipids are connected to the diet factor
# design = matrix(c(0,0,1,
#                   0,0,1,
#                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#note: the tau parameter is the regularization parameter
wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
plotIndiv(wrap.result.rgcca)



#############
### plsda ###
#############

data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 3)

# with pch and ggplot2 style - note: no need to input the group color as it is done internally in a supervised approach
plotIndiv(plsda.breast, comp = c(1, 2),  pch = 17, ind.names = FALSE,
          rep.space = "XY-variate", plot.ellipse = TRUE, style = "ggplot2")

# with ind names and lattice style
plotIndiv(plsda.breast, ind.names = TRUE, comp = c(1, 2), plot.indiv = FALSE, 
          rep.space = "XY-variate", plot.ellipse = TRUE, style = "lattice", cex = c(1, 1))


#############
### splsda ###
#############

data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

splsda.breast <- splsda(X, Y,keepX=c(10,10),ncomp=2)

# default option
plotIndiv(splsda.breast)

# default option with no ind name: pch and color are set automatically
plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2))

# default option with no ind name: pch and color are set automatically, with legend
plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2), add.legend = TRUE)


plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), plot.indiv = FALSE, rep.space = "XY-variate", plot.ellipse = TRUE, style = "ggplot2", cex = c(1, 1))
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), plot.indiv = FALSE, rep.space = "XY-variate", plot.ellipse = TRUE, style = "lattice", cex = c(1, 1))

plotIndiv(splsda.breast,comp=c(2,3))#erreur-ok

