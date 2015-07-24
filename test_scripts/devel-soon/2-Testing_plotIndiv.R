# -----------------------------------------------------------------------------------
# Testing_plotIndiv.R
# Author:    KA Le Cao 
# Date started:  23/07/2015
# Last updated:  
# Objective: only code for help file
# Latest update: 
# -----------------------------------------------------------------------------------


rm(list=ls())

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# Florian
sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R

# KA
sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R



## plot of individuals for objects of class 'rcc' 
# ----------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# default, only in the X space
plotIndiv(nutri.res) 

# ellipse with respect to genotype in the XY space, names also indicate genotype
plotIndiv(nutri.res, rep.space= 'XY-variate', plot.ellipse = TRUE, ellipse.level = 0.9, group = nutrimouse$genotype, ind.names = nutrimouse$genotype)

# ellipse with respect to genotype in the XY space, with legend
plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype, add.legend = TRUE)


# lattice style
plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype, add.legend = TRUE, style = 'lattice')

# classic style, in the Y space
plotIndiv(nutri.res, rep.space= 'Y-variate', group = nutrimouse$genotype, add.legend = TRUE, style = 'graphics')


## plot of individuals for objects of class 'pls' or 'spls'  
# ----------------------------------------------------   
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

#default
plotIndiv(toxicity.spls)

# in the Y space, colors indicate time of necropsy, text is the dose
plotIndiv(toxicity.spls, rep.space= 'Y-variate', group = liver.toxicity$treatment[, 'Time.Group'], ind.names = liver.toxicity$treatment[, 'Dose.Group'], add.legend = TRUE)


## plot of individuals for objects of class 'plsda' or 'splsda'  
# ----------------------------------------------------   
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

splsda.breast <- splsda(X, Y,keepX=c(10,10),ncomp=2)

# default option: note the outcome color is included by default!
plotIndiv(splsda.breast)

# default option with no ind name: pch and color are set automatically
plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2))

# default option with no ind name: pch and color are set automatically, with legend
plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2), add.legend = TRUE)

# playing with style
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), plot.indiv = FALSE, plot.ellipse = TRUE, style = "ggplot2", cex = c(1, 1))
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), plot.indiv = FALSE, plot.ellipse = TRUE, style = "lattice", cex = c(1, 1))


## variable representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------
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

# default style: one panel for each block
plotIndiv(nutrimouse.sgcca)

# for the block 'lipid' with ellipse plots and legend, different styles
plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", main = 'my plot')
plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", main = 'my plot')
plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, add.legend = TRUE, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", main = 'my plot')


## variable representation for objects of class 'sgccda' 
# ----------------------------------------------------
# Note: the code differs from above as we use a 'supervised' GCCA analysis
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)

nutrimouse.sgccda1 <- wrapper.sgccda(blocks = data,
                                     Y = Y,
                                     design = design1,
                                     ncomp = c(2, 2),
                                     keep.blocks = list(c(10,10), c(15,15)),
                                     scheme = "centroid",
                                     verbose = FALSE,
                                     bias = FALSE)


# plotIndiv
# ----------
# displaying all blocks. bu default colors correspond to outcome Y
plotIndiv(nutrimouse.sgccda1)

# displaying only 2 blocks
plotIndiv(nutrimouse.sgccda1, blocks = c(1,2), group = nutrimouse$diet)

# with some ellipse, legend and title
plotIndiv(nutrimouse.sgccda1, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE, add.legend = TRUE, main = 'my sample plot')


