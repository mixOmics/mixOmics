# -----------------------------------------------------------------------------------
# Testing_plotIndiv.R
# Author:    KA Le Cao 
# Date started:  23/07/2015
# Last updated:  
# Objective: only code for help file
# Latest update: 
# -----------------------------------------------------------------------------------


## note for later patch!!: add.legend not looking great on lattice or graphics

rm(list=ls())

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# Florian
##sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R

library(mixOmics, lib.loc = 'MyR/')

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
#changing the colors with argument col and ellipse will be plotted according to the color
plotIndiv(nutri.res, col= as.numeric(nutrimouse$diet), plot.ellipse = TRUE)

# or we can specify the argument group for plotting the ellipse according to group
plotIndiv(nutri.res, col= as.numeric(nutrimouse$diet), 
          plot.ellipse = TRUE, group = nutrimouse$genotype)


# plotting the samples in the XY space, with names indicating genotype
plotIndiv(nutri.res, rep.space= 'XY-variate', plot.ellipse = TRUE, ellipse.level = 0.9, 
          group = nutrimouse$genotype, ind.names = nutrimouse$genotype)

# ellipse with respect to genotype in the XY space, with legend according to group argument
plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype, add.legend = TRUE)


# lattice style, with legend according to group argument
plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype, 
          style = 'lattice')

# classic style, in the Y space
plotIndiv(nutri.res, rep.space= 'Y-variate', group = nutrimouse$genotype, 
          style = 'graphics')



# note: legend not looking that great for those styles, omitted from help file
# lattice style, with legend according to group argument
plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype, 
          add.legend = TRUE, style = 'lattice')

# classic style, in the Y space
plotIndiv(nutri.res, rep.space= 'Y-variate', group = nutrimouse$genotype, 
          add.legend = TRUE, style = 'graphics')


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
plotIndiv(toxicity.spls, rep.space= 'Y-variate', group = liver.toxicity$treatment[, 'Time.Group'], 
          ind.names = liver.toxicity$treatment[, 'Dose.Group'], add.legend = TRUE)

# in the Y space, colors indicate time of necropsy, text is the dose, 
# changing the color per group, ellipse plots
plotIndiv(toxicity.spls, rep.space= 'Y-variate', group = liver.toxicity$treatment[, 'Time.Group'], 
          ind.names = liver.toxicity$treatment[, 'Dose.Group'], add.legend = TRUE,
          col.per.group = c(1:4), plot.ellipse = TRUE)


## plot of individuals for objects of class 'plsda' or 'splsda'  
# ----------------------------------------------------   
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

splsda.breast <- splsda(X, Y,keepX=c(10,10),ncomp=2)

# default option: note the outcome color is included by default as it is a supervised approach
plotIndiv(splsda.breast)

# default option with no ind name: pch and color are set automatically
plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2))

# default option with no ind name: pch and color are set automatically, with legend
plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2), add.legend = TRUE)

# playing with style
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), plot.indiv = FALSE, plot.ellipse = TRUE, style = "ggplot2", cex = c(1, 1))
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), plot.indiv = FALSE, plot.ellipse = TRUE, style = "lattice", cex = c(1, 1))


## plot of individuals for objects of class 'sgcca' (or 'rgcca')
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
plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE, 
          plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", 
          main = 'my plot')

plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", main = 'my plot')
plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", main = 'my plot')


## plot of individuals for objects of class 'sgccda' 
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
                                     keep = list(c(10,10), c(15,15)),
                                     scheme = "centroid",
                                     verbose = FALSE,
                                     bias = FALSE)

# displaying all blocks. by default colors correspond to outcome Y
plotIndiv(nutrimouse.sgccda1)

# displaying only 2 blocks
plotIndiv(nutrimouse.sgccda1, blocks = c(1,2), group = nutrimouse$diet)

# with some ellipse, legend and title
plotIndiv(nutrimouse.sgccda1, blocks = c(1,2), group = nutrimouse$diet, 
          plot.ellipse = TRUE, add.legend = TRUE, main = 'my sample plot')


