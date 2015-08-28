# -----------------------------------------------------------------------------------
# Testing_plotVar.R
# Author:    KA Le Cao 
# Date started:  28/08/2015
# Last updated:  
# Objective: 
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
##sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R


library(mixOmics, lib.loc = '../../MyR/')
# KA
#sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R
#source('../../../plotVar.R')

#library(mixOmics)

# =====================================
# testing examples first
# =====================================

## variable representation for objects of class 'rcc'
# ----------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

plotVar(nutri.res) #(default)

# playing with the style
plotVar(nutri.res, style = 'lattice') #(default)

# changing x and y labels
plotVar(nutri.res, comp = c(1,3), cutoff = 0.5, 
        X.label = 'PC1', Y.label = 'PC3')

# one correlation circle plot per data set
plotVar(nutri.res, comp = c(1,2), cutoff = 0.5, 
        overlap = FALSE)


# with pch symbols
plotVar(nutri.res, comp = c(1,2), pch = c(16,2))



## variable representation for objects of class 'pls' or 'spls'
# ----------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

# default shows only the variables selected on the plotted components
plotVar(toxicity.spls)

# shows only the variables selected on the plotted components
plotVar(toxicity.spls, comp = c(1,3))

# shows only the variables selected on the selected components
plotVar(toxicity.spls, comp.select = c(1:3))


# change variable names
new.names = list(paste('gene', 1:ncol(X)), paste('clinic', 1:ncol(Y)))
plotVar(toxicity.spls, overlap = FALSE, var.names = new.names)

# prefilter even further and use of pch
plotVar(toxicity.spls, comp.select = c(1:3), cutoff = 0.8, pch = c(15,16))

# change colors
plotVar(toxicity.spls, col = color.mixo(3:4))

my.col = list(c(rep(1, ncol(X))), c(rep(3,ncol(Y))))
plotVar(toxicity.spls, col = my.col)


## variable representation for objects of class 'splsda'
# ----------------------------------------------------
\dontrun{
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- as.factor(liver.toxicity$treatment[, 4])
  
  ncomp <- 2
  keepX <- rep(20, ncomp)
  
  splsda.liver <- splsda(X, Y, ncomp = ncomp, keepX = keepX)
  # use of pch symbols
  plotVar(splsda.liver, pch = 16, col = 3)
}



## variable representation for objects of class 'sgcca' 
# ----------------------------------------------------
## see example in ??wrapper.sgcca
data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
# set up the data as list
data = list(nutrimouse$gene, nutrimouse$lipid,Y)

# set up the design matrix:
# with this design, gene expression and lipids are connected to the diet factor
# design = matrix(c(0,0,1,
#                   0,0,1,
#                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


#note: the penalty parameters will need to be tuned
wrap.result.sgcca = wrapper.sgcca(blocks = data, design = design, penalty = c(.3,.3, 1),
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
#wrap.result.sgcca


# showing 2 blocks, with variables selected on comp 1 for block 1 and comp 1 for block 2
plotVar(wrap.result.sgcca, comp = c(1,2), 
        blocks = c(1,2), comp.select = c(1,1), 
        overlap = FALSE,
        main = 'Variables selected on component 1 only')


# displaying variables selected on comp 2 for block 1 and comp 2 for block 2
plotVar(wrap.result.sgcca, comp = c(1,2), blocks = c(1,2), comp.select = c(2,2), 
        main = 'Variables selected on component 2 only')


## variable representation for objects of class 'rgcca'
# ----------------------------------------------------
data(nutrimouse)
# need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
Y = unmap(nutrimouse$diet)

data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, 
                byrow = TRUE, dimnames = list(names(data), names(data)))

nutrimouse.rgcca <- wrapper.rgcca(blocks = data,
                                  design = design,
                                  tau = "optimal",
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid",
                                  verbose = FALSE)

# changing cex
plotVar(nutrimouse.rgcca, comp = c(1,2), blocks = c(1,2), cex = c(1.5, 1.5))
# changing font
plotVar(nutrimouse.rgcca, comp = c(1,2), blocks = c(1,2), font = c(1,3))


# set up the data as list
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
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
wrap.result.rgcca = wrapper.rgcca(blocks = data, design = design, tau = c(1, 1, 0),
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
#wrap.result.rgcca
plotVar(wrap.result.rgcca, comp = c(1,2), blocks = c(1,2))
}

