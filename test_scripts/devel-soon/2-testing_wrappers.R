# -----------------------------------------------------------------------------------
# Testing-wrappers.R
# Author:    KA Le Cao 
# Date started:  21/07/2015
# Last updated:  
# Objective: testing for the multi block module with wrappers (follows the example.R file)
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


# ==========================================
# RGCCA example with nutrimouse: only on 2 data sets!
# ==========================================
data(nutrimouse)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)

nutrimouse.rgcca1 <- wrapper.rgcca(blocks = data,
                                          design = design1,
                                          tau = "optimal",
                                          ncomp = c(2, 2),
                                          scheme = "centroid",
                                          verbose = FALSE)
# have a look at the loading vectors
head(nutrimouse.rgcca1$loadings[[1]])
head(nutrimouse.rgcca1$loadings[[2]])


# plotIndiv
# ----------
plotIndiv(nutrimouse.rgcca1, blocks = NULL)

# yups, all good
plotIndiv(nutrimouse.rgcca1, blocks = c(1,2), group = nutrimouse$diet)

# with some ellipse
plotIndiv(nutrimouse.rgcca1, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)
# and legend
plotIndiv(nutrimouse.rgcca1, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE, add.legend = TRUE, main = 'my sample plot')

# this should throw an error, as in rgcca we project the samples into their own space
plotIndiv(nutrimouse.rgcca1, group = nutrimouse$diet, rep.space = "XY-variate")
# like this:
plotIndiv(nutrimouse.rgcca1, group = nutrimouse$diet, style = "lattice", rep.space = "X-variate")
plotIndiv(nutrimouse.rgcca1, group = nutrimouse$diet, style = "graphics", rep.space = "X-variate")






# plotVar: needs to update the .Rd file (KA)
# ------
# both data sets at once
plotVar(nutrimouse.rgcca1, block = c(1,2), col = c('green', 'blue'), cex = c(1,1))
plotVar(nutrimouse.rgcca1, block = c(1,2), col = color.mixo(c(1,2)), cex = c(1,1))

# one data set: wil throw out a warning message as it expects at least 2 data sets.
plotVar(nutrimouse.rgcca1, block = c(2))
# one data set
plotVar(nutrimouse.rgcca1, block = c(1))



# ==========================================
# sGCCA example with nutrimouse
# ==========================================
rm(list=ls())
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R

data(nutrimouse)
Y = unmap(nutrimouse$diet) # here need to unmap as we are not using a supervised method
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# 1 - using the function with a penalty
# ---------------------------------------
nutrimouse.sgcca1 <- wrapper.sgcca(blocks = data,
                                          design = design1,
                                          penalty = c(0.3, 0.5, 1),
                                          ncomp = c(2, 2, 1),
                                          scheme = "centroid",
                                          verbose = FALSE, 
                                          bias = FALSE)

nutrimouse.sgcca1$loadings[[1]]
nutrimouse.sgcca1$loadings[[2]]
nutrimouse.sgcca1$loadings[[3]]


# plotIndiv
# ----------
# will throw an error as Y is deflated only once and per default ncomp = 2 on each block
plotIndiv(nutrimouse.sgcca1, blocks = NULL)

# yups, all good
plotIndiv(nutrimouse.sgcca1, blocks = c(1,2), group = nutrimouse$diet)

# with some ellipse
plotIndiv(nutrimouse.sgcca1, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)
# and legend
plotIndiv(nutrimouse.sgcca1, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE, add.legend = TRUE, main = 'my sample plot')

# plotVar: needs to update the .Rd file (KA)
# ------
# both data sets at once
plotVar(nutrimouse.sgcca1, block = c(1,2), col = color.mixo(c(1,2)), cex = c(1,1))

# one data set: wil throw out a warning message as it expects at least 2 data sets.
plotVar(nutrimouse.sgcca1, block = c(2), cex = 1)
# one data set
plotVar(nutrimouse.sgcca1, block = c(1), cex = 1)


# variables selected
# --------
# variables selected from block1 on comp 1
selectVar(nutrimouse.sgcca1, block = 1, comp = 1)

# variables selected from block2 on comp 1
selectVar(nutrimouse.sgcca1, block = 2, comp = 1)


# 2 - as an alternative: using the function with keep.block
# ---------------------------------------
data(nutrimouse)
Y = unmap(nutrimouse$diet) # here need to unmap as we are not using a supervised method
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


# using keep.blocks: should be a list for each block, and the number of component per block
list.keep = list(keep.gene = c(30,30), keep.lipid = c(10,10), keepY = ncol(Y))  # for the outcome we need to keep all variables, i.e. ncol(unmap(outcome))

nutrimouse.sgcca.keep <- wrapper.sgcca(blocks = data,
                                   design = design1,
                                   keep.blocks = list.keep,
                                   ncomp = c(2, 2, 1),
                                   scheme = "centroid",
                                   verbose = FALSE, 
                                   bias = FALSE)

nutrimouse.sgcca.keep$loadings[[1]]
nutrimouse.sgcca.keep$loadings[[2]]
nutrimouse.sgcca.keep$loadings[[3]]


# plotIndiv
# ----------
# will throw an error as Y is deflated only once and per default ncomp = 2 on each block
plotIndiv(nutrimouse.sgcca.keep, blocks = NULL)

# yups, all good
plotIndiv(nutrimouse.sgcca.keep, blocks = c(1,2), group = nutrimouse$diet)

# with some ellipse
plotIndiv(nutrimouse.sgcca.keep, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)
# and legend
plotIndiv(nutrimouse.sgcca.keep, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE, add.legend = TRUE, main = 'my sample plot')

# plotVar: needs to update the .Rd file (KA)
# ------
# both data sets at once
plotVar(nutrimouse.sgcca.keep, block = c(1,2), col = color.mixo(c(1,2)), cex = c(1,1))

# one data set: wil throw out a warning message as it expects at least 2 data sets.
plotVar(nutrimouse.sgcca.keep, block = c(2), cex = 1)
# one data set
plotVar(nutrimouse.sgcca.keep, block = c(1), cex = 1)


# variables selected
# --------
# variables selected from block1 on comp 1
selectVar(nutrimouse.sgcca.keep, block = 1, comp = 1)

# variables selected from block2 on comp 1
selectVar(nutrimouse.sgcca.keep, block = 2, comp = 1)




# ==========================================
# sGCCA-DA example with nutrimouse
# ==========================================
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutrimouse.sgccda1 <- wrapper.sgccda(blocks = data,
                                     Y = Y,
                                     design = design1,
                                     ncomp = c(2, 2),
                                     keep.blocks = list(c(10,10), c(15,15)),
                                     scheme = "centroid",
                                     verbose = FALSE,
                                     bias = FALSE)
# FR: errors:need to put Y as an unmapped matrix
# see my changes in wrappers.R, l 59-63 but did not work
Design matrix has changed to include Y
Error: 'block 3' must be a numeric matrix.

## I stoppped here, but maybe selectVar does not work

# plotIndiv
# ----------
# should throw an error as we request all blocks to be displayed and Y has only 1 comp
plotIndiv(nutrimouse.rgcca1, blocks = NULL)

# yups, all good
plotIndiv(nutrimouse.rgcca1, blocks = c(1,2), group = nutrimouse$diet)

# with some ellipse
plotIndiv(nutrimouse.rgcca1, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)
# and legend
plotIndiv(nutrimouse.rgcca1, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE, add.legend = TRUE, main = 'my sample plot')

head(nutrimouse.rgcca1$loadings[[1]])
head(nutrimouse.rgcca1$loadings[[2]])
head(nutrimouse.rgcca1$loadings[[3]])


# plotVar: needs to update the .Rd file (KA)
# ------
# both data sets at once
plotVar(nutrimouse.rgcca1, block = c(1,2), col = c('green', 'blue'), cex = c(1,1))
plotVar(nutrimouse.rgcca1, block = c(1,2), col = color.mixo(c(1,2)), cex = c(1,1))

# one data set
plotVar(nutrimouse.rgcca1, block = c(2))
# one data set
plotVar(nutrimouse.rgcca1, block = c(1))



# variables selected?
# --------
# FR: not working
selectVar(nutrimouse.rgcca1, block = 1, comp = 1)
# Error in UseMethod("selectVar") : 
#   no applicable method for 'selectVar' applied to an object of class "rgcca"


