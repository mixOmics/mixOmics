# -----------------------------------------------------------------------------------
# Testing-wrapper.R
# Author:    B Gautier, KA Le Cao 
# Date started:  28/07/2015
# Last updated:  
# Objective:
# Latest update: 
# -----------------------------------------------------------------------------------

library(mixOmics)
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

# KA
sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R


### RGCCA # 
# -----------

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

# blocks should specify the block data set where the sample plot can be performed (ideally when there are >= 2 components!)
# we indicate the diet variable colors.
plotIndiv(nutrimouse.rgcca, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)

# have a look at the looadings
head(nutrimouse.rgcca$loadings[[1]])
head(nutrimouse.rgcca$loadings[[2]])
head(nutrimouse.rgcca$loadings[[3]])


## sGCCA
# -------------
# same data as above but sparse approach

# version 1 using the penalisation penalty criterion
# ---
nutrimouse.sgcca <- wrapper.sgcca(blocks = data,
                                   design = design,
                                   penalty = c(0.3, 0.5, 1),
                                   ncomp = c(2, 2, 1),
                                   scheme = "centroid",
                                   verbose = FALSE, 
                                   bias = FALSE)

# blocks should specify the block data set where the sample plot can be performed (ideally when there are >= 2 components!)
# we indicate the diet variable colors.
plotIndiv(nutrimouse.sgcca, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgcca, comp = 1, block = 1)
selectVar(nutrimouse.sgcca, comp = 1, block = 2)

# variable plot on the selected variables
plotVar(nutrimouse.sgcca, col = color.mixo(1:2), cex = c(2,2))

# version 2 using the keep.block penalty criterion (number of variables to keep)
# it is a list per block and per component, need to specify all variables for the Y 'outcome' here 
# (see below for sgccda code, which is more appropriate)
# ----
nutrimouse.sgcca <- wrapper.sgcca(blocks = data,
                                  design = design,
                                  ncomp = c(2, 2, 1),
                                  # for keep: each element of the list corresponds to a block 
                                  # and is of length the # comp per block
                                  keep = list(c(10,10), c(15,15), c(ncol(Y))),
                                  scheme = "centroid",
                                  verbose = FALSE, 
                                  bias = FALSE)

# blocks should specify the block data set where the sample plot can be performed (ideally when there are >= 2 components!)
# we indicate the diet variable colors.
plotIndiv(nutrimouse.sgcca, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgcca, comp = 1, block = 1)
selectVar(nutrimouse.sgcca, comp = 1, block = 2)


## sGCC-DA
# -------------
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutrimouse.sgccda <- wrapper.sgccda(blocks = data,
                                    Y = Y,
                                    design = design,
                                    keep = list(c(10,10), c(15,15)),
                                    ncomp = c(2, 2, 1),
                                    scheme = "centroid",
                                    verbose = FALSE,
                                    bias = FALSE)


# blocks should specify the block data set where the sample plot can be performed (ideally when there are >= 2 components!)
# we indicate the diet variable colors.
plotIndiv(nutrimouse.sgccda, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgccda, comp = 1, block = 1)
selectVar(nutrimouse.sgccda, comp = 1, block = 2)

# variable plot on the selected variables
plotVar(nutrimouse.sgccda, col = color.mixo(1:2), cex = c(2,2))

