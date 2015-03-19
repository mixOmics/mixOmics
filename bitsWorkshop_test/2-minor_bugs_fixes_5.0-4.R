# ----------------------------------------
# bug fixes for update 5.0-4
# date: 18/02/2015
# latest update: 
# ----------------------------------------

setwd("~/Documents/k.lecao/Packages/mixOmics/GIT/package-mixomics")

# ------ notes for me to compile the package (if need be) on a terminal
R CMD build --resave-data mixOmics
R CMD INSTALL -l MyR/ mixOmics_5.0-4.tar.gz
R CMD check mixOmics --as-cran --timings
# ---------------------------------


# now in R, load the package
detach("package:mixOmics", unload=TRUE)

library(mixOmics, lib.loc = 'MyR/')
# check that version 5.0-4 is loaded
sessionInfo() # ok mixOmics 5.0-4

# ============================================================
# 1 - small fix: network
# =============================================================
# 1- network
# network function: v <-
#   network default: red and green

# in the network function, I have only changed the default color to color.GreenRed (instead of blue, red as
# previously set)
source('mixOmics/R/network.R')
library(igraph)  # for some reason needed to upload that to use graph.data.frame

# testing the help file
help(network)

## network representation for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
jpeg('graphics/example1-network.jpeg')
network(nutri.res, comp = 1:3, threshold = 0.6)
dev.off()

# ?! bug to fix?: when setting interactive = TRUE:
# Error in network.default(simMat, ...) : attempt to apply non-function# or: error is: 
# Error in par(def.par) : 
# invalid value specified for graphical parameter "pin"



jpeg('graphics/example2-network.jpeg')
network(nutri.res, comp = 1:3, threshold = 0.45,
        color.node = c("mistyrose", "lightcyan"),
        shape.node = c("circle", "rectangle"), 
        color.edge = color.jet(100),
        lty.edge = c("solid", "solid"), lwd.edge = c(2, 2), 
        show.edge.labels = FALSE)
dev.off()


## network representation for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
## Not run: 
jpeg('graphics/example3-network.jpeg')
network(toxicity.spls, comp = 1:3, threshold = 0.8, 
        X.names = NULL, Y.names = NULL, keep.var = TRUE,
        color.node = c("mistyrose", "lightcyan"),
        shape.node = c("rectangle", "circle"),
        color.edge = color.spectral(100),
        lty.edge = c("solid", "solid"), lwd.edge = c(1, 1), 
        show.edge.labels = FALSE, interactive = FALSE)
dev.off()
## End(Not run)

# ============================================================
# 2 - small fix: splsda.Rd: email from Carlos
# 
# Carlos email
# . I was using the example with the breast.tumors dataset available at http://127.0.0.1:30821/library/mixOmics/html/plsda.html
# 
# and I think there's an error in the legend plotting, since the points that belong to the before treatment appear as after treatment and viceversa.
# 
# I guess it's because when the color palette is specified with
# palette(c("red", "blue"))
# 
# the "After" points are assigned red and the "Before" points are colored with blue. And in the legend the color blue shows first so it is labeled as "After".
# 
# For example the point 23 is colored blue and labeled as After treatment in the legend, but in the raw data it's assigned as before treatment.

# =============================================================
## First example
data(breast.tumors)
X <- breast.tumors$gene.exp
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(breast.tumors$sample$treatment)

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))

col.breast <- color.mixo(as.numeric(Y))

# individual names appear
plotIndiv(res, ind.names = Y, col = col.breast)
legend('topright', c('Before', 'After'), pch = c(16, 16), 
       col = unique(col.breast), 
       title = "Treatment")

# check:color match
data.frame(col.breast, Y)


## Second example
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))

col.rat <- color.mixo(Y)
# individual name is set to the treatment
plotIndiv(splsda.liver, col = col.rat, ind.names = Y)



# ============================================================
# 3 - small fix: sgcca, email from Wenbo (?)
# I think there is a bug in the plotVar function because it is plotting variables not selected for block 2. The code I used to plot selected variables in block 2 is 
# 
# plotVar(wrap.result.sgcca, comp = c(1,2), block = 2, ncomp.select = 1, labels = TRUE, cex=0.5, col = 'blue').
# 
# If I look at the loadings for that block and comp, some of those variables plotted have zero loading values.
# =============================================================

help(wrapper.sgcca)

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
wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.3, 1), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
wrap.result.sgcca
#did the algo converge?
wrap.result.sgcca$crit  # yes



help(plotVar)
## variable representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------

# send to wenbo with updated code:
source('mixOmics/R/select.var.R')
source('mixOmics/R/plotVar.R')
# also updated the color.mixo in plotVar, so would need source('mixOmics/R/color.mixo.R')

#variables selected on component 1 for the two blocs
select.var(wrap.result.sgcca, comp = 1, block = c(1,2))$name.var

#variables selected on component 2 for each block 
select.var(wrap.result.sgcca, comp = 2, block = c(1,2))$name.var


plotVar(wrap.result.sgcca, comp = c(1,2), block = c(1,2), ncomp.select = c(1,1), labels = TRUE)
title(main = c('Variables selected on component 1 only'))
plotVar(wrap.result.sgcca, comp = c(1,2), block = c(1,2), ncomp.select = c(2,2), labels = TRUE)
title(main = c('Variables selected on component 2 only'))
# -> this one shows the variables selected on both components
plotVar(wrap.result.sgcca, comp = c(1,2), block = c(1,2), labels = TRUE)
title(main = c('Variables selected on components 1 and 2'))


# end send to Wenbo

## variable representation for objects of class 'rgcca'
# ----------------------------------------------------
data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
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
wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
#wrap.result.rgcca
plotVar(wrap.result.rgcca, comp = c(1,2), block = c(1,2))
plotVar(wrap.result.rgcca, comp = c(1,2), block = c(1,2), labels = TRUE)


