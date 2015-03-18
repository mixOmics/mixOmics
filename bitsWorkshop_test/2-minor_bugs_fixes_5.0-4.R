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
# 1 - small fix: splsda.Rd: email from Carlos
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


