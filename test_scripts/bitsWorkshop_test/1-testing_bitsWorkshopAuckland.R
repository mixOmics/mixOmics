# ----------------------------------------
# testing few bits for the Auckland workshop
# date: 17/02/2015
# latest update: fixing s.match
# ----------------------------------------

setwd("~/Documents/k.lecao/Packages/mixOmics/GIT/package-mixomics")

# ------ notes for me to compile the package (if need be) on a terminal
R CMD build --resave-data mixOmics
R CMD check mixOmics_5.0-4.tar.gz --as-cran --timings
R CMD INSTALL -l MyR/ mixOmics_5.0-4.tar.gz
# ---------------------------------


# now in R, load the package
detach("package:mixOmics", unload=TRUE)

library(mixOmics, lib.loc = 'MyR/')
# check that version 5.0-4 is loaded
sessionInfo() # ok mixOmics 5.0-4

# ---------------------
#

source('mixOmics/R/s.match.R')


# =======================
#    s.match
# =======================
# during the workshop in Toulouse, it appeared that the s.match did not adjust for a change of signs.
# below are the script from the workshop
# look especially at example with rgcca to test several scenarios

# =============================== simu in graphical outputs =================

###################################################
### code chunk number 48: graph_output.Rnw:240-242
###################################################

source('bitsWorkshop_test/Simulation_data.R')

simu.res = pls(X, Y, ncomp = 3, mode = "canonical")


###################################################
### code chunk number 49: graph_output.Rnw:255-269
###################################################
ind <- 7
col.simu <- rep('black', nrow(X))
col.simu[ind] <- 'red'

par(mfrow = c(1, 3), cex = 1.2)
plotIndiv(simu.res, rep.space = 'X-variate', col = col.simu)
title(main = 'Simulated data, PLS X-subspace')

plotIndiv(simu.res, rep.space = 'Y-variate', col = col.simu)
title(main = 'Simulated data, PLS Y-subspace')

plotIndiv(simu.res, rep.space = 'XY-variate', col = col.simu)
title(main = 'Simulated data, PLS XY-subspace')
par(mfrow = c(1, 1))


###################################################
### code chunk number 50: graph_output.Rnw:303-309
###################################################
s.match(simu.res$variates$X[, c(1, 2)], simu.res$variates$Y[, c(1, 2)],
        col = col.simu)

abline(
  v = c(simu.res$variates$X[ind, 1], simu.res$variates$Y[ind, 1]),
  h = c(simu.res$variates$X[ind, 2], simu.res$variates$Y[ind, 2]),
  col = 3:4, lty = 2)


# ================================ nutrimouse rcc ===========

data(nutrimouse)
X <- nutrimouse$lipid
dim(X)
Y <- nutrimouse$gene
dim(Y)
#just a head:
head(cbind(rownames(X), rownames(Y)))

load('bitsWorkshop_test/save-cv.score.RData')

nutrimouse.rcc <- rcc(X,Y,ncomp = 3, lambda1 = cv.score$opt.lambda1, 
                      lambda2 = cv.score$opt.lambda2)

col.nutri <- as.numeric(nutrimouse$genotype) + 2

###################################################
### code chunk number 79: NT3
###################################################
s.match(nutrimouse.rcc$variates$X[,1:2], 
        nutrimouse.rcc$variates$Y[,1:2], col = col.nutri)
title('Nutrimouse, arrow plot')



# ===============================rGCCA ========================
data(nutrimouse)
data = list(nutrimouse$gene, nutrimouse$lipid, nutrimouse$gene)  # !! change workshop here

design1 = matrix(c(0,1,1,
                   1,0,1,
                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
design1


nutrimouse.rgcca1 = wrapper.rgcca(data = data, design = design1, 
                                  tau = 'optimal', 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
nutrimouse.rgcca1

col.diet <- as.numeric(nutrimouse$diet) 
pch.geno <- c(rep(16, 20), rep(17, 20))



###################################################
### code chunk number 165: rgcca1-s.match
###################################################

# this one should be ok
s.match(nutrimouse.rgcca1$variates[[1]], nutrimouse.rgcca1$variates[[2]], col = col.diet)

# this is when we simulate one component having a negative sign (here I have tested several combination)
# by changing the - sign
variates.gene.D1 = cbind(-nutrimouse.rgcca1$variates[[1]][,1],nutrimouse.rgcca1$variates[[1]][,2])
variates.lipid.D1 = cbind(-nutrimouse.rgcca1$variates[[2]][,1],-nutrimouse.rgcca1$variates[[2]][,2])

cor(variates.gene.D1, variates.lipid.D1)

plot(variates.gene.D1, col = col.diet, pch = pch.geno)
plot(variates.lipid.D1, col = col.diet, pch = pch.geno)
s.match(variates.gene.D1, variates.lipid.D1, col = col.diet)



