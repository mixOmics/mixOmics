# -----------------------------------------------------------------------------------
# Bugreport_plotIndiv.R
# Author:    KA Le Cao 
# Date started:  22/09/2015
# Last updated:  22/09/2015
# Objective: testing plotIndiv and reporting bugs
# Latest update: 
# -----------------------------------------------------------------------------------

# KA
#sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R
library(mixOmics)
source('plotIndivz.r')  # udpate FB from 19/09/2015
library(ellipse)

# --------------------------------------
# 1 - adding ylim and xlim in a graphics plot?
# --------------------------------------
data(linnerud) 
X <- linnerud$exercise 
Y <- linnerud$physiological 
linn.pls <- pls(X, Y, ncomp = 2, mode = "regression")
indiv1 <- c(200, 40, 60) 
indiv2 <- c(190, 45, 45) 
newdata <- rbind(indiv1, indiv2) 
colnames(newdata) <- colnames(X) 
newdata
pred <- predict(linn.pls, newdata)
pred$variates  # here, pb as prediction will be outside the plot for indiv 2
plotIndiv(linn.pls, comp = 1:2, rep.space = "Y-variate", style = 'graphics') 

# and ylim does not work
plotIndiv(linn.pls, comp = 1:2, rep.space = "Y-variate", style = 'graphics', ylim = c(-3,3)) 

points(pred$variates[, 1], pred$variates[, 2], pch = 19, cex = 1.2) 
text(pred$variates[, 1], pred$variates[, 2], c("new ind.1", "new ind.2"), pos = 3)


# -------------------------
# 2 - when pch specified, need to put ind.names = FALSE internally
#-----------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# or we can specify the argument group for plotting the ellipse according to group
# legend will indicate the groups specified in argument 'group'
plotIndiv(nutri.res, col= as.numeric(nutrimouse$diet), 
          plot.ellipse = TRUE, group = nutrimouse$genotype, 
          add.legend = TRUE, 
          #ind.names = FALSE, # i.e: no need to add this arg when pch is filled
          pch = 19)


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


# for the block 'lipid' with ellipse plots and legend and pch, different styles
# see bug report
plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE, 
          plot.ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", 
          ind.name = FALSE,
          pch = c(15:19),
          main = 'my plot')

# ----------------
# 3 - legend for graphics style? (lattice is ok)
# --------------
# classic style, in the Y space, is the legend far off here when you zoom in? see bug report
plotIndiv(nutri.res, rep.space= 'Y-variate', group = nutrimouse$genotype, 
          add.legend = TRUE, style = 'graphics')



