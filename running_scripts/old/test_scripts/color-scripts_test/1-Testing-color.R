# -----------------------------------------------------------------------------------
# Testing-color.R
# Author:    KA Le Cao 
# Date started:  26/02/2015
# Last updated:  16/03/2015
# Objective: testing an update of the new mixOmics palette
# Latest update: breaking up all functions into different files, update of the .Rd
# -----------------------------------------------------------------------------------

library(mixOmics)

source('../mixOmics/R/color.jet.R')
# -----------------------
# existing jet colors
# ----------------------
jpeg('../graphics/color.jet.jpeg')
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
  image(matrix(z, ncol = 1), col = color.jet(n), 
        xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
  box()
  par(usr = c(-1, 1, -1, 1))  	
  axis(1, at = c(-1, 0, 1))
}
dev.off()

source('../mixOmics/R/color.spectral.R')
# -----------------------
#new: spectral colors
# ----------------------
jpeg('../graphics/color.spectral.jpeg')
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
  image(matrix(z, ncol = 1), col = color.spectral(n), 
        xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
  box()
  par(usr = c(-1, 1, -1, 1))    
  axis(1, at = c(-1, 0, 1))
}
dev.off()

source('../mixOmics/R/color.GreenRed.R')
# -----------------------
#new: GreenRed colors
# ----------------------
jpeg('../graphics/color.GreenRed.jpeg')
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
  image(matrix(z, ncol = 1), col = color.GreenRed(n), 
        xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
  box()
  par(usr = c(-1, 1, -1, 1))    
  axis(1, at = c(-1, 0, 1))
}
dev.off()


# # --------------------------------
# new: including the mixOmics colors
# # -------------------------------

# -------------------
# testing the function mixo.colors
# -------------------
source('../mixOmics/R/color.mixo.R')

library(mixOmics)
help(pca)

# -----------
# multidrug data
# ----------
data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)

my.colors = color.mixo(as.numeric(factor(multidrug$cell.line$Class)))
# samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 0.5, 
          col = my.colors)

# outputs the 3rd color
color.mixo(3)

# -------------
# nutrimouse data
# -----------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

my.colors = color.mixo(as.numeric(nutrimouse$diet))
my.pch = ifelse(nutrimouse$genotype == 'wt', 16, 17)
plotIndiv(nutri.res, ind.names = FALSE, col = my.colors, pch = my.pch, cex = 1.5)

# -------------
# liver.toxicity
# ------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))


my.colors = color.mixo(as.numeric(as.factor(liver.toxicity$treatment$Time.Group)))
list.pch = c(15, 16, 17, 19)
my.pch = list.pch[as.numeric(as.factor(liver.toxicity$treatment$Dose.Group))]


pdf('../graphics/color.plotIndiv.pdf')
plotIndiv(toxicity.spls, comp = 1:2, ind.names = FALSE,
          rep.space = "X-variate", col = my.colors, pch = my.pch, cex = 1.5,
          X.label = 'Component 1', Y.label = 'Component 2'
          )
dev.off()

# 3D plot
# -------
list.pch2 = c('s', 't', 'c', 'o')
my.pch2 = list.pch2[as.numeric(as.factor(liver.toxicity$treatment$Dose.Group))]

plot3dIndiv(toxicity.spls, ind.names = FALSE, pch = my.pch2,
            col = my.colors, cex = 15)
library(rgl)
rgl.postscript('../graphics/color-plotIndiv3dplot.pdf', fmt = 'pdf')


# -------------
# networks
# ------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

pdf('../graphics/color-network.pdf')
network(toxicity.spls, comp = 1:3, threshold = 0.6, 
        X.names = NULL, Y.names = NULL, keep.var = TRUE,
        color.node = color.mixo(c(1, 2)),
        shape.node = c("rectangle", "circle"),
        color.edge =jet.colors(100),
        lty.edge = c("solid", "solid"), lwd.edge = c(1, 1), 
        show.edge.labels = FALSE, interactive = FALSE)
dev.off()

# -------------
# plotVar
# ------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(20, 20, 20), 
                      keepY = c(10, 10, 10))

pdf('../graphics/color-plotVar.pdf')
plotVar(toxicity.spls, keep.var = TRUE, Y.label = TRUE, cex = c(1,0.6), col = color.mixo(c(1,2)))  
dev.off()


plot3dVar(toxicity.spls, var.label = FALSE, Y.label = TRUE, keep.var = TRUE, 
          rad.in = 0.5, cex = c(1,0.6),
          col = color.mixo(c(1,2)))
library(rgl)
rgl.postscript('../graphics/color-plotVar3dplot.pdf', fmt = 'pdf')




