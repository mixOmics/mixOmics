# -----------------------------------------------------------------------------------
# Testing-color.R
# Author:    KA Le Cao 
# Date started:  26/02/2015
# Last updated:  12/03/2015
# Objective: testing an update of the new mixOmics palette
# Latest update: 
# -----------------------------------------------------------------------------------

library(mixOmics)

# -----------------------
# existing jet colors
# ----------------------
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
  image(matrix(z, ncol = 1), col = jet.colors(n), 
        xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
  box()
  par(usr = c(-1, 1, -1, 1))  	
  axis(1, at = c(-1, 0, 1))
}

# 
# # --------------------------------
# # including the mixOmics colors
# # -------------------------------
# mixo.gray = gray.colors(1, start = 0.76, gamma = 1)
# 
# mixo.col = c('#388ECC', # mixOmics logo blue
#              '#F68B33', # mixOmics logo orange
#              mixo.gray, # mixOmics logo grey
#              '#009E73', # shiny dark green
#              '#CC79A7', # shiny purple/pink
#              '#F0E442' #shiny yellow
#              #'#D55E00', #shiny dark orange
#              #'#0072B2', #shiny dark blue
#              #'#999999',  # shiny grey
#              #'#E69F00', # shiny orange
#              #'#56B4E9' #Shiny blue
#              )
# 
# par(mfrow=c(1,1))            
# plot( 1:length(mixo.col), pch = 18, col = mixo.col, cex = 5)            


# -------------------
# testing the function mixo.colors
# -------------------
source('color.palettes.R')

library(mixOmics)
help(pca)

# -----------
# multidrug data
# ----------
data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)

my.colors = mixo.colors(as.numeric(factor(multidrug$cell.line$Class)))
# samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 0.5, 
          col = my.colors)

# outputs the 3rd color
mixo.colors(3)

# -------------
# nutrimouse data
# -----------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

my.colors = mixo.colors(as.numeric(nutrimouse$diet))
my.pch = ifelse(nutrimouse$genotype == 'wt', 16, 17)
plotIndiv(nutri.res, ind.names = FALSE, col = my.colors, pch = my.pch, cex = 1.5)

# legend('topleft', c("WT", "PPAR"), pch = c(16, 17), 
#        col = unique(my.colors), text.col = c("blue", "red"),
#        cex = 1, pt.cex = c(1.2, 1.2))


# -------------
# liver.toxicity
# ------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))


my.colors = mixo.colors(as.numeric(as.factor(liver.toxicity$treatment$Time.Group)))
list.pch = c(15, 16, 17, 19)
my.pch = list.pch[as.numeric(as.factor(liver.toxicity$treatment$Dose.Group))]

plotIndiv(toxicity.spls, comp = 1:2, ind.names = FALSE,
          rep.space = "X-variate", col = my.colors, pch = my.pch, cex = 1.5,
          X.label = 'Component 1', Y.label = 'Component 2'
          )


# 3D plot
# -------
list.pch2 = c('s', 't', 'c', 'o')
my.pch2 = list.pch2[as.numeric(as.factor(liver.toxicity$treatment$Dose.Group))]

plot3dIndiv(toxicity.spls, ind.names = FALSE, pch = my.pch2,
            col = my.colors, cex = 15)
library(rgl)
rgl.postscript('example-plotIndiv3dplot.pdf', fmt = 'pdf')


# -------------
# networks
# ------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

pdf('example-network.pdf')
network(toxicity.spls, comp = 1:3, threshold = 0.6, 
        X.names = NULL, Y.names = NULL, keep.var = TRUE,
        color.node = mixo.colors(c(1, 2)),
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

pdf('example-plotVar.pdf')
plotVar(toxicity.spls, keep.var = TRUE, Y.label = TRUE, cex = c(1,0.6), col = mixo.colors(c(1,2)))  
dev.off()


plot3dVar(toxicity.spls, var.label = FALSE, Y.label = TRUE, keep.var = TRUE, 
          rad.in = 0.5, cex = c(1,0.6),
          col = mixo.colors(c(1,2)))
library(rgl)
rgl.postscript('example-plotVar3dplot.pdf', fmt = 'pdf')




