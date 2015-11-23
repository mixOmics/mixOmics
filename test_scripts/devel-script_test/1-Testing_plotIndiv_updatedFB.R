
# -----------------------------------------------------------------------------------
# 0-Testing_plotIndiv.R
# Author:    KA Le Cao 
# Date started:  22/09/2015
# Last updated:  
# Objective: testing update FB dated Sept 2015
# Latest update: 
# -----------------------------------------------------------------------------------

require(lattice)
require(ellipse)
require(ggplot2)

library(mixOmics)
#library(R.utils)

#sourceDirectory(path='C:/Users/Yan/Downloads/mixOmics')

source('plotIndivz.r')  # udpate FB from 19/09/2015


### PLS, sPLS

data(liver.toxicity)

X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3,
                      keepX = c(50, 50, 50), keepY = c(10, 10, 10))


col <- c("blue", "red", "darkgreen", "darkviolet")
pch <- c(15, 16, 17, 18)

cex <- c(1, 1.5, 2, 2.5)



plotIndiv(toxicity.spls, comp = c(1,2), ind.names = rep(c("a","b","c","d"),each=16),
          rep.space = "XY-variate",abline.line=TRUE,group=rep(c(1,2,3,4), each=16))

# test d'erreurs (KA)
plotIndiv(toxicity.spls, comp = c(1,2), ind.names = rep(c("a","b","c","d"),each=4),
          rep.space = "XY-variate",abline.line=TRUE,group=rep(c(1,2,3,4), each=16))
plotIndiv(toxicity.spls, comp = c(1,2), ind.names = rep(c("a","b","c","d"),each=16),
          rep.space = "XY-variate",abline.line=TRUE,group=rep(c(1,2,3,4), each=4))

# voir bug report avec ind.names et pch (prio sur pch)
plotIndiv(toxicity.spls, comp = c(2,1),
          rep.space = "X-variate", X.label="X1",ind.names=FALSE, pch=pch,
          group=rep(c(1,2,3,4), each=16))

# ajout legende: ok
plotIndiv(toxicity.spls,rep.space = "Y-variate",cex=c(1, 1.5, 2, 2.5),
          ind.names=FALSE,style="lattice", plot.ellipse=TRUE,
          group=rep(c(1,2,3,4), each=16), add.legend = TRUE)



plotIndiv(toxicity.spls, comp = c(1,3),style="lattice",
          rep.space="XY-variate")

groupe=rep(c("groupe 1","groupe 2","groupe 3","groupe 4"),each = 16)

#?? weird, pch is by default when group is specified? is this new? 
plotIndiv(toxicity.spls,add.legend=TRUE,abline.line = TRUE,group=groupe)

plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=FALSE,col=col,style="lattice",group=groupe) 

plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=TRUE,style="lattice",group=groupe)
plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=TRUE,style="graphics",group=groupe) #####

plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=FALSE,group=groupe)


plotIndiv(toxicity.spls,ind.names=FALSE,group=groupe,style="graphics")
abline(0,0)
points(5,2)

#adding ylim: not working
plotIndiv(toxicity.spls,ind.names=FALSE,group=groupe,style="graphics", ylim = c(-20, 6))


# show.var: does not work for spls
plotIndiv(toxicity.spls, show.var = TRUE)
# ===================
# rcc
# ==================
data(nutrimouse)

X <- nutrimouse$lipid
Y <- nutrimouse$gene

col=c("red","blue")
pch=c(0,1)
cex=c(1,2)

nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)

plotIndiv(nutri.res,group=rep(c("groupe 1", "groupe 2"),each=20))

plotIndiv(nutri.res,col=col,group=rep(c("groupe 1", "groupe 2"),each=20))

plotIndiv(nutri.res,col=col,group=rep(c("groupe 1", "groupe 2"),each=20),ind.names = F) 

nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 10)
plotIndiv(nutri.res,col=col,style="lattice",add.legend=TRUE,plot.ellipse=TRUE,
          group=rep(c("groupe 1", "groupe 2"),each=20)) ### plot.ellipse

nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 10, lambda2 = 0.008)
plotIndiv(nutri.res,main="plot",Y.label="Y1",comp=c(2,3))


plotIndiv(nutri.res,style="graphics")
abline(0,0)
plotIndiv(nutri.res,style="graphics",add.legend=T)
points(0,0)
plotIndiv(nutri.res,style="graphics")
abline(0,0)


# ---------------
# PLS DA
# ---------------
data(breast.tumors)
X = breast.tumors$gene.exp
Y = breast.tumors$sample$treatment

breast.plsda = plsda(X, Y, ncomp = 2)

plotIndiv(breast.plsda, ind.names = TRUE,style="lattice",add.legend= TRUE)

col=unique(color.mixo(as.numeric(map(breast.plsda$ind.mat))))
col # only 2 colors provided but then arg is recycled. Wrong example?
plotIndiv(breast.plsda,cex=0.5,col=col, ind.names = TRUE,style="graphics", plot.ellipse=TRUE)

plotIndiv(breast.plsda, ind.names = FALSE,ellipse.level=0.5,plot.ellipse=TRUE) 

# ----------------
# sGCCA
# ---------------

diet = unmap(nutrimouse$diet)
blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = c(3, 2, 1))

#blocks
plotIndiv(nutri.sgcca) #### 
plotIndiv(nutri.sgcca, blocks = "lipid", ind.names = TRUE, group = nutrimouse$diet) 

plotIndiv(nutri.sgcca, blocks = "gene", ind.names = nutrimouse$diet, 
          group=nutrimouse$genotype,add.legend=TRUE)


plotIndiv(nutri.sgcca, blocks = "gene", ind.names = FALSE,add.legend=TRUE,
          style="graphics",Y.label="Dim")

plotIndiv(nutri.sgcca, blocks = "gene", ind.names = 1:40,
          style="lattice",col="blue",X.label="Dim",abline.line=TRUE)
plotIndiv(nutri.sgcca, blocks = c("gene","lipid"),
          X.label="Dim",abline.line=TRUE,
          group = nutrimouse$genotype, 
          ind.names = nutrimouse$diet,
          add.legend=TRUE)


# -----------------
# sPCA
# ---------------
# FB: variance cannot be displayed with sPCA, and so it should not appear as NA.
# can you disable for spCA

liver.spca= spca(liver.toxicity$gene, keepX = c(10,10, 10))

plotIndiv(liver.spca, ind.names=FALSE,
          group=liver.toxicity$treatment[, 3],plot.ellipse=TRUE,style="lattice")

plotIndiv(liver.spca, ind.names=FALSE,add.legend=TRUE,plot.ellipse=TRUE,
          group=liver.toxicity$treatment[, 3],style="graphics") ###


plotIndiv(liver.spca, ind.names = FALSE, 
          group=liver.toxicity$treatment[, 3])


# -----------------
# PCA (KA added)
# ---------------
# FB: variance should appear as: 'Expl. var = '
# FB: instead of 'Dimension' what should appear is 'Component' (similar to the other functions)

liver.pca= pca(liver.toxicity$gene)

plotIndiv(liver.pca, ind.names=FALSE,
          group=liver.toxicity$treatment[, 3],plot.ellipse=TRUE,style="lattice")


# -----------------
# IPCA (KA added), needs more tests
# ---------------
# FB: instead of 'Dimension' what should appear is 'Component' (similar to the other functions)

liver.ipca= ipca(liver.toxicity$gene)

plotIndiv(liver.ipca, ind.names=FALSE,
          group=liver.toxicity$treatment[, 3],plot.ellipse=TRUE,style="lattice")



# -----------------
# sIPCA (KA added), needs more tests
# ---------------
# FB: instead of 'Dimension' what should appear is 'Component' (similar to the other functions)


liver.sipca= sipca(liver.toxicity$gene, keepX = c(10,10, 10))

plotIndiv(liver.sipca, ind.names=FALSE,
          group=liver.toxicity$treatment[, 3],plot.ellipse=TRUE,style="lattice")

