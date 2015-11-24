require(lattice)
require(ellipse)
require(ggplot2)

library(mixOmics)


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


plotIndiv(toxicity.spls, comp = c(2,1),
          rep.space = "X-variate", X.label="X1",ind.names=FALSE,pch=pch,
          group=rep(c(1,2,3,4), each=16))


plotIndiv(toxicity.spls,rep.space = "Y-variate",cex=c(1, 1.5, 2, 2.5),
          ind.names=FALSE,style="lattice", plot.ellipse=TRUE,
          group=rep(c(1,2,3,4), each=16))



plotIndiv(toxicity.spls, comp = c(1,3),style="lattice",
          rep.space="XY-variate")

groupe=rep(c("groupe 1","groupe 2","groupe 3","groupe 4"),each = 16)

plotIndiv(toxicity.spls,add.legend=TRUE,abline.line = TRUE,group=groupe)

plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=FALSE,col=col,style="lattice",group=groupe) 

plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=TRUE,style="lattice",group=groupe)
plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=TRUE,style="graphics",group=groupe) #####

plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=FALSE,group=groupe)


plotIndiv(toxicity.spls,ind.names=FALSE,group=groupe,style="graphics")
abline(0,0)
points(5,2)



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

data(breast.tumors)
X = breast.tumors$gene.exp
Y = breast.tumors$sample$treatment

breast.plsda = plsda(X, Y, ncomp = 3)
plotIndiv(breast.plsda, ind.names = TRUE,style="lattice",add.legend=FALSE)

col=unique(color.mixo(as.numeric(map(breast.plsda$ind.mat))))


plotIndiv(breast.plsda,cex=2,col=col, ind.names = TRUE,style="graphics",plot.ellipse=TRUE)

plotIndiv(breast.plsda, ind.names = FALSE,ellipse.level=0.5,plot.ellipse=TRUE) 



diet = unmap(nutrimouse$diet)
blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutri.sgcca <- wrapper.sgcca(blocks, ncomp = c(3, 2, 1))####
nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = c(3, 3, 3))

#blocks
plotIndiv(nutri.sgcca) #### 
plotIndiv(nutri.sgcca, blocks = "lipid", ind.names = TRUE) 

plotIndiv(nutri.sgcca, blocks = "gene", ind.names = nutrimouse$diet, 
          group=nutrimouse$diet,add.legend=TRUE)


plotIndiv(nutri.sgcca, blocks = "gene", ind.names = FALSE,add.legend=TRUE,
          style="graphics",Y.label="Dim")
plotIndiv(nutri.sgcca, blocks = "gene", ind.names = 1:40,
          style="lattice",col="blue",X.label="Dim",abline.line=TRUE)
plotIndiv(nutri.sgcca, blocks = c("gene","lipid"),
          col="blue",X.label="Dim",abline.line=TRUE,add.legend=TRUE)



liver.spca= spca(liver.toxicity$gene,ncomp=3)

plotIndiv(liver.spca, ind.names=FALSE,
          group=liver.toxicity$treatment[, 3],plot.ellipse=TRUE,style="lattice")
plotIndiv(liver.spca, ind.names=FALSE,add.legend=TRUE,plot.ellipse=TRUE,
          group=liver.toxicity$treatment[, 3],style="graphics") ###


plotIndiv(liver.spca, ind.names = FALSE, 
          group=liver.toxicity$treatment[, 3])

data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
plotIndiv(pca.res)
plotIndiv(pca.res,X.label="X")
plotIndiv(pca.res,Y.label="Y")
plotIndiv(pca.res,Z.label="Z",style="3d")

plotIndiv(toxicity.spls,ylim=list(c(1,5),c(1,5)),style="graphics")
plotIndiv(toxicity.spls,ylim=list(c(1,5),c(1,5)),style="ggplot2")
plotIndiv(toxicity.spls,xlim=list(c(1,5),c(1,5)),style="lattice")
plotIndiv(nutri.sgcca,xlim=list(c(1,5),c(4,5),c(3,6)),style="graphics")
plotIndiv(nutri.sgcca,ylim=list(c(1,5),c(4,5),c(3,6)),style="lattice")
plotIndiv(nutri.sgcca,ylim=list(c(1,5),c(4,5),c(3,6)),style="ggplot2")
