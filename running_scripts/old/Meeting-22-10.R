library(ellipse)
library(mixOmics)

#########Meeting 22/10

########Issue 41

#### bug to fix in plotIndiv with graphics and Explained variance for PCA (only): 'PC1: 11.36 %'

data(liver.toxicity)

X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3,
                      keepX = c(50, 50, 50), keepY = c(10, 10, 10))


col <- c("blue", "red", "darkgreen", "darkviolet")
pch <- c(15, 16, 17, 18)

cex <- c(1, 1.5, 2, 2.5)
groupe=rep(c("groupe 1","groupe 2","groupe 3","groupe 4"),each = 16)

plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=TRUE,style="graphics",group=groupe,main='title') 
plotIndiv(toxicity.spls,add.legend=TRUE,ind.names=TRUE,style="graphics",group=groupe,main='title',plot.ellipse=T)


data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
plotIndiv(pca.res,style="graphics",main='title')

plotIndiv(pca.res,style="ggplot2",main='title')
plotIndiv(pca.res,style="ggplot2")


##Fix bug plot.pca
# FB can you make sure the 'main' title appears? thks

plot(pca.count, main = 'OTU raw count')


###issue 37

###plotContrib pour la PLS et autres methodes ### not push on bitbucket

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

plotContrib(toxicity.spls,group=liver.toxicity$treatment[, 'Time.Group'])
plotContrib(toxicity.spls,blocks='Y',group=liver.toxicity$treatment[, 'Time.Group'])
plotContrib(toxicity.spls,blocks='Y',group=liver.toxicity$treatment[, 'Time.Group'],legend.color = c(1:4))


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

plotContrib(nutrimouse.sgcca, group = nutrimouse$diet,  blocks = "lipid")##not sort but correct in selectVar
plotContrib(nutrimouse.sgcca, group = nutrimouse$diet,  blocks = "lipid", 
          main = 'my plot')

Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)

nutrimouse.sgccda1 <- wrapper.sgccda(blocks = data,
                                     Y = Y,
                                     design = design1,
                                     ncomp = c(2, 2),
                                     keep = list(c(10,10), c(15,15)),
                                     scheme = "centroid",
                                     verbose = FALSE,
                                     bias = FALSE)

plotContrib(nutrimouse.sgccda1, blocks = 2) ##not sort but correct in selectVar


######issue 37 and 43


###check plotContrib:OK

### plotArrow

data(liver.toxicity)

X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3,
                      keepX = c(50, 50, 50), keepY = c(10, 10, 10))


col <- c("blue", "red", "darkgreen", "darkviolet")
pch <- c(15, 16, 17, 18)

cex <- c(1, 1.5, 2, 2.5)



plotArrow(toxicity.spls, comp = c(1,2), ind.names = rep(c("a","b","c","d"),each=16),abline.line=FALSE)


plotArrow(toxicity.spls, comp = c(2,1),ind.names=FALSE,pch=pch,
          group=rep(c(1,2,3,4), each=16))


plotArrow(toxicity.spls,cex=c(1, 1.5, 2, 2.5),
          ind.names=FALSE,
          group=rep(c(1,2,3,4), each=16))



plotArrow(toxicity.spls, comp = c(1,3),
          ind.names=TRUE)
plotArrow(toxicity.spls, comp = c(1,3),label=TRUE) ## not push


groupe=rep(c("groupe 1","groupe 2","groupe 3","groupe 4"),each = 16)

plotArrow(toxicity.spls,add.legend=TRUE,abline.line=FALSE,group=groupe)

plotArrow(toxicity.spls,add.legend=TRUE,ind.names=FALSE,col=col,group=groupe) 


plotArrow(toxicity.spls,add.legend=TRUE,ind.names=FALSE,group=groupe)


plotArrow(toxicity.spls,ind.names=FALSE,group=groupe,abline.line = FALSE)
abline(0,0)
points(5,2)



data(nutrimouse)

X <- nutrimouse$lipid
Y <- nutrimouse$gene

col=c("red","blue")
pch=c(0,1)
cex=c(1,2)

nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)

plotArrow(nutri.res,group=rep(c("groupe 1", "groupe 2"),each=20))
plotArrow(nutri.res,col=col,group=rep(c("groupe 1", "groupe 2"),each=20))
plotArrow(nutri.res,col=col,group=rep(c("groupe 1", "groupe 2"),each=20),ind.names = T) 

nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 10)
plotArrow(nutri.res,col=col,add.legend=TRUE,
          group=rep(c("groupe 1", "groupe 2"),each=20)) 

nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 10, lambda2 = 0.008)
plotArrow(nutri.res,main="plot",comp=c(2,3))


plotArrow(nutri.res)
abline(0,0)
plotArrow(nutri.res,add.legend=T)
points(0,0)



diet = unmap(nutrimouse$diet)
blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = c(3, 3, 1))

#blocks
plotArrow(nutri.sgcca) #### 
plotArrow(nutri.sgcca,  ind.names = TRUE) 

plotArrow(nutri.sgcca, ind.names = nutrimouse$diet, 
          group=nutrimouse$diet,add.legend=TRUE)


plotArrow(nutri.sgcca,  ind.names = FALSE,add.legend=TRUE)
plotArrow(nutri.sgcca,  ind.names = 1:40,col="blue",abline.line=FALSE)
plotArrow(nutri.sgcca, 
          col="blue",abline.line=FALSE,add.legend=TRUE)




plotArrow(toxicity.spls,ylim=c(1,5))
plotArrow(nutri.sgcca,xlim=c(1,5))
