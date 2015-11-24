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



plotArrow(toxicity.spls, comp = c(1,2), ind.names = rep(c("a","b","c","d"),each=16),abline.line=FALSE)


plotArrow(toxicity.spls, comp = c(2,1),ind.names=FALSE,pch=pch,
          group=rep(c(1,2,3,4), each=16))


plotArrow(toxicity.spls,cex=c(1, 1.5, 2, 2.5),
          ind.names=FALSE,
          group=rep(c(1,2,3,4), each=16))



plotArrow(toxicity.spls, comp = c(1,3))

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
plotArrow(nutri.sgcca,notblocks="diet") #### 
plotArrow(nutri.sgcca,  ind.names = TRUE,notblocks="diet") 

plotArrow(nutri.sgcca, notblocks="diet", ind.names = nutrimouse$diet, 
          group=nutrimouse$diet,add.legend=TRUE)


plotArrow(nutri.sgcca, notblocks="diet", ind.names = FALSE,add.legend=TRUE)
plotArrow(nutri.sgcca, notblocks="diet", ind.names = 1:40,col="blue",abline.line=FALSE)
plotArrow(nutri.sgcca, notblocks="diet",
          col="blue",abline.line=FALSE,add.legend=TRUE)




plotArrow(toxicity.spls,ylim=c(1,5))
plotArrow(nutri.sgcca,notblocks="diet",xlim=c(1,5))
