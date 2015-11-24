require(lattice)
require(ellipse)
require(ggplot2)

library(mixOmics)
library(R.utils)

sourceDirectory(path='C:/Users/Yan/Downloads/mixOmics')



############# 3d

####### SPLS
data(liver.toxicity)

X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3,
                      keepX = c(30, 30, 30), keepY = c(10, 10, 10))

plotVar(toxicity.spls, style = "3d")

# test comp
plotVar(toxicity.spls,comp=c(1,3,3), style = "3d")
plotVar(toxicity.spls,comp=c(3,2,1), style = "3d")
plotVar(toxicity.spls,comp=c(1,2,3), style = "3d")

# test X.var.names, Y.var.names
plotVar(toxicity.spls,var.names = c(TRUE,FALSE), style = "3d")
plotVar(toxicity.spls,var.names = c(FALSE,TRUE), style = "3d")

# test xlab, ylab
plotVar(toxicity.spls,X.label="alpha",Y.label="beta",Z.label="gamma", style = "3d")

# test col
plotVar(toxicity.spls,col=c("green","red"), style = "3d")
col.x=rep(c("blue","red"),each=1558)
col.y=rep("green",10)
plotVar(toxicity.spls,col=list(col.x,col.y),style = "3d")
plotVar(toxicity.spls,col=list(col.x,col.y),X.var.names = TRUE, style = "3d")
plotVar(toxicity.spls,col=list(col.x,col.y),Y.var.names = TRUE, style = "3d")

# test cex
plotVar(toxicity.spls,cex=c(0.5,1.5), style = "3d")
cex.x=rep(c(0.8,1.2),1558)
cex.y=rep(1.5,10)
plotVar(toxicity.spls,cex=list(cex.x,cex.y), style = "3d")
plotVar(toxicity.spls,col=list(col.x,col.y),cex=list(cex.x,cex.y), style = "3d")

# test font
plotVar(toxicity.spls,font=c(1,2),var.names = c(FALSE,TRUE), style = "3d")
plotVar(toxicity.spls,font=c(1,3),var.names = c(FALSE,TRUE), style = "3d")
plotVar(toxicity.spls,font=c(2,2),var.names = c(TRUE,FALSE), style = "3d")
plotVar(toxicity.spls,font=c(4,2),var.names = c(TRUE,FALSE), style = "3d")
# test rad.in
plotVar(toxicity.spls,rad.in=0.1, style = "3d")
plotVar(toxicity.spls,rad.in=0.9, style = "3d")
plotVar(toxicity.spls,rad.in=0, style = "3d")
plotVar(toxicity.spls,rad.in=1, style = "3d")



# test add.legend
plotVar(toxicity.spls,add.legend=FALSE, style = "3d")
plotVar(toxicity.spls,add.legend=TRUE, style = "3d")
plotVar(toxicity.spls,add.legend=list(points=list(col=c("orange","cyan")),title="Protides",pch=c(0,1)),style="3d")


########### RCC
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.res <- rcc(X, Y)

plotVar(linn.res, style = "3d")

# test comp
plotVar(linn.res,comp=c(1,2), style = "3d")
plotVar(linn.res,comp=c(2,1), style = "3d")

# test col
plotVar(linn.res,col=c("green","red"), style = "3d")
col.x=rep(c("blue","red"),c(1,2))
col.y=rep("green",3)
plotVar(linn.res,col=list(col.x,col.y), style = "3d")
plotVar(linn.res,col=list(col.x,col.y),X.var.names = TRUE, style = "3d")
plotVar(linn.res,col=list(col.x,col.y),Y.var.names = TRUE, style = "3d")



######### SPCA

data(liver.toxicity)
rat.spca <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))

plotVar(rat.spca, style = "3d")
plotVar(rat.spca,add.legend = FALSE, style = "3d")
plotVar(rat.spca,col=list(rep(c("green","red"),each=1558)),style="3d")


######### SGCCA
data(nutrimouse)
diet = unmap(nutrimouse$diet)
blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = c(3, 3, 3))

plotVar(nutri.sgcca, style = "3d")

# test comp
plotVar(nutri.sgcca,comp=c(1,3,3), style = "3d")
plotVar(nutri.sgcca,comp=c(2,3,1), style = "3d")
plotVar(nutri.sgcca,comp=c(3,2,1), style = "3d")

#test blocks
plotVar(nutri.sgcca, blocks = c(1,3), style = "3d")
plotVar(nutri.sgcca, blocks= 2:1, style = "3d")
plotVar(nutri.sgcca,comp=c(2,1,3), blocks= c("diet","lipid"), style = "3d")
plotVar(nutri.sgcca,blocks=2:3, style = "3d")
plotVar(nutri.sgcca,blocks=1:3, style = "3d")

#test var.names
plotVar(nutri.sgcca,var.names = c(FALSE,TRUE), style = "3d")
plotVar(nutri.sgcca,var.names = c(FALSE,TRUE),blocks=2:3, style = "3d")
plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,TRUE),blocks=1:3, style = "3d")

)


