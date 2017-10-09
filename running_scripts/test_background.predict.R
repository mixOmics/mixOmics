# created on 01/03/17
# Author: Florian Rohart
#purpose: test the background.predict function

opar <- par(no.readonly = TRUE)

#library(mixOmics)
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

splsda.breast <- splsda(X, Y,keepX=c(10,10),ncomp=2)

background = background.predict(splsda.breast,comp.predicted=2, dist = "centroids.dist")

plotIndiv(splsda.breast, background = background, style = "ggplot2")
plotIndiv(splsda.breast, background = background, style = "graphics")
plotIndiv(splsda.breast, background = background, style = "lattice")

#
plotIndiv(splsda.breast, background = background,legend=T, ellipse=T, title = "background not covering enough by default with ellipse")
# increase xlim and ylim
background = background.predict(splsda.breast,comp.predicted=2, dist = "centroids.dist", xlim = c(-6,7), ylim = c(-6,6))
plotIndiv(splsda.breast, background = background,legend=T, ellipse=T, title="now all good")

background = background.predict(splsda.breast,comp.predicted=2, dist = "centroids.dist", xlim = c(-6,7), ylim = c(-6,6), resolution = 500)
plotIndiv(splsda.breast, background = background,legend=T, ellipse=T, title="increasing resolution")



set.seed(123)

data(liver.toxicity)
X=liver.toxicity$gene
Y=as.factor(liver.toxicity$treatment[, 4])


plsda.liver <- plsda(X, Y, ncomp = 4)

background = background.predict(plsda.liver,comp.predicted=2, dist = "centroids.dist")
plotIndiv(plsda.liver, background = background,legend=T,title = "good background")
plotIndiv(plsda.liver, background = background,legend=T,ellipse=T, title = "background not covering enough by default with ellipse")

# increase limits
background = background.predict(plsda.liver,comp.predicted=2, dist = "centroids.dist", xlim = c(-70,80), ylim=c(-60,80))
plotIndiv(plsda.liver, background = background,legend=T, ellipse=T, title="now all good")
plotIndiv(plsda.liver, background = background,legend=T, ellipse=T, style="lattice")
plotIndiv(plsda.liver, background = background,legend=T, ellipse=T, style="graphics")


plotIndiv(plsda.liver, background = background,legend=T, col.per.group = c(2,3,4,1),title = "change col.per.group")
plotIndiv(plsda.liver, background = background,legend=T, col.per.group = c(2,4,1,3),title = "change col.per.group")


par(mfrow=c(2,3))
for(j in 1:2)
{
    out = background.predict(plsda.liver, dist = "max.dist", comp.predicted = j)
    plotIndiv(plsda.liver, background = out,legend=T,title = paste0("max.dist, pred.comp=",j))
    
    out = background.predict(plsda.liver, dist = "centroids.dist", comp.predicted = j)
    plotIndiv(plsda.liver, background = out,legend=T,title = paste0("centroids.dist, pred.comp=",j))
    
    out = background.predict(plsda.liver, dist = "mahalanobis.dist", comp.predicted = j)
    plotIndiv(plsda.liver, background = out,legend=T,title = paste0("mahalanobis.dist, pred.comp=",j))

}

par(opar)
par(mfrow=c(1,1))
