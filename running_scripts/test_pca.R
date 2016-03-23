#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

data(multidrug)

## this data set contains missing values, therefore
## the 'prcomp' function cannot be applied
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
par(mar=c(9,5,4,0))
plot(pca.res)
pca.res
par(mar=c(5.1, 4.1, 4.1, 2.1))

biplot(pca.res, xlabs = multidrug$cell.line$Class, cex = 0.7)

# samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
group = as.numeric(as.factor(multidrug$cell.line$Class)))


plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, group = as.numeric(as.factor(multidrug$cell.line$Class)),plot.ellipse=T)

plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, col = as.numeric(as.factor(multidrug$cell.line$Class)))

plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, col = rep(1:2,30))
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, group = as.numeric(as.factor(multidrug$cell.line$Class)),col = rep(1:10,6),plot.ellipse=TRUE)


plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
group = as.numeric(as.factor(multidrug$cell.line$Class)),main="bla")


plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
group = as.numeric(as.factor(multidrug$cell.line$Class)),add.legend=TRUE)


plotIndiv(pca.res, ind.names = F,#multidrug$cell.line$Class,
group = as.numeric(as.factor(multidrug$cell.line$Class)),add.legend=TRUE)


plotIndiv(pca.res, cex = 0.2,
    col = as.numeric(as.factor(multidrug$cell.line$Class)),style="3d")

# variables representation
plotVar(pca.res)

plotVar(pca.res, rad.in = 0.5, cex = 0.5,style="3d")


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################
if(additional.test==TRUE)
{
    
    plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,    group = as.numeric(as.factor(multidrug$cell.line$Class)),style="lattice",plot.centroid=TRUE,plot.ellipse=TRUE)
    
    plotIndiv(pca.res, style="graphics")
    plotIndiv(pca.res, style="ggplot2")


    plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
    group = as.numeric(as.factor(multidrug$cell.line$Class)),xlim=c(-10,10))


    plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
    group = as.numeric(as.factor(multidrug$cell.line$Class)),plot.centroid=TRUE)
    
    plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
    group = as.numeric(as.factor(multidrug$cell.line$Class)),plot.ellipse=TRUE)


    plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
    group = as.numeric(as.factor(multidrug$cell.line$Class)),main="bla")
    
    
    plotIndiv(pca.res, ind.names = multidrug$cell.line$Class,
    group = as.numeric(as.factor(multidrug$cell.line$Class)),add.legend=TRUE)
    
    
    plotIndiv(pca.res, ind.names = F,#multidrug$cell.line$Class,
    group = as.numeric(as.factor(multidrug$cell.line$Class)),add.legend=TRUE)

}
par(opar)
