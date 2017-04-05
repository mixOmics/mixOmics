#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

data(nutrimouse)
X = nutrimouse$gene
spls.nut = spls(X=X, Y= X, ncomp = 3)#, keepX = rep(10, 3), keepY = rep(10, 3))



acp = pca(X, scale=T)
sacp = spca(X)
ind.match = names(acp)%in% names(sacp)
ind.match2 = match(names(acp)[ind.match], names(sacp))
all.equal(acp[ind.match],sacp[ind.match2])



data(liver.toxicity)
spca.rat <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))
spca.rat

## variable representation
plotVar(spca.rat,style="3d")
#rgl.close()

plotVar(spca.rat, cex = 3)

## samples representation
plotIndiv(spca.rat, ind.names = liver.toxicity$treatment[, 3],
group = as.numeric(liver.toxicity$treatment[, 3]))
plotIndiv(spca.rat, cex = 0.7,
col = as.numeric(liver.toxicity$treatment[, 3]),style="3d")



#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    
    data(liver.toxicity)
    
    #test ncomp
    liver.spca<- spca(liver.toxicity$gene, ncomp = 1)
    #plotIndiv(liver.spca)
    liver.spca<- spca(liver.toxicity$gene, ncomp = 2)
    plotIndiv(liver.spca)
    liver.spca<- spca(liver.toxicity$gene, ncomp = 3)
    plotIndiv(liver.spca)
    plotIndiv(liver.spca,comp=c(1,3))
    #liver.spca<- spca(liver.toxicity$gene, ncomp = 70) # error ncomp>64
    
    #test scale
    liver.spca<- spca(liver.toxicity$gene, scale=TRUE)
    plotIndiv(liver.spca)
    liver.spca<- spca(liver.toxicity$gene, scale=FALSE)
    plotIndiv(liver.spca)
    
    #test center
    liver.spca<- spca(liver.toxicity$gene, center=TRUE)
    plotIndiv(liver.spca)
    liver.spca<- spca(liver.toxicity$gene, center=FALSE)
    plotIndiv(liver.spca)
    
    #test keepX
    liver.spca<- spca(liver.toxicity$gene, keepX = rep(50, 2))
    plotIndiv(liver.spca)
    liver.spca<- spca(liver.toxicity$gene, keepX = rep(20, 2))
    plotIndiv(liver.spca)
    liver.spca<- spca(liver.toxicity$gene, keepX = rep(5, 2))
    plotIndiv(liver.spca)
}
par(opar)
