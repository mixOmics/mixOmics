#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
library(mixOmics)

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
plotIndiv(spca.rat, cex = 0.01,
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
    liver.spca<- spca(liver.toxicity$gene, keepX = rep(50, 3))
    plotIndiv(liver.spca)
    liver.spca<- spca(liver.toxicity$gene, keepX = rep(20, 3))
    plotIndiv(liver.spca)
    liver.spca<- spca(liver.toxicity$gene, keepX = rep(5, 3))
    plotIndiv(liver.spca)
}