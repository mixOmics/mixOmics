#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

data(liver.toxicity)

# implement IPCA on a microarray dataset
ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode="deflation")
ipca.res

# samples representation
plotIndiv(ipca.res, ind.names = as.character(liver.toxicity$treatment[, 4]),
group = as.numeric(as.factor(liver.toxicity$treatment[, 4])))
cat("\n\n")

plotIndiv(ipca.res, cex = 0.5,col = as.numeric(as.factor(liver.toxicity$treatment[, 4])),style="3d")

plotIndiv(ipca.res, style="3d",col = as.numeric(as.factor(liver.toxicity$treatment[, 4])))

# variables representation
plotVar(ipca.res, cex = 0.5)
cat("\n\n")

plotVar(ipca.res, rad.in = 0.5, cex = 0.5,style="3d")


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    
    data(multidrug)
    data(liver.toxicity)
    
    #test ncomp
    #liver.ipca=ipca(multidrug$ABC.trans)###
    liver.ipca=ipca(liver.toxicity$gene)
    
    #liver.ipca=pca(liver.toxicity,ncomp=2)####
    liver.ipca=pca(liver.toxicity$gene,ncomp=2)
    plotIndiv(liver.ipca)
    plot(liver.ipca)
    liver.ipca=pca(liver.toxicity$gene,ncomp=63)
    plotIndiv(liver.ipca,comp=c(62,63))
    plot(liver.ipca)
    
    #test mode
    liver.ipca=ipca(liver.toxicity$gene,mode="deflation")
    plotIndiv(liver.ipca)
    
    liver.ipca=ipca(liver.toxicity$gene,mode="parallel") #??
    plotIndiv(liver.ipca)
    
    #test fun
    liver.ipca=ipca(liver.toxicity$gene,fun="logcosh")
    plotIndiv(liver.ipca)
    
    liver.ipca=ipca(liver.toxicity$gene,fun="exp")
    plotIndiv(liver.ipca)
    
    #test w.init
    mat=matrix(1:9,3,3)
    liver.ipca=ipca(liver.toxicity$gene,w.init=mat,ncomp=3)
    plotIndiv(liver.ipca)
    
    mat=matrix(9:17,3,3)
    liver.ipca=ipca(liver.toxicity$gene,w.init=mat,ncomp=3)
    plotIndiv(liver.ipca)
    
    mat=matrix(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),3,3)
    liver.ipca=ipca(liver.toxicity$gene,w.init=mat,ncomp=3)
    plotIndiv(liver.ipca)
    
    ## error: w.init is a zero matrix
    #mat=matrix(rep(0,9),3,3)
    #liver.ipca=ipca(liver.toxicity$gene,w.init=mat,ncomp=3) ###
    #plotIndiv(liver.ipca)
    
    #error: w.init contains NA values
    #mat=matrix(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),3,3)
    #mat[2,1]=NA
    #liver.ipca=ipca(liver.toxicity$gene,w.init=mat,ncomp=3)###
    #plotIndiv(liver.ipca)
    
}
par(opar)

