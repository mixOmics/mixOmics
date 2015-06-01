# Author : F.Rohart
# created 18-08-2014
# last modified 18-08-2014
#
# perform the meta.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects


wrapper.meta.block.spls <- function(X,
Y,
indY,
study,
ncomp=rep(2,length(X)),
keepX.constraint,
keepY.constraint,
keepX,
keepY,
design,
scheme,
mode,
scale = TRUE,
bias,
init ,
tol = 1e-06,
verbose,
max.iter = 500,
near.zero.var = FALSE)
{
    
    
    result <- wrapper.sparse.meta.block(X=X,Y=Y,indY=indY,study=study,ncomp=ncomp,keepX.constraint=keepX.constraint,
    keepY.constraint=keepY.constraint,keepX=keepX,keepY=keepY,design=design,scheme=scheme,mode=mode,scale=scale,
    bias=bias,init=init,tol=tol,verbose=verbose,max.iter=max.iter,near.zero.var=near.zero.var)
    
    
    
    cl = match.call()
    cl[[1]] = as.name("meta.block.spls")
    
    if(missing(indY))
    {
        keepX=result$keepA[-(length(X)+1)]
        keepY=result$keepA[length(X)+1][[1]]
        keepX.constraint=result$keepA.constraint[-(length(X)+1)]
        keepY.constraint=result$keepA.constraint[length(X)+1]
    }else{
        keepX=result$keepA[-indY]
        keepY=result$keepA[indY][[1]]
        keepX.constraint=result$keepA.constraint[-indY]
        keepY.constraint=result$keepA.constraint[indY][[1]]
    }
    
    out=list(call=cl,X=result$X,Y=result$Y[[1]],ncomp=result$ncomp,mode=result$mode,study=result$study,
    keepX=keepX,keepY=keepY,keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,
    variates=result$variates,loadings=result$loadings,variates.partial=result$variates.partial,loadings.partial=result$loadings.partial,
    names=result$names,tol=result$tol,iter=result$iter,nzv=result$nzv,scale=scale)
  
    if(!missing(ncomp))   out$ncomp=ncomp


    class(out) = "meta.block.spls"
    return(invisible(out))
    
}



