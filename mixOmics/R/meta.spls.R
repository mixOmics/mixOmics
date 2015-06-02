# Author : F.Rohart
# created 22-04-2015
# last modified 22-04-2015
#
# perform the meta.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects


wrapper.meta.spls <- function(X, Y, ncomp = 2, mode = c("regression", "canonical", "invariant", "classic"), study,
keepX=rep(ncol(X), ncomp), keepY=rep(ncol(Y), ncomp), keepX.constraint=NULL,keepY.constraint=NULL, max.iter = 500,
tol = 1e-06, near.zero.var = FALSE,scale = FALSE)
{
    
    
    #-- validation des arguments --#
    # most of the checks are done in the wrapper.meta.spls.hybrid function
    X = as.matrix(X)
    Y = as.matrix(Y)

    result <- wrapper.meta.spls.hybrid(X=X,Y=Y,ncomp=ncomp,scale=scale,near.zero.var=near.zero.var,study=study,mode=mode,
    keepX=keepX,keepY=keepY,keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,max.iter=max.iter,tol=tol)
    
    
    cl = match.call()
    cl[[1]] = as.name("meta.spls")
    
    out=list(call=cl,X=result$X[[1]],Y=result$Y[[1]],ncomp=result$ncomp,study=result$study,mode=result$mode,keepX=result$keepA[[1]],keepY=result$keepA[[2]],
    keepX.constraint=result$keepA.constraint[[1]],keepY.constraint=result$keepA.constraint[[2]],
    variates=result$variates,loadings=result$loadings,variates.partial=result$variates.partial,loadings.partial=result$loadings.partial,
    names=result$names,tol=result$tol,iter=result$iter,nzv=result$nzv,scale=scale)
    
    class(out) = "meta.spls"
    return(invisible(out))
 
    
    
    
    
}
