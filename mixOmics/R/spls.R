# Author : F.Rohart
# created 22-04-2015
# last modified 22-04-2015
#
# perform the meta.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects


wrapper.spls <- function(X, Y, ncomp = 2, mode = c("regression", "canonical", "invariant", "classic"), study,
keepX , keepY ,keepX.constraint,keepY.constraint, max.iter = 500, tol = 1e-06, near.zero.var = FALSE,scale = TRUE)
{
    
    
    #-- validation des arguments --#
    # most of the checks are done in the wrapper.meta.spls.hybrid function
    X = as.matrix(X)
    Y = as.matrix(Y)

    if (!is.numeric(X) || !is.numeric(Y))
    stop("'X' and/or 'Y' must be a numeric matrix.")

    if(missing(keepX))
    keepX=rep(ncol(X), ncomp)

    if(missing(keepY))
    keepY=rep(ncol(Y), ncomp)

    result <- wrapper.meta.spls.hybrid(X=X,Y=Y,ncomp=ncomp,scale=scale,near.zero.var=near.zero.var,study=study,
    keepX=keepX,keepY=keepY,keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,max.iter=max.iter,tol=tol)
    
    
    cl = match.call()
    cl[[1]] = as.name("spls")
    result$call = cl
    
    class(result) = "spls"
    return(invisible(result))
    
    
    
    
    
}
