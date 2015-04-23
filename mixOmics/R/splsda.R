# Author : F.Rohart
# created 22-04-2015
# last modified 22-04-2015
#
# perform the meta.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects


wrapper.splsda <- function(X, Y, ncomp = 2, mode = c("regression", "canonical", "invariant", "classic"), study,
keepX = rep(ncol(X), ncomp),keepX.constraint=list(), max.iter = 500, tol = 1e-06, near.zero.var = FALSE,scale = TRUE)
{
    
    
    #-- validation des arguments --#
    # most of the checks are done in the wrapper.meta.spls.hybrid function
    X = as.matrix(X)
    
    if (is.null(dim(Y)))
    {
        Y = as.factor(Y)
    }  else {
        stop("'Y' should be a factor or a class vector.")
    }
    
    Y.mat=unmap(Y)

    result <- wrapper.meta.spls.hybrid(X=X,Y=Y.mat,ncomp=ncomp,scale=scale,near.zero.var=near.zero.var,study=study,
    keepX=keepX,keepX.constraint=keepX.constraint,max.iter=max.iter,tol=tol)
    
    
    cl = match.call()
    cl[[1]] = as.name("splsda")
    result$call = cl
    result$Y=Y
    result$ind.mat = Y.mat; result$names$Y = levels(Y)
    row.names(result$variates$Y) = row.names(X); row.names(result$loadings$Y) = paste0("Y", c(1 : nlevels(Y)))

    class(result) = "splsda"
    return(invisible(result))
    
    
    
    
    
}
