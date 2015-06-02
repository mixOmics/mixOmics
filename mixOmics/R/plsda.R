# Author : F.Rohart
# created 22-04-2015
# last modified 22-04-2015
#
# perform the meta.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects


wrapper.plsda <- function(X, Y, ncomp = 2, mode = c("regression", "canonical", "invariant", "classic"),
max.iter = 500, tol = 1e-06, near.zero.var = FALSE,scale = TRUE)
{
    
    
    #-- validation des arguments --#
    # most of the checks are done in the wrapper.meta.spls.hybrid function
    
    if (is.null(Y))
    stop("'Y' has to be something else than NULL.")
    
    if (is.null(dim(Y)))
    {
        Y = as.factor(Y)
    }  else {
        stop("'Y' should be a factor or a class vector.")
    }
    
    Y.mat=unmap(Y)
    colnames(Y.mat) = paste0("Y", 1:ncol(Y.mat))

    X = as.matrix(X)

    result <- wrapper.meta.spls.hybrid(X=X,Y=Y.mat,ncomp=ncomp,scale=scale,near.zero.var=near.zero.var,mode=mode,
    max.iter=max.iter,tol=tol)
        
    cl = match.call()
    cl[[1]] = as.name("plsda")

    out=list(call=cl,X=result$X[[1]],Y=Y,ind.mat=result$Y[[1]],ncomp=result$ncomp,mode=result$mode,variates=result$variates,loadings=result$loadings,
        names=result$names,tol=result$tol,iter=result$iter,nzv=result$nzv,scale=scale)
        out$names$Y = levels(Y)
        row.names(out$variates$Y) = row.names(out$variates$X)
        row.names(out$loadings$Y) = paste0("Y", c(1 : nlevels(Y)))
   
    class(out) = "plsda"
    return(invisible(out))
    



}
