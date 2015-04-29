# Author : F.Rohart
# created 18-08-2014
# last modified 18-08-2014
#
# perform the meta.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects


wrapper.meta.spls.hybrid <- function(X,
Y,
study,
ncomp=2,
keepX.constraint,
keepY.constraint,
keepX,
keepY,
scale=FALSE,
near.zero.var=FALSE,
max.iter = 500,
tol = 1e-06)
{
    
 
    
    #-- validation des arguments --#
    # most of the checks are done in the meta.block.spls function
    X = as.matrix(X)
    Y = as.matrix(Y)
    
    if(missing(keepX.constraint))
    {
        if(missing(keepX))
        {
            keepX=rep(ncol(X),ncomp)
        }
        keepX.constraint=list()
    }else{
        if(length(keepX.constraint)>ncomp)
        stop(paste0("you should have length(keepX.constraint) lower or equal to ",ncomp,"."))

        if(missing(keepX))
        {
            keepX=rep(ncol(X),ncomp-length(keepX.constraint))
        }
    }
    
    if(missing(keepY.constraint))
    {
        if(missing(keepY))
        {
            keepY=rep(ncol(Y),ncomp)
        }
        keepY.constraint=list()
    }else{
        if(length(keepY.constraint)>ncomp)
        stop(paste0("you should have length(keepY.constraint) lower or equal to ",ncomp,"."))

       if(missing(keepY))
        {
            keepY=rep(ncol(Y),ncomp-length(keepY.constraint))
        }
    }

    #set the default study factor
    if(missing(study))
    {
        study=factor(rep(1,nrow(X)))
    }else{
        study=as.factor(study)
    }
    if(length(study)!=nrow(X)) stop(paste0("'study' must be a factor of length ",nrow(X),"."))
    
    
    design = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)
    
    # match keepX.constraint and the colnames of X in order for keepX.constraint to be a list of character
    # safety if keepX.constraint contains a mixed of character/numeric. It should one or the other, not a mix
    if(length(keepX.constraint)>0)
    {
        if(!is.numeric(unlist(keepX.constraint)))
        {
            ind=match(unlist(keepX.constraint),colnames(X))
            if(sum(is.na(ind))>0) stop("'keepX.constraint' must contains a subset of colnames(X) or the position of the X-variables you wish to keep.")
        }
        X.indice=X[,unlist(keepX.constraint),drop=FALSE]
        keepX.constraint=relist(colnames(X.indice),skeleton=keepX.constraint)
    }
    
    # same for keepY.constraint
    if(length(keepY.constraint)>0)
    {
        if(!is.numeric(unlist(keepY.constraint)))
        {
            ind=match(unlist(keepY.constraint),colnames(Y))
            if(sum(is.na(ind))>0) stop("'keepY.constraint' must contains a subset of colnames(Y) or the position of the Y-variables you wish to keep.")
        }
        Y.indice=Y[,unlist(keepY.constraint),drop=FALSE]
        keepY.constraint=relist(colnames(Y.indice),skeleton=keepY.constraint)
    }
    
    
    # we need numbers in keepX.constraint from now on
    keepX.constraint= lapply(keepX.constraint,function(x){match(x,colnames(X))})
    keepY.constraint= lapply(keepY.constraint,function(x){match(x,colnames(Y))})
    
    # we need a vector of length ncomp for keepA
    # update keepX and keepY
    #keepX=c(unlist(lapply(keepX.constraint,length)),keepX) #of length ncomp, can contains 0
    #keepY=c(unlist(lapply(keepY.constraint,length)),keepY) #of length ncomp, can contains 0
    
    # A: list of matrices
    # indY: integer, pointer to one of the matrices of A
    # design: design matrix, links between matrices. Diagonal must be 0
    # ncomp: vector of ncomp, per matrix
    # scheme: a function "g", refer to the article (thanks Benoit)
    # scale: do you want to scale ? mean is done by default and cannot be changed (so far)
    # bias: scale the data with n or n-1
    # init: one of "svd" or "random", initialisation of the algorithm
    # tol: nobody cares about this
    # verbose: show the progress of the algorithm
    # mode: canonical, classic, invariant, regression
    # max.iter: nobody cares about this
    # study: factor for each matrix of A, must be a vector
    # keepA: keepX of spls for each matrix of A. must be a list. Each entry must be of the same length (max ncomp)
    # keepA.constraint: keepX.constraint, which variables are kept on the first num.comp-1 components. It is a list of characters
    # near.zero.var: do you want to remove variables with very small variance
    
    check=Check.entry.pls(X, Y, ncomp, keepX, keepY,keepX.constraint,keepY.constraint) # to have the warnings relative to X and Y, instead of blocks
    X=check$X
    Y=check$Y
    ncomp=check$ncomp

    result <- sparse.meta.block(A = list(X = X, Y = Y), indY = 2, mode = "regression", ncomp = c(ncomp, ncomp), tol = tol, max.iter = max.iter,
    design = design, keepA = list(keepX,keepY),keepA.constraint = list(keepX.constraint,keepY.constraint),
    scale = scale, scheme = "centroid",init="svd", study = study,near.zero.var=near.zero.var)
    
    result$ncomp = ncomp

    class(result) = c("meta.spls.hybrid")
    return(invisible(result))
    
}



