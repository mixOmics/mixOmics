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
mode,
scale=FALSE,
near.zero.var=FALSE,
max.iter = 500,
tol = 1e-06)
{
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || length(ncomp)>1)
    stop("invalid number of variates, 'ncomp'.")
    
    
    
    #-- validation des arguments --#
   
    #set the default study factor
    if(missing(study))
    {
        study=factor(rep(1,nrow(X)))
    }else{
        study=as.factor(study)
    }
    if(length(study)!=nrow(X)) stop(paste0("'study' must be a factor of length ",nrow(X),"."))
    
    if(any(table(study)<=1)) stop("At least one study has only one sample, please consider removing before calling the function again")
    if(any(table(study)<5)) warning("At least one study has less than 5 samples, mean centering might not do as expected")
    
    
    design = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)
    
    
    


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
    
    check=Check.entry.pls(X, Y, ncomp, keepX, keepY,keepX.constraint,keepY.constraint,mode,near.zero.var=near.zero.var) # to have the warnings relative to X and Y, instead of blocks
    X=check$X
    Y=check$Y
    ncomp=check$ncomp
    mode=check$mode
    keepX.constraint=check$keepX.constraint
    keepY.constraint=check$keepY.constraint
    keepX=check$keepX
    keepY=check$keepY
    nzv.A=check$nzv.A


    result <- sparse.meta.block(A = list(X = X, Y = Y), indY = 2, mode = mode, ncomp = c(ncomp, ncomp), tol = tol, max.iter = max.iter,
    design = design, keepA = list(keepX,keepY),keepA.constraint = list(keepX.constraint,keepY.constraint),
    scale = scale, scheme = "centroid",init="svd", study = study)#,near.zero.var=near.zero.var)
    
    result$ncomp = ncomp
    if(near.zero.var)
    result$nzv=nzv.A

    class(result) = c("meta.spls.hybrid")
    return(invisible(result))
    
}



