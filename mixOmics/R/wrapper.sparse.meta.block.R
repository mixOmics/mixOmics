# Author : F.Rohart
# created 18-08-2014
# last modified 18-08-2014
#
# perform the meta.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects


wrapper.sparse.meta.block <- function(  X,
Y,
indY,
study,
ncomp,
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
    
    
    if(missing(scheme)) scheme= "centroid"
    if(missing(bias)) bias= FALSE
    if(missing(verbose)) verbose= FALSE
    if(missing(mode)) mode="regression"
    #print("bla")
    
    #need to check if Y or indY is a factor, unmap it and then do the checks (no other factors etc)
    #print(missing(ncomp))

    
    check=Check.entry.wrapper.sparse.meta.block(X=X,Y=Y,indY=indY,ncomp=ncomp,keepX=keepX,
    keepX.constraint=keepX.constraint,keepY=keepY,keepY.constraint=keepY.constraint,
    study=study,design=design,init=init,scheme=scheme,
    scale=scale,bias=bias,near.zero.var=near.zero.var,mode=mode,tol=tol,
    max.iter=max.iter,verbose=verbose)
    #   print("bla2")
    
    A=check$A
    indY=check$indY
    study=check$study
    design=check$design
    ncomp=check$ncomp
    keepA=check$keepA
    keepA.constraint=check$keepA.constraint
    init=check$init
    nzv.A=check$nzv.A

    #print(missing(ncomp))

    #  print(ncomp)
    #print(length(A))
    #message("A block analysis is being performed")
    
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
    
    
    result=sparse.meta.block(A=A,indY=indY,design=design,ncomp=ncomp,scheme = scheme,
    scale =scale,  bias = bias,init = init, tol = tol, verbose = verbose,tau=NULL,
    mode = mode, max.iter = max.iter,study = study, keepA = keepA,
    keepA.constraint = keepA.constraint)#, near.zero.var = near.zero.var)
    
    
    
    result$ncomp = ncomp
    if(near.zero.var)
    result$nzv=nzv.A

    class(result) = c("sparse.meta.block")
    return(invisible(result))
    

}



