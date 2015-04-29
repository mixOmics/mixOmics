# Author : F.Rohart
# created 18-08-2014
# last modified 18-08-2014
#
# perform the meta.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects


wrapper.block.splsda <- function(X,
Y,
indY,
ncomp=rep(2,length(X)),
keepX.constraint,
keepX,
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
    if(!missing(Y))
    {
        if (is.null(dim(Y))) {
            Y = as.factor(Y)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        Y.input=Y
        Y=unmap(Y)
        
    }else if(!missing(indY))
    {
        temp=X[[indY]] #not called Y to not be an input of the wrapper.sparse.meta.block
        if (is.null(dim(temp))) {
            temp = as.factor(temp)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        Y.input=temp
        X[[indY]]=unmap(temp)
    }else if(missing(indY))
    {
        stop("Either 'Y' or 'indY' is needed")
        
    }

    
    result <- wrapper.sparse.meta.block(X=X,Y=Y,indY=indY,ncomp=ncomp,keepX.constraint=keepX.constraint,
    keepX=keepX,design=design,scheme=scheme,mode=mode,scale=scale,
    bias=bias,init=init,tol=tol,verbose=verbose,max.iter=max.iter,near.zero.var=near.zero.var)
    
    
    
    
    cl = match.call()
    cl[[1]] = as.name("block.splsda")
    
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
    
    
    out=list(call=cl,X=result$X,Y=Y.input,ind.mat=result$Y[[1]],ncomp=result$ncomp,mode=result$mode,keepX=keepX,keepY=keepY,
    keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,
    variates=result$variates,loadings=result$loadings,names=result$names,
    tol=result$tol,iter=result$iter,nzv=result$nzv,scale=scale)
    
    if(!missing(ncomp))   out$ncomp=ncomp

    class(out) = "block.splsda"
    return(invisible(out))
    
}



