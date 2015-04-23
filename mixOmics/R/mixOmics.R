# Author : F.Rohart
# created 18-08-2014
# last modified 18-08-2014
#
# perform the meta.pls on a subset of variables on one only dimension, deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses, not for pls/spls as they would be overlapping batch effects

mixOmics=function(  X,
                    Y,
                    indY, #only use if Y not provided
                    ncomp,
                    keepX, #sparse
                    keepX.constraint, #hybrid
                    keepY, #sparse
                    keepY.constraint, #hybrid
                    study, #meta
                    design, #block
                    lambda=NULL, #sgcca for soft_threshold (superseeded by keepX
                    tau,# rgcca
                    scheme= "centroid", #block
                    scale=FALSE,
                    bias=FALSE,
                    near.zero.var=FALSE,
                    mode="regression",
                    tol= 1e-06,
                    max.iter=500,
                    verbose=FALSE)

{
    
    if(missing(Y) & missing(indY))
    {stop("'Y' or 'indY' are needed; if you just have one dataset, please use the pca() function")}
    
    if(is.list(X))# either rgcca, sgcca, meta.block
    {
        message("A block analysis is being performed")
  
        result=meta.block.spls(A=list(X,Y),indY=indY,design=design,lambda=lambda,ncomp=ncomp,scheme = scheme, scale =scale,  bias = bias,
        init = "svd", tol = tol, verbose = verbose,
        mode = mode, sparse = sparse, max.iter = max.iter,study = study, keepA = keepX,
        keepA.constraint = keepX.constraint, near.zero.var = near.zero.var)
        
        
    }else{#either pls,spls, plsda, splsda
        if(missing(Y))
        stop("Y is missing")
        if(!is.numeric(Y)) Y=as.factor(Y)
        
        
        
        
        if(is.factor(Y))#either plsda, splsda
        {
            Check.entry.pls.single(X, ncomp, keepX,keepX.constraint) # to have the warnings relative to X and Y, instead of blocks
            if(length(Y)!=nrow(X)) {stop("unequal number of rows in 'X' and 'Y'.")}
            if(missing(keepX)) #plsda, meta.plsda
            {
                message("a Partial Least Squares - Discriminant Analysis is being performed (PLS-DA)")
                res=wrapper.plsda(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                
            }else{#splsda, meta.splsda
                message("a sparse Partial Least Squares - Discriminant Analysis is being performed (sPLS-DA)")
                res=wrapper.splsda(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,keepX=keepX,keepX.constraint=keepX.constraint,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
            }
            
            
        }else{
            Check.entry.pls(X, Y, ncomp, keepX, keepY,keepX.constraint,keepY.constraint) # to have the warnings relative to X and Y, instead of blocks

            if(missing(keepX)) #pls, meta.pls
            {
                message("a Partial Least Squares is being performed (PLS)")
                res=wrapper.pls(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,
                max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                
            }else{#spls, meta.spls
                
                message("a sparse Partial Least Squares is being performed (sPLS)")
                res=wrapper.spls(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,keepX=keepX,keepY=keepY,
                keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,max.iter = max.iter, tol = tol,
                near.zero.var = near.zero.var,scale = scale)


                
            }
            
        }
        
        
    }
    
    class(res)=c("mixOmics",class(res))
    return(invisible(res))
        
        
    }