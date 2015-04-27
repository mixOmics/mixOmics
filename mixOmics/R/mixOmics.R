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
tau=NULL,# rgcca, number between 0,1 or "optimal"
init,
scheme= "centroid", #block
scale=FALSE,
bias=FALSE,
near.zero.var=FALSE,
mode="regression",
tol= 1e-06,
max.iter=500,
verbose=FALSE)

{
    if(is.list(X))# either rgcca, sgcca, meta.block
    {
        
        if(is.null(tau)) #RGCCA/meta
        {
            message("A block analysis is being performed")
            
            result=meta.block.spls(A=list(X,Y),indY=indY,design=design,ncomp=ncomp,scheme = scheme,
            scale =scale,  bias = bias,init = init, tol = tol, verbose = verbose, tau = tau,
            mode = mode, sparse = sparse, max.iter = max.iter,study = study, keepA = keepX,
            keepA.constraint = keepX.constraint, near.zero.var = near.zero.var)
        }else{        #if(no tau) SGCCA/meta
            
            message("A block analysis is being performed")
            
            result=meta.block.spls(A=list(X,Y),indY=indY,design=design,ncomp=ncomp,scheme = scheme,
            scale =scale,  bias = bias,init = init, tol = tol, verbose = verbose, tau = tau,
            mode = mode, sparse = sparse, max.iter = max.iter,study = study, keepA = keepX,
            keepA.constraint = keepX.constraint, near.zero.var = near.zero.var)
            
        }
        
    }else{#either pls,spls, plsda, splsda
        if(missing(Y))
        stop("Y is missing")
        if(!is.numeric(Y)) Y=as.factor(Y)
        
        if(is.factor(Y))#either plsda, splsda
        {
            
            #Check.entry.pls.single(X, ncomp, keepX,keepX.constraint) # to have the warnings relative to X and Y, instead of blocks
            if(length(Y)!=nrow(X)) {stop("unequal number of rows in 'X' and 'Y'.")}
            if(missing(keepX)) #plsda, meta.plsda
            {
                if(missing(study))
                {
                    message("a Partial Least Squares - Discriminant Analysis is being performed (PLS-DA)")
                    res=wrapper.plsda(X=X, Y=Y, ncomp = ncomp, mode = mode,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                }else{# meta
                    message("a meta Partial Least Squares - Discriminant Analysis is being performed (meta.PLS-DA)")
                    res=wrapper.meta.plsda(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
                
            }else{#splsda, meta.splsda
                if(missing(study))
                {
                    message("a sparse Partial Least Squares - Discriminant Analysis is being performed (sPLS-DA)")
                    res=wrapper.splsda(X=X, Y=Y, ncomp = ncomp, mode = mode, keepX=keepX,keepX.constraint=keepX.constraint,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                }else{# meta
                    message("a meta sparse Partial Least Squares - Discriminant Analysis is being performed (meta.sPLS-DA)")
                    res=wrapper.meta.splsda(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,keepX=keepX,keepX.constraint=keepX.constraint,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
            }
            
        }else{ #pls or spls
            
            
            #Check.entry.pls(X, Y, ncomp, keepX, keepY,keepX.constraint,keepY.constraint) # to have the warnings relative to X and Y, instead of blocks
            
            if(missing(keepX)) #pls, meta.pls
            {
                if(missing(study))
                {
                    message("a Partial Least Squares is being performed (PLS)")
                    res=wrapper.pls(X=X, Y=Y, ncomp = ncomp, mode = mode,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                }else{ # meta
                    message("a meta Partial Least Squares is being performed (meta.PLS)")
                    res=wrapper.meta.pls(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
            }else{
                if(missing(study))
                {
                    message("a sparse Partial Least Squares is being performed (sPLS)")
                    res=wrapper.spls(X=X, Y=Y, ncomp = ncomp, mode = mode, keepX=keepX,keepY=keepY,
                    keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,max.iter = max.iter, tol = tol,
                    near.zero.var = near.zero.var,scale = scale)
                }else{
                    message("a meta sparse Partial Least Squares is being performed (meta.sPLS)")
                    res=wrapper.meta.spls(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,keepX=keepX,keepY=keepY,
                    keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,max.iter = max.iter, tol = tol,
                    near.zero.var = near.zero.var,scale = scale)
                    
                }
            }
            
            
        }
    }
    class(res)=c("mixOmics",class(res))
    return(invisible(res))
    
    
}