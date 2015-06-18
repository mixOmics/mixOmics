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
scheme, #block
scale,
bias,
near.zero.var=FALSE,
mode,
tol= 1e-06,
max.iter=500,
verbose=FALSE)

{
    if(is.list(X) & !is.data.frame(X))# either rgcca, sgcca,sgcca-DA, meta.block, meta.block-DA
    {
        
        #need to check if Y or indY is a factor, unmap it and then do the checks (no other factors etc)
        if((missing(indY)& missing(Y)) & is.null(tau))
        stop("Either 'Y', 'indY' or 'tau' is needed")
        
        if(is.null(tau)) # SGCCA/meta
        {
            
            isfactorY=FALSE


            if(!missing(Y))
            {
                if(is.list(Y) & !is.data.frame(X)) stop("Y must be a matrix or a factor")
            
                if (is.factor(Y)) {
                    #Y = as.factor(Y)
                    isfactorY=TRUE
                }
                
            }else if(!missing(indY))
            {
                temp=X[[indY]] #not called Y to not be an input of the wrappers
                if (is.factor(temp)) {
                    #temp = as.factor(temp)
                    isfactorY=TRUE
                }
            }else if(missing(indY))
            {
                stop("Either 'Y' or 'indY' is needed")
                
            }
            


            
            if(isfactorY)# either block.plsda/block.splsda/meta.block.plsda/meta.block.splsda
            {
                
                if(missing(keepX) & missing(keepX.constraint))
                {
                    if(missing(study)) #block.plsda
                    {
                        if(missing(scale)) scale=FALSE
                        message("a block Partial Least Squares - Discriminant Analysis is being performed (block PLS-DA)")
                        res=wrapper.block.plsda(X=X, Y=Y, indY=indY, ncomp = ncomp,design=design,scheme=scheme,
                        mode = mode,scale=scale,bias=bias,init=init,tol=tol,verbose=verbose,max.iter=max.iter,near.zero.var=near.zero.var)
                        
                    }else{# meta.block.plsda
                        if(missing(scale)) scale=FALSE
                        message("a meta block Partial Least Squares - Discriminant Analysis is being performed (meta.block.PLS-DA)")
                        res=wrapper.meta.block.plsda(X=X, Y=Y, indY=indY,study=study, ncomp = ncomp,design=design,scheme=scheme,
                        mode = mode,scale=scale,bias=bias,init=init,tol=tol,verbose=verbose,max.iter=max.iter,near.zero.var=near.zero.var)
                    }
                    
                    
                }else{
                    if(missing(study))# block.splsda
                    {
                        if(missing(scale)) scale=FALSE
                        message("a block sparse Partial Least Squares - Discriminant Analysis is being performed (block.sPLS-DA)")
                        res=wrapper.meta.block.splsda(X=X, Y=Y, indY=indY, ncomp = ncomp,keepX=keepX,keepX.constraint=keepX.constraint,
                        design=design,scheme=scheme,mode= mode,scale=scale,bias=bias,init=init,tol=tol,verbose=verbose,
                        max.iter=max.iter,near.zero.var=near.zero.var)
                        
                        
                    }else{# meta.block.splsda
                        if(missing(scale)) scale=FALSE
                        message("a meta block sparse Partial Least Squares - Discriminant Analysis is being performed (meta.block.sPLS-DA)")
                        res=wrapper.meta.block.splsda(X=X, Y=Y, indY=indY, ncomp = ncomp,study=study,keepX=keepX,keepX.constraint=keepX.constraint,
                        design=design,scheme=scheme,mode= mode,scale=scale,bias=bias,init=init,tol=tol,verbose=verbose,
                        max.iter=max.iter,near.zero.var=near.zero.var)
                    }
                    
                }
                
            }else{ # either block.pls/block.spls/meta.block.pls/meta.block.spls
                
                
                if(missing(keepX) & missing(keepX.constraint))
                {
                    if(missing(study)) #block.pls
                    {
                        if(missing(scale)) scale=FALSE
                        message("a block Partial Least Squares is being performed (block PLS)")
                        res=wrapper.block.pls(X=X, Y=Y, indY=indY, ncomp = ncomp,design=design,scheme=scheme,
                        mode = mode,scale=scale,bias=bias,init=init,tol=tol,verbose=verbose,max.iter=max.iter,near.zero.var=near.zero.var)
                        
                    }else{# meta.block.pls
                        if(missing(scale)) scale=FALSE
                        message("a meta block Partial Least Squares is being performed (meta.block.PLS)")
                        res=wrapper.meta.block.pls(X=X, Y=Y, indY=indY,study=study, ncomp = ncomp,design=design,scheme=scheme,
                        mode = mode,scale=scale,bias=bias,init=init,tol=tol,verbose=verbose,max.iter=max.iter,near.zero.var=near.zero.var)
                    }
                    
                    
                }else{
                    if(missing(study))# block.spls
                    {
                        if(missing(scale)) scale=FALSE
                        message("a block sparse Partial Least Squares is being performed (block.sPLS)")
                        res=wrapper.block.spls(X=X, Y=Y, indY=indY, ncomp = ncomp,keepX=keepX,keepX.constraint=keepX.constraint,
                        design=design,scheme=scheme,mode= mode,scale=scale,bias=bias,init=init,tol=tol,verbose=verbose,
                        max.iter=max.iter,near.zero.var=near.zero.var)
                        
                        
                    }else{# meta.block.spls
                        if(missing(scale)) scale=FALSE
                        message("a meta block sparse Partial Least Squares is being performed (meta.block.sPLS)")
                        res=wrapper.meta.block.spls(X=X, Y=Y, indY=indY, ncomp = ncomp,study=study,keepX=keepX,keepX.constraint=keepX.constraint,
                        design=design,scheme=scheme,mode= mode,scale=scale,bias=bias,init=init,tol=tol,verbose=verbose,
                        max.iter=max.iter,near.zero.var=near.zero.var)
                        
                    }
                    
                }
                
            }
            
        }else{ # RGCCA
            
            
            if(!missing(study)) {message("'study' is not used")}
            
            if(missing(keepX) & missing(keepX.constraint)) #RGCCA
            {
                message("A RGCCA analysis is being performed")
                if(missing(scale)) scale=FALSE
                res=wrapper.rgcca(X=X,design=design,tau=tau,ncomp = ncomp,mode=mode,
                max.iter=max.iter,scheme = scheme,scale = scale,init = init,bias = bias,tol = tol,verbose = verbose)
                
            }else{ #sparse RGCCA
                if(missing(scale)) scale=FALSE
                message("A sparse RGCCA analysis is being performed")
                res=wrapper.sparse.rgcca(X=X,design=design,tau=tau,mode=mode,ncomp = ncomp,keepX=keepX,keepX.constraint=keepX.constraint,
                max.iter=max.iter,scheme = scheme,scale = scale,init = init,bias = bias,tol = tol,verbose = verbose)
                
                
            }
        }
        
        
        
        #end if(is.list(X))
    }else{#either pls,spls, plsda, splsda or meta. pls/spls/plsda/splsda
        if(missing(Y))
        stop("Y is missing")
        if(is.list(Y) & !is.data.frame(X)) stop("Y must be a matrix or a factor")
        #if(!is.numeric(Y) & ncol(Y)==1) Y=as.factor(Y)
        
        if(missing(mode)) mode="regression"
        #check for unused inputs (scheme, etc etc)
        if(!is.null(tau) | !missing(design) | !missing(init) | !missing(scheme) | !missing(bias) | !missing(verbose))
        {
            if(!is.null(tau)) {message("'tau' is not used")}
            if(!missing(design)) {message("'design' is not used")}
            if(!missing(init)) {message("'init' is not used")}
            if(!missing(scheme)) {message("'scheme' is not used")}
            if(!missing(bias)) {message("'bias' is not used")}
            if(!missing(verbose)) {message("'verbose' is not used")}
            stop("unused input parameters")
        }
        
        
        if(is.factor(Y))#either plsda, splsda
        {
            
            #Check.entry.pls.single(X, ncomp, keepX,keepX.constraint) # to have the warnings relative to X and Y, instead of blocks
            if(length(Y)!=nrow(X)) {stop("unequal number of rows in 'X' and 'Y'.")}
            if(missing(keepX) & missing(keepX.constraint) & missing(keepY) & missing(keepY.constraint))  #plsda, meta.plsda
            {
                if(missing(study))
                {
                    if(missing(scale)) scale=TRUE
                    message("a Partial Least Squares - Discriminant Analysis is being performed (PLS-DA)")
                    res=wrapper.plsda(X=X, Y=Y, ncomp = ncomp, mode = mode,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                }else{# meta
                    if(missing(scale)) scale=FALSE
                    message("a meta Partial Least Squares - Discriminant Analysis is being performed (meta.PLS-DA)")
                    res=wrapper.meta.plsda(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
                
            }else{#splsda, meta.splsda
                if(missing(study))
                {
                    if(missing(scale)) scale=TRUE
                    message("a sparse Partial Least Squares - Discriminant Analysis is being performed (sPLS-DA)")
                    res=wrapper.splsda(X=X, Y=Y, ncomp = ncomp, mode = mode, keepX=keepX,keepX.constraint=keepX.constraint,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                }else{# meta
                    if(missing(scale)) scale=FALSE
                    message("a meta sparse Partial Least Squares - Discriminant Analysis is being performed (meta.sPLS-DA)")
                    res=wrapper.meta.splsda(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,keepX=keepX,keepX.constraint=keepX.constraint,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
            }
            
        }else{ #pls or spls
            
            
            #Check.entry.pls(X, Y, ncomp, keepX, keepY,keepX.constraint,keepY.constraint) # to have the warnings relative to X and Y, instead of blocks
            
            if(missing(keepX) & missing(keepX.constraint) & missing(keepY) & missing(keepY.constraint))  #pls, meta.pls
            {
                if(missing(study))
                {
                    if(missing(scale)) scale=TRUE
                    message("a Partial Least Squares is being performed (PLS)")
                    res=wrapper.pls(X=X, Y=Y, ncomp = ncomp, mode = mode,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                }else{ # meta
                    if(missing(scale)) scale=FALSE
                    message("a meta Partial Least Squares is being performed (meta.PLS)")
                    res=wrapper.meta.pls(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,
                    max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
            }else{
                if(missing(study))
                {
                    if(missing(scale)) scale=TRUE
                    message("a sparse Partial Least Squares is being performed (sPLS)")
                    res=wrapper.spls(X=X, Y=Y, ncomp = ncomp, mode = mode, keepX=keepX,keepY=keepY,
                    keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,max.iter = max.iter, tol = tol,
                    near.zero.var = near.zero.var,scale = scale)
                }else{
                    if(missing(scale)) scale=FALSE
                    message("a meta sparse Partial Least Squares is being performed (meta.sPLS)")
                    res=wrapper.meta.spls(X=X, Y=Y, ncomp = ncomp, mode = mode, study=study,keepX=keepX,keepY=keepY,
                    keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,max.iter = max.iter, tol = tol,
                    near.zero.var = near.zero.var,scale = scale)
                    
                }
            }
            
            
        }
    }
    cl = match.call()
    res$call=cl
    class(res)=c("mixOmics",class(res))
    return(invisible(res))
    
    
}