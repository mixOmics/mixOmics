#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2015
# last modified: 25-02-2016
#
# Copyright (C) 2015
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


# --------------------------------------
# Check.entry.single
# --------------------------------------
Check.entry.single = function(X,  ncomp, keepX, keepX.constraint,q)
{

    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop(paste0("'X[[",q,"]]' must be a numeric matrix."))
    
    X = as.matrix(X)
    
    if (!is.numeric(X) )
    stop(paste0("'X[[",q,"]]'  must be a numeric matrix."))
    
    N = nrow(X)
    P = ncol(X)
    
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop(paste0("invalid number of variates 'ncomp' for matrix 'X[[",q,"]]'."))
    
    ncomp = round(ncomp)
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]]=X.names
    }
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X)  = ind.names
    }
    
    if(length(unique(rownames(X)))!=nrow(X)) stop("samples should have a unique identifier/rowname")
    if(length(unique(X.names))!=P) stop("Unique indentifier is needed for the columns of X")
    
    return(list(X=X,ncomp=ncomp,X.names=X.names,ind.names=ind.names))
}



# --------------------------------------
# Check.entry.pls
# --------------------------------------

Check.entry.pls = function(X, Y, ncomp, keepX, keepY, keepX.constraint,keepY.constraint,mode,verbose, near.zero.var,max.iter,tol,logratio,DA,multilevel)
{
    
    if(missing(mode)) mode="regression"
    
    if(length(mode)>1) mode=mode[1]
    
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
    {stop("Choose one of the four following modes: canonical, invariant, classic or regression")}
    
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")
    
    X = as.matrix(X)
    
    if (!(logratio %in% c("none", "CLR")))
    {stop("Choose one of the two following logratio transformation: none or CLR")}
    
    if(!is.null(multilevel))
    {
        #multilevel analysis: withinVariation and then pls-like
        # f it's DA analysis, Y is ignored and we look in 'multilevel' input parameter
        if(DA)
        {
            if(!is.null(Y))
            {
                message("Multilevel Analysis will be performed based on 'multilevel' input, 'Y' is ignored")
            }
            Y=multilevel
        }else{
            if ((nrow(X) != nrow(multilevel)))
            stop("unequal number of rows in 'X' and 'multilevel'.")
            Y = as.matrix(Y)
            
            if (!is.numeric(X) || !is.numeric(Y))
            stop("'X' and/or 'Y' must be a numeric matrix.")
            
        }
    }else{
        Y = as.matrix(Y)
        
        if (!is.numeric(X) || !is.numeric(Y))
        stop("'X' and/or 'Y' must be a numeric matrix.")
        
    }
    

    N = nrow(X)
    Q = ncol(Y)
    P= ncol(X)
    
    if ((N != nrow(Y)))
    stop("unequal number of rows in 'X' and 'Y'.")
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || length(ncomp)>1)
    stop("invalid number of variates, 'ncomp'.")
    
    
    ncomp = round(ncomp)
    if(ncomp > P)
    {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", P, ".")
        ncomp = P
    }
    
    if(!is.numeric(tol) | tol<=0)
    stop("tol must be non negative")
    if(!is.numeric(max.iter) | max.iter<=0)
    stop("max.iter must be non negative")
    
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]]=X.names
    }
    
    if (dim(Y)[2] == 1) Y.names = "Y"
    if (dim(Y)[2] > 1)
    {
        Y.names = dimnames(Y)[[2]]
        if (is.null(Y.names))
        {
            Y.names = paste("Y", 1:Q, sep = "")
            dimnames(Y)[[2]]=Y.names
        }
    }
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = dimnames(Y)[[1]]
        rownames(X) = ind.names
    }
    
    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X) = rownames(Y) = ind.names
    }
    
    rownames(X) = rownames(Y) = ind.names
    
    if(length(unique(X.names))!=P) stop("Unique indentifier is needed for the columns of X")
    if(length(unique(Y.names))!=Q) stop("Unique indentifier is needed for the columns of Y")
    
    
    # check on keepX and keepX.constraint
    if(missing(keepX.constraint))
    {
        if(missing(keepX))
        {
            keepX=rep(P,ncomp)
        }else{
            if(length(keepX)<ncomp) {keepX=c(keepX,rep(P,ncomp-length(keepX)))} #complete the keepX already provided
        }
        keepX.constraint=list()
    }else{
        if(length(keepX.constraint)>ncomp)
        stop(paste0("you should have length(keepX.constraint) lower or equal to ",ncomp,"."))
        
        if(missing(keepX))
        {
            keepX=rep(P,ncomp-length(keepX.constraint))
        }else{
            
            if((length(keepX.constraint)+length(keepX)) <ncomp)
            {keepX=c(keepX,rep(P,ncomp-length(keepX)-length(keepX.constraint)))}
            if((length(keepX.constraint)+length(keepX))>ncomp) stop(paste0("length (keepX.constraint) + length(keepX) should be lower than ncomp"))
            
        }
    }
    
    # check on keepY and keepY.constraint
    if(missing(keepY.constraint))
    {
        if(missing(keepY))
        {
            keepY=rep(Q,ncomp)
        }else{
            if(length(keepY)<ncomp) {keepY=c(keepY,rep(Q,ncomp-length(keepY)))} #complete the keepY already provided
        }
        keepY.constraint=list()
    }else{
        if(length(keepY.constraint)>ncomp)
        stop(paste0("you should have length(keepY.constraint) lower or equal to ",ncomp,"."))
        
        if(missing(keepY))
        {
            keepY=rep(Q,ncomp-length(keepY.constraint))
        }else{
            
            if((length(keepY.constraint)+length(keepY)) <ncomp)
            {keepY=c(keepY,rep(Q,ncomp-length(keepY)-length(keepY.constraint)))}
            if((length(keepY.constraint)+length(keepY))>ncomp) stop(paste0("length (keepY.constraint) + length(keepY) should be lower than ncomp"))
        }
    }
    
    
    
    if (any(keepX<0))
    stop("each component of 'keepX' must be non negative ")
    if (any(keepY<0))
    stop("each component of 'keepY' must be non negative ")
    
    if (any(keepX > ncol(X)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    if (any(keepY > ncol(Y)))
    stop("each component of 'keepY' must be lower or equal than ", Q, ".")
    
    if (is.numeric(unlist(keepX.constraint)) && any(unlist(keepX.constraint) > ncol(X)))
    stop("each entry of 'keepX.constraint' must be lower or equal than ", P, ".")
    if ( is.numeric(unlist(keepY.constraint)) && any(unlist(keepY.constraint) > ncol(Y)))
    stop("each entry of 'keepY.constraint' must be lower or equal than ", Q, ".")
    
    
    if(!is.logical(near.zero.var))
    stop("near.zero.var must be either TRUE or FALSE")
    
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
    

    # at this stage keepA.constraint need to be character, to remove easily variables with near zero variance
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = nearZeroVar(X)
        
        if (length(nzv.A$Position) > 0)
        {
            names.remove.X=colnames(X)[nzv.A$Position]
            X = X[, -nzv.A$Position,drop=FALSE]
            #if (verbose)
            warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
            if(ncol(X) == 0) {stop("No more variables in X")}
            
            # at this stage, keepA.constraint need to be numbers
            if(length(keepX.constraint)>0)
            {
                #remove the variables from keepA.constraint if removed by near.zero.var
                keepX.constraint=match.keepX.constraint(names.remove.X,keepX.constraint)
            }
            #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
            if(any(keepX>ncol(X)))
            {
                ind=which(keepX>ncol(X))
                keepX[ind]=ncol(X)
            }
        }
        
    }else{nzv.A=NULL}
    
    # we need numbers in keepX.constraint from now on
    keepX.constraint= lapply(keepX.constraint,function(x){match(x,colnames(X))})
    keepY.constraint= lapply(keepY.constraint,function(x){match(x,colnames(Y))})
    
    return(list(X=X,Y=Y,ncomp=ncomp,X.names=X.names,Y.names=Y.names,ind.names=ind.names,mode=mode,keepX.constraint=keepX.constraint,
    keepY.constraint=keepY.constraint,keepX=keepX,keepY=keepY,nzv.A=nzv.A))
}

Check.entry.mint.block.spls = function(A, indY, design ,ncomp , scheme , scale ,  bias,
init , tol , verbose,mode , max.iter,study , keepA, keepA.constraint)
{
    
    if (length(A) < 1)
    stop("A is a list with at least 2 blocks.")
    
    if (is.null(indY)) {
        if (mode != "canonical")
        stop("Only canonical deflation can be done with indY empty")
    } else if (!abs(indY - round(indY) < 1e-25)) {
        stop ("indY must be an integer")
    } else if (indY > length(A)) {
        stop ("indY must point to a block of A")
    } else if (!(mode %in% c("canonical", "invariant", "classic", "regression"))) {
        stop("Choose one of the four following modes: canonical, invariant, classic or regression")
    }
    
    x=unlist(lapply(A,nrow))
    if(!isTRUE(all.equal( max(x) ,min(x))))
    stop("The samplesize mmust be the same for all blocks")
    
    
    #set the default study factor
    if(missing(study))
    {
        study=factor(rep(1,nrow(A[[1]])))
    }else{
        study=as.factor(study)
    }
    if(length(study)!=nrow(A[[1]])) stop(paste0("'study' must be a factor of length ",x[1],"."))
    
    #check dimnames, ncomp, keepA and keepA.constraint per block of A
    for(q in 1:length(A))
    {
        check=Check.entry.single(A[[q]], ncomp[q], keepA[[q]], keepA.constraint[[q]],q)
        A[[q]]=check$X
        ncomp[q]=check$ncomp
    }
    
    
    if (!(scheme %in% c("horst", "factorial","centroid"))) {
        stop("Choose one of the three following schemes: horst, centroid or factorial")
    } else {
        if (verbose)
        cat("Computation of the SGCCA block components based on the", scheme, "scheme \n")
    }
    
    if(!init%in%c("svd","svd.single"))
    stop("init should be one of 'svd' or 'svd.single'")
    
    study=as.factor(study)
    
    if(tol<=0)
    stop("tol must be non negative")
    
    if(!is.logical(verbose))
    stop("verbose must be either TRUE or FALSE")
    if(!is.logical(scale))
    stop("scale must be either TRUE or FALSE")
    if(!is.logical(bias))
    stop("bias must be either TRUE or FALSE")
    
    
    return(list(A=A,ncomp=ncomp,study=study))
    
}

# --------------------------------------
# Check.entry.wrapper.sparse.mint.block
# --------------------------------------


Check.entry.wrapper.sparse.mint.block = function(X,
Y,
indY, #only use if Y not provided
ncomp,
keepX, #sparse
keepX.constraint, #hybrid
keepY, #sparse
keepY.constraint, #hybrid
study, #mint
design, #block
init,
scheme, #block
scale,
bias,
near.zero.var,
mode,
tol,
max.iter,
verbose)
{
    #need to give the default values of mint.block.spls to mixOmics
    
    if(!is.list(X))
    {stop("X must be a list")}
    
    if((missing(indY)& missing(Y)) )
    stop("Either 'Y' or 'indY' is needed")
    
    if(missing(ncomp)) {ncomp = rep(1, length(X))}
    
    #check length(ncomp)=length(A)
    if(length(ncomp)!=length(X)) stop("'ncomp' must be a vector of length the number of blocks in X")
    
    #check dimnames and ncomp per block of A
    for(q in 1:length(X))
    {
        check=Check.entry.single(X[[q]], ncomp[q],q=q)
        X[[q]]=check$X
        ncomp[q]=check$ncomp
    }
    
    
    #check ncomp[q]<ncol(X[[q]])
    for(q in 1:length(X))
    {
        ncomp[q] = round(ncomp[q])
        if(ncomp[q] > ncol(X[[q]]))
        {
            warning(paste0("Reset maximum number of variates 'ncomp[",q,"]' to ncol(X[[",q,"]])= ", ncol(X[[q]]), "."))
            ncomp[q] = ncol(X[[q]])
        }
    }
    
    #add names to the blocks if no names or not unique name for each block
    if(length(unique(names(X)))!=length(X))
    names(X)=paste0("block",1:length(X))
    
    
    # construction of keepA and keepA.constraint
    check=check.keepA.and.keepA.constraint(X=X,keepX=keepX,keepX.constraint=keepX.constraint,ncomp=ncomp)
    keepA=check$keepA
    keepA.constraint=check$keepA.constraint
    
    
    # =====================================================
    # mint.block.spls_iteration algo
    # =====================================================
    if(!missing(Y))# Y is not missing, we don't care about indY
    {
        
        if(!missing(indY)){warnings("'Y' and 'indY' are provided, 'Y' is used.")}
        
        if(is.list(Y))
        stop("Y must be a matrix")
        
        Y=as.matrix(Y)
        if (!is.numeric(Y) )
        stop("'Y' must be a numeric matrix.")
        
        #check dimnames and ncomp per block of A
        check=Check.entry.single(Y, max(ncomp),q=1)
        Y=check$X
        
        
        if(missing(keepY.constraint))
        {
            if(missing(keepY)){keepY=rep(ncol(Y),max(ncomp))}
            
            keepA.constraint[[length(X)+1]]=list() #keepY.constraint
        }else{
            if ( is.numeric(unlist(keepY.constraint)) && any(unlist(keepY.constraint) > ncol(Y)))
            stop("each entry of 'keepY.constraint' must be lower or equal to ", ncol(Y), ".")
            
            # check that max(ncomp)>=length(keepY.constraint)
            if(length(keepY.constraint)>max(ncomp))
            stop(paste0("you should have length(keepY.constraint) lower or equal to ",max(ncomp),"."))
            
            if(missing(keepY)) {keepY=rep(ncol(Y),max(ncomp)-length(keepY.constraint))}  #if not missing keepY.constraint but missing keepY, we complete keepY pls-like to have length(keepY.constraint)+length(keepY)=ncomp
            
            #keepA[[length(X)+1]]
        }
        
        #check that keepY is ok
        if(is.list(keepY))
        stop("'keepY' must be a numeric vector")
        if (any(keepY > ncol(Y)))
        stop(paste0("each component of 'keepY' must be lower or equal to ", ncol(Y), "."))
        if(length(keepY)>max(ncomp))
        stop("length of 'keepY' must be lower or equal to ", max(ncomp), ".")
        
        keepA[[length(X)+1]]=keepY
        
        # build the list A and indY
        A=X
        A[[length(A)+1]]=Y
        names(A)[length(A)]="Y"
        indY=length(A)
        
        if (mode == "canonical")
        ncomp = c(ncomp, min(ncomp, ncol(Y) - 1))
        if (mode == "regression")
        ncomp = c(ncomp, max(ncomp))
        #adjust ncomp for Y
        
        
    }else{        #missing(Y) but indY not missing
        A=X #list provided as input
        
    }
    
    # --------------------------------------------------------------------------------
    # at this stage, we have A, indY, keepA, keepA.constraint, ncomp verified
    # --------------------------------------------------------------------------------
    
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
    {stop("Choose one of the four following modes: canonical, invariant, classic or regression")}
    
    #set the default study factor
    if(missing(study))
    {
        study=factor(rep(1,nrow(A[[1]])))
    }else{
        study=as.factor(study)
    }
    if(length(study)!=nrow(A[[1]])) stop(paste0("'study' must be a factor of length ",nrow(A[[1]]),"."))
    
    if(any(table(study)<=1)) stop("At least one study has only one sample, please consider removing before calling the function again")
    if(any(table(study)<5)) warning("At least one study has less than 5 samples, mean centering might not do as expected")
    
    if(missing(init)) init="svd"
    
    if(!init%in%c("svd","svd.single"))
    stop("init should be one of 'svd' or 'svd.single'")
    
    
    
    if (!abs(indY - round(indY) < 1e-25)) {stop ("indY must be an integer")}
    if (indY > length(A)) {stop ("indY must point to a block of A")}
    
    
    
    # =====================================================
    # with or without tau (RGGCA or mint.block.spls algo)
    # =====================================================
    x=unlist(lapply(A,nrow))
    if(!isTRUE(all.equal( max(x) ,min(x))))
    stop("The samplesize must be the same for all blocks")
    
    #check scheme
    if (!(scheme %in% c("horst", "factorial","centroid"))) {
        stop("Choose one of the three following schemes: horst, centroid or factorial")
    } else {
        if (verbose)
        cat("Computation of the SGCCA block components based on the", scheme, "scheme \n")
    }
    
    if(missing(design))
    {design=1 - diag(length(A))}
    
    #check design matrix
    if(nrow(design)!=ncol(design))
    stop(paste0("'design' must be a square matrix."))
    if(nrow(design)!=length(A))
    stop(paste0("'design' must be a square matrix with",length(A),"columns."))
    
    if(!is.numeric(tol) | tol<=0)
    stop("tol must be non negative")
    if(!is.numeric(max.iter) | max.iter<=0)
    stop("max.iter must be non negative")
    
    if(!is.logical(verbose))
    stop("verbose must be either TRUE or FALSE")
    if(!is.logical(scale))
    stop("scale must be either TRUE or FALSE")
    if(!is.logical(bias))
    stop("bias must be either TRUE or FALSE")
    if(!is.logical(near.zero.var))
    stop("near.zero.var must be either TRUE or FALSE")
    
    
    # at this stage keepA.constraint need to be character, to remove easily variables with near zero variance
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(A,nearZeroVar)
        for(q in 1:length(A))
        {
            if (length(nzv.A[[q]]$Position) > 0)
            {
                names.remove.X=colnames(A[[q]])[nzv.A[[q]]$Position]
                A[[q]] = A[[q]][, -nzv.A[[q]]$Position,drop=FALSE]
                if (verbose)
                warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
                if(ncol(A[[q]]) == 0) {stop(paste0("No more variables in",A[[q]]))}
                
                # at this stage, keepA.constraint need to be numbers
                if(length(keepA.constraint[[q]])>0)
                {
                    #remove the variables from keepA.constraint if removed by near.zero.var
                    keepA.constraint[[q]]=match.keepX.constraint(names.remove.X,keepA.constraint[[q]])
                    # replace character by numbers
                    keepA.constraint[[q]]= lapply(keepA.constraint[[q]],function(x){match(x,colnames(A[[q]]))})
                }
                #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
                if(any(keepA[[q]]>ncol(A[[q]])))
                {
                    ind=which(keepA[[q]]>ncol(A[[q]]))
                    keepA[[q]][ind]=ncol(A[[q]])
                }
            }
            
        }
    }else{nzv.A=NULL}
    
    
    
    return(list(A=A,ncomp=ncomp,study=study,keepA=keepA,keepA.constraint=keepA.constraint,
    indY=indY,design=design,init=init,nzv.A=nzv.A))
    
}


# --------------------------------------
# Check.entry.sgcca
# --------------------------------------


Check.entry.sgcca = function(X,
design,
ncomp ,
scheme,
mode,
scale,
init,
bias,
tol,
verbose,
max.iter,
near.zero.var,
keepX,
keepX.constraint)
{
    #need to give the default values of mint.block.spls to mixOmics
    
    if(!is.list(X))
    {stop("X must be a list of at list two matrices")}
    
    if(length(X)<2)
    {stop("X must be a list of at list two matrices")}
    
    if(missing(ncomp)) {ncomp = rep(1, length(X))}
    
    #check dimnames and ncomp per block of A
    for(q in 1:length(X))
    {
        check=Check.entry.single(X[[q]], ncomp[q],q=q)
        X[[q]]=check$X
        ncomp[q]=check$ncomp
    }
    
    
    #check length(ncomp)=length(A)
    if(length(ncomp)!=length(X)) stop("'ncomp' must be a vector of length the number of blocks in X")
    #check ncomp[q]<ncol(X[[q]])
    for(q in 1:length(X))
    {
        ncomp[q] = round(ncomp[q])
        if(ncomp[q] > ncol(X[[q]]))
        {
            warning(paste0("Reset maximum number of variates 'ncomp[",q,"]' to ncol(X[[",q,"]])= ", ncol(X[[q]]), "."))
            ncomp[q] = ncol(X[[q]])
        }
    }
    
    #add names to the blocks if no names or not unique name for each block
    if(length(unique(names(X)))!=length(X))
    names(X)=paste0("block",1:length(X))
    
    A=X#input
    
    if(missing(init)) init="svd"
    
    if(!init%in%c("svd","svd.single"))
    stop("init should be one of 'svd' or 'svd.single'")
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
    {stop("Choose one of the four following modes: canonical, invariant, classic or regression")}
    
    
    # =====================================================
    # with or without tau (RGGCA or mint.block.spls algo)
    # =====================================================
    
    x=unlist(lapply(A,nrow))
    if(!isTRUE(all.equal( max(x) ,min(x))))
    stop("The samplesize must be the same for all blocks")
    
    
    
    #check scheme
    if(missing(scheme)) scheme= "centroid"
    if (!(scheme %in% c("horst", "factorial","centroid"))) {
        stop("Choose one of the three following schemes: horst, centroid or factorial")
    } else {
        if (verbose)
        cat("Computation of the SGCCA block components based on the", scheme, "scheme \n")
    }
    
    
    if(missing(design))
    {design=1 - diag(length(A))}
    
    #check design matrix
    if(nrow(design)!=ncol(design))
    stop(paste0("'design' must be a square matrix."))
    if(nrow(design)!=length(A))
    stop(paste0("'design' must be a square matrix with",length(A),"columns."))
    
    
    if(missing(bias)) bias= FALSE
    if(missing(verbose)) verbose= FALSE
    
    if(tol<=0)
    stop("tol must be non negative")
    
    if(max.iter<=0)
    stop("max.iter must be non negative")
    
    if(!is.logical(verbose))
    stop("verbose must be either TRUE or FALSE")
    if(!is.logical(scale))
    stop("scale must be either TRUE or FALSE")
    if(!is.logical(bias))
    stop("bias must be either TRUE or FALSE")
    if(!is.logical(near.zero.var))
    stop("near.zero.var must be either TRUE or FALSE")
    
    
    # construction of keepA and keepA.constraint
    check=check.keepA.and.keepA.constraint(X=A,keepX=keepX,keepX.constraint=keepX.constraint,ncomp=ncomp)
    keepA=check$keepA
    keepA.constraint=check$keepA.constraint
    
    
    
    # at this stage keepA.constraint need to be character, to remove easily variables with near zero variance
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(A,nearZeroVar)
        for(q in 1:length(A))
        {
            if (length(nzv.A[[q]]$Position) > 0)
            {
                names.remove.X=colnames(A[[q]])[nzv.A[[q]]$Position]
                A[[q]] = A[[q]][, -nzv.A[[q]]$Position,drop=FALSE]
                if (verbose)
                warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
                if(ncol(A[[q]]) == 0) {stop(paste0("No more variables in",A[[q]]))}
                
                # at this stage, keepA.constraint need to be numbers
                if(length(keepA.constraint[[q]])>0)
                {
                    #remove the variables from keepA.constraint if removed by near.zero.var
                    keepA.constraint[[q]]=match.keepX.constraint(names.remove.X,keepA.constraint[[q]])
                    # replace character by numbers
                    keepA.constraint[[q]]= lapply(keepA.constraint[[q]],function(x){match(x,colnames(A[[q]]))})
                }
                #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
                if(any(keepA[[q]]>ncol(A[[q]])))
                {
                    ind=which(keepA[[q]]>ncol(A[[q]]))
                    keepA[[q]][ind]=ncol(A[[q]])
                }
            }
            
        }
    }else{nzv.A=NULL}
    
    
    return(list(A=A,ncomp=ncomp,design=design,init=init,scheme=scheme,verbose=verbose,bias=bias,nzv.A=nzv.A,
    keepA=keepA,keepA.constraint=keepA.constraint))
    
}



# --------------------------------------
# Check.entry.rgcca
# --------------------------------------


Check.entry.rgcca = function(X,
design ,
tau ,
ncomp,
scheme,
mode,
scale,
init,
bias,
tol,
verbose,
max.iter,
near.zero.var)
{
    #need to give the default values of mint.block.spls to mixOmics
    
    if(!is.list(X))
    {stop("X must be a list of at list two matrices")}
    
    if(length(X)<2)
    {stop("X must be a list of at list two matrices")}
    
    if(is.null(tau))
    stop("'tau' is needed")
    
    if(missing(ncomp)) {ncomp = rep(1, length(X))}
    
    #check dimnames and ncomp per block of A
    for(q in 1:length(X))
    {
        check=Check.entry.single(X[[q]], ncomp[q],q=q)
        X[[q]]=check$X
        ncomp[q]=check$ncomp
    }
    
    
    #check length(ncomp)=length(A)
    if(length(ncomp)!=length(X)) stop("'ncomp' must be a vector of length the number of blocks in X")
    #check ncomp[q]<ncol(X[[q]])
    for(q in 1:length(X))
    {
        ncomp[q] = round(ncomp[q])
        if(ncomp[q] > ncol(X[[q]]))
        {
            warning(paste0("Reset maximum number of variates 'ncomp[",q,"]' to ncol(X[[",q,"]])= ", ncol(X[[q]]), "."))
            ncomp[q] = ncol(X[[q]])
        }
    }
    
    #add names to the blocks if no names or not unique name for each block
    if(length(unique(names(X)))!=length(X))
    names(X)=paste0("block",1:length(X))
    
    A=X#input
    
    
    
    if(is.numeric(tau))
    {
        if(any(tau<0) | any(tau>1)) stop("'tau' contains either values between 0 and 1, or 'optimal'.")
        if(is.vector(tau))
        {
            if(length(tau)!=length(A)) stop(paste0("'tau' must be of length ",length(A),"."))
            tau = matrix(rep(tau, max(ncomp)), nrow = max(ncomp), ncol = length(tau), byrow = TRUE)
        }
    }else{
        if(tau!="optimal") stop("'tau' contains either values between 0 and 1, or 'optimal'.")
    }
    
    if(missing(init)) init="svd.single"
    
    if (init != "svd.single") stop("init should be 'svd.single'.")
    if(missing(mode)) mode="canonical"
    if (mode != "canonical") stop("Only canonical deflation can be done when 'tau' is provided. Try again with mode='canonical'")
    
    
    
    # =====================================================
    # with or without tau (RGGCA or mint.block.spls algo)
    # =====================================================
    
    x=unlist(lapply(A,nrow))
    if(!isTRUE(all.equal( max(x) ,min(x))))
    stop("The samplesize must be the same for all blocks")
    
    
    
    #check scheme
    if(missing(scheme)) scheme= "centroid"
    if (!(scheme %in% c("horst", "factorial","centroid"))) {
        stop("Choose one of the three following schemes: horst, centroid or factorial")
    } else {
        if (verbose)
        cat("Computation of the SGCCA block components based on the", scheme, "scheme \n")
    }
    
    
    if(missing(design))
    {design=1 - diag(length(A))}
    
    #check design matrix
    if(nrow(design)!=ncol(design))
    stop(paste0("'design' must be a square matrix."))
    if(nrow(design)!=length(A))
    stop(paste0("'design' must be a square matrix with",length(A),"columns."))
    
    
    if(missing(bias)) bias= FALSE
    if(missing(verbose)) verbose= FALSE
    if(missing(near.zero.var)) near.zero.var= FALSE
    
    if(tol<=0)
    stop("tol must be non negative")
    
    if(max.iter<=0)
    stop("max.iter must be non negative")
    
    if(!is.logical(verbose))
    stop("verbose must be either TRUE or FALSE")
    if(!is.logical(scale))
    stop("scale must be either TRUE or FALSE")
    if(!is.logical(bias))
    stop("bias must be either TRUE or FALSE")
    if(!is.logical(near.zero.var))
    stop("near.zero.var must be either TRUE or FALSE")
    
    # at this stage keepA.constraint need to be character, to remove easily variables with near zero variance
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(A,nearZeroVar)
        for(q in 1:length(A))
        {
            if (length(nzv.A[[q]]$Position) > 0)
            {
                names.remove.X=colnames(A[[q]])[nzv.A[[q]]$Position]
                A[[q]] = A[[q]][, -nzv.A[[q]]$Position,drop=FALSE]
                if (verbose)
                warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
                if(ncol(A[[q]]) == 0) {stop(paste0("No more variables in",A[[q]]))}
                
            }
            
        }
    }else{nzv.A=NULL}
    
    
    return(list(A=A,ncomp=ncomp,design=design,init=init,scheme=scheme,verbose=verbose,bias=bias,nzv.A=nzv.A))
    
}



# --------------------------------------
# check keepA and keepA.constraint
# --------------------------------------

check.keepA.and.keepA.constraint=function(X,keepX,keepX.constraint,ncomp)
{
    # X:data
    # keepA
    # keepA.constraint
    
    keepA=list()
    keepA.constraint=list()
    if(missing(keepX.constraint) || length(keepX.constraint)==0)
    {
        if(missing(keepX) || length(keepX)==0)
        {
            #if both keepX.constraint and keepX are missing, pls-like: keepX=ncol(X)
            for(q in 1:length(X)) {keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))} #keepX
            names(keepA)=names(X)
        }else{
            
            if( !is.list(keepX))
            stop("'keepX' must be a list")
            
            for(q in 1:length(X))
            {
                
                if(q<=length(keepX))
                {
                    #checking entries of keepX
                    if(is.list(keepX[[q]]))
                    stop(paste0("keepX[[",q,"]]' must be a vector"))
                    if (any(keepX[[q]] > ncol(X[[q]])))
                    stop(paste0("each component of 'keepX[[",q,"]]' must be lower or equal to ncol(X[[",q,"]])=",ncol(X[[q]]),"."))
                    if (any(keepX[[q]]<0))
                    stop(paste0("each component of 'keepX[[",q,"]]' must be non negative."))
                    if(length(keepX[[q]])>ncomp[q])
                    stop(paste0("length of 'keepX[[",q,"]]' must be lower or equal to ncomp[",q,"]=",ncomp[q], "."))
                    
                    keepA[[q]]=keepX[[q]]
                    if(length(keepA[[q]])<max(ncomp)) {keepA[[q]]=c(keepA[[q]],rep(ncol(X[[q]]),max(ncomp)-length(keepA[[q]])))} #complete the keepX already provided
                }else{
                    keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))##
                }
            }
            
        }
        
        for(q in 1:length(X)) {keepA.constraint[[q]]=list()} #keepX.constraint
        #print(keepA)
        
    }else{
        
        
        #check entries keepX.constraint
        if(length(keepX.constraint)>length(X))
        stop("length(keepX.constraint) is higher than the number of blocks in X")
        
        if(length(keepX.constraint)<length(X))
        for(q in (length(keepX.constraint)+1):length(X))
        keepX.constraint[[q]]=list() #artificially creating an entry tht we won't use, so we can access all entries from 1 to length(X)
        
        #check names of keepX.constraint, gives the name of the blocks
        if(length(unique(names(keepX.constraint)))!=length(keepX.constraint) | sum(is.na(match(names(keepX.constraint),names(X)))))
        {names(keepX.constraint)=names(X)[1:length(keepX.constraint)]}
        
        
        # check that ncomp>=length(keepX.constraint)
        for(q in 1:length(X))
        {
            if(length(keepX.constraint[[q]])>ncomp[q])
            stop(paste0("you should have length(keepX.constraint[[",q,"]]) lower or equal to ncomp[",q,"]=",ncomp[q],"."))
            #if(!is.list(keepX.constraint[[q]]))
            #stop(paste0("'keepX.constraint[[",q,"]]' must be a list of length ",ncomp[q]," of variables to keep on each of the ncomp[",q,"]=",ncomp[q]," components for block ",q,"."))
        }
        
        if(missing(keepX) || length(keepX)==0)
        {
            #if not missing keepX.constraint but missing keepX, we complete keepX pls-like to have length(keepX.constraint)+length(keepX)=ncomp
            for(q in 1:length(X))
            {
                keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepX.constraint[[q]]))
                
            }
            names(keepA)=names(X)
        }else{ #complete keepA so that length(keepX.constraint[[q]])+length(keepX[[q]])=max(ncomp)
            
            if( !is.list(keepX))
            stop("'keepX' must be a list")
            
            for(q in 1:length(X))
            {
                
                #check the entries provided before completed by pls-like
                if(q<=length(keepX)& q<=length(keepX.constraint))
                if((length(keepX.constraint[[q]])+length(keepX[[q]]))>ncomp[q]) stop(paste0("length (keepX.constraint[[",q,"]]) + length(keepX[[",q,"]]) = ",(length(keepX.constraint[[q]])+length(keepX[[q]])),"; it should be lower or equal to ncomp[",q,"]=",ncomp[q], "."))
                if(q<=length(keepX))
                {
                    #checking entries of keepX
                    if(is.list(keepX[[q]]))
                    stop(paste0("keepX[[",q,"]]' must be a vector"))
                    if (any(keepX[[q]] > ncol(X[[q]])))
                    stop(paste0("each component of 'keepX[[",q,"]]' must be lower or equal to ncol(X[[",q,"]])=",ncol(X[[q]]),"."))
                    if (any(keepX[[q]]<0))
                    stop(paste0("each component of 'keepX[[",q,"]]' must be non negative."))
                    if(length(keepX[[q]])>ncomp[q])
                    stop(paste0("length of 'keepX[[",q,"]]' must be lower or equal to ncomp[",q,"]=",ncomp[q], "."))
                    
                    keepA[[q]]=keepX[[q]]
                    if(length(keepA[[q]])<max(ncomp)) {keepA[[q]]=c(keepA[[q]],rep(ncol(X[[q]]),max(ncomp)-length(keepA[[q]])-length(keepX.constraint[[q]])))} #complete the keepX already provided
                    
                }else{
                    keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepX.constraint[[q]]))
                }
            }
            
        }
        
        
        # match keepX.constraint and the colnames of X in order for keepX.constraint to be a list of character
        # safety if keepX.constraint contains a mixed of character/numeric. It should one or the other, not a mix
        for(q in 1:length(X))
        {
            if(length(keepX.constraint[[q]])>0)
            {
                if (is.numeric(unlist(keepX.constraint[[q]])) && any(unlist(keepX.constraint[[q]]) > ncol(X[[q]])))
                stop(paste0("each entry of 'keepX.constraint",q,"' must be lower or equal than ", ncol(X[[q]]), "."))
                
                if(!is.numeric(unlist(keepX.constraint[[q]])))
                {
                    ind=match(unlist(keepX.constraint[[q]]),colnames(X[[q]]))
                    if(sum(is.na(ind))>0) stop("'keepX.constraint' must contains a subset of colnames(X) or the position of the X-variables you wish to keep.")
                }
                X.indice=X[[q]][,unlist(keepX.constraint[[q]]),drop=FALSE]
                keepX.constraint[[q]]=relist(colnames(X.indice),skeleton=keepX.constraint[[q]])
            }
            
            # we need numbers in keepX.constraint from now on
            keepX.constraint[[q]]= sapply(keepX.constraint[[q]],function(x){match(x,colnames(X[[q]]))})
        }
        
        keepA.constraint=keepX.constraint
        
    }
    #print("constraint")
    #print(keepA.constraint)
    #print("keepA")
    #print(keepA)
    if(FALSE)
    {
        #check that keepX is ok
        for(q in 1:length(X))
        {
            if(is.list(keepA[[q]]))
            stop(paste0("keepX[[",q,"]]' must be a vector"))
            if (any(keepA[[q]] > ncol(X[[q]])))
            stop(paste0("each component of 'keepX[[",q,"]]' must be lower or equal to ncol(X[[",q,"]])=",ncol(X[[q]]),"."))
            if(length(keepA[[q]])>max(ncomp))
            stop(paste0("length of 'keepX[[",q,"]]' must be lower or equal to ncomp[",q,"]=",ncomp[q], "."))
        }
        
        
        #check that length(keepA.constraint)+length(keepA)=ncomp for all blocks
        for(q in 1:length(X))
        {
            if (any(keepA[[q]] > ncol(X[[q]])))
            stop(paste0("each component of 'keepA[[",q,"]]' must be lower or equal to ", ncol(X[[q]]), "."))
            #check keepX and keepX.constraint
            if((length(keepA.constraint[[q]])+length(keepA[[q]]))!=max(ncomp)) stop(paste0("length (keepX.constraint[[",q,"]]) + length(keepX[[",q,"]])=",(length(keepA.constraint[[q]])+length(keepA[[q]])),"; it should be max(ncomp)."))
        }
    }
    #print(keepA)
    #print("constr")
    #print(keepA.constraint)
    
    
    
    
    names(keepA)=names(X)
    
    return(list(keepA=keepA,keepA.constraint=keepA.constraint))
}


