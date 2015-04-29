# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created: 22-04-2015
# last modification: 22-04-2015
# Copyright (C) 2014
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



# --------------------------------------
# Check.entry.single
# --------------------------------------
Check.entry.single = function(X,  ncomp, keepX, keepX.constraint,q)
{
    
    
    #  if(length(levels(study)) == 1)  # Aida
    #  stop("\nstudys must have more than one level")      #WHY?
    
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
    
    # if (any(keepX > ncol(X)))
    #stop(paste0("each component of 'keepA[[",q,"]]' must be lower or equal to ", P, "."))
    ##check keepX and keepX.constraint
    #if((length(keepX.constraint)+length(keepX))!=ncomp) stop(paste0("length (keepA.constraint[[",q,"]]) + length(keepA[[",q,"]]) should be ncomp"))
    #-- add colnames/rownames if missing --#
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
# Check.entry.meta.block.spls
# --------------------------------------

Check.entry.meta.block.spls = function(A, indY, design ,ncomp , scheme , scale ,  bias,
init , tol , verbose,mode, sparse , max.iter,study , keepA, keepA.constraint)
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
    if(!is.logical(sparse))
    stop("sparse must be either TRUE or FALSE")
    
    
    return(list(A=A,ncomp=ncomp,study=study))
    
}



# --------------------------------------
# Check.entry.pls
# --------------------------------------

Check.entry.pls = function(X, Y, ncomp, keepX, keepY, keepX.constraint,keepY.constraint)
{
    
    
    #  if(length(levels(study)) == 1)  # Aida
    #  stop("\nstudys must have more than one level")      #WHY?
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")
    
    X = as.matrix(X)
    
    Y = as.matrix(Y)
    
    if (!is.numeric(X) || !is.numeric(Y))
    stop("'X' and/or 'Y' must be a numeric matrix.")
    
    N = nrow(X)
    Q = ncol(Y)
    P= ncol(X)
    
    if ((N != nrow(Y)))
    stop("unequal number of rows in 'X' and 'Y'.")
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    
    ncomp = round(ncomp)
    if(ncomp > P)
    {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", P, ".")
        ncomp = P
    }
    
    
    if(missing(keepX)) # so the function can be used for pls as well
    {keepX=rep(P,ncomp)}
    if(missing(keepY)) # so the function can be used for pls as well
    {keepY=rep(Q,ncomp)}
    
    if(missing(keepX.constraint))
    keepX.constraint=list()
    if(missing(keepY.constraint))
    keepY.constraint=list()
    
    if((length(keepX.constraint)+length(keepX))!=ncomp) stop(paste0("length (keepX.constraint) + length(keepX) should be ncomp"))
    if((length(keepY.constraint)+length(keepY))!=ncomp) stop(paste0("length (keepY.constraint) + length(keepY) should be ncomp"))
    
    
    if (any(keepX > ncol(X)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    
    if (length(keepY) != ncomp)
    stop("length of 'keepX' must be equal to ", ncomp, ".")
    if (any(keepY > ncol(Y)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    
    
    
    
    #-- initialisation des matrices --#
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
    
    return(list(X=X,Y=Y,ncomp=ncomp,X.names=X.names,Y.names=Y.names,ind.names=ind.names))
}

# --------------------------------------
# Check.entry.pls.single
# --------------------------------------

Check.entry.pls.single = function(X, ncomp, keepX,keepX.constraint)
{
    
    #  if(length(levels(study)) == 1)  # Aida
    #  stop("\nstudys must have more than one level")      #WHY?
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")
    
    X = as.matrix(X)
    
    
    if (!is.numeric(X) )
    stop("'X' must be a numeric matrix.")
    
    N = nrow(X)
    P= ncol(X)
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    ncomp = round(ncomp)
    if(ncomp > P)
    {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", P, ".")
        ncomp = P
    }
    
    
    if(missing(keepX)) # so the function can be used for pls as well
    {keepX=rep(P,ncomp)}
    
    if(missing(keepX.constraint))
    keepX.constraint=list()
    
    if((length(keepX.constraint)+length(keepX))!=ncomp) stop(paste0("length (keepX.constraint) + length(keepX) should be ncomp"))
    
    
    if (any(keepX > ncol(X)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    
    
    
    #-- add colnames/rownames if missing --#
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
    
    
    if(length(unique(X.names))!=P) stop("Unique indentifier is needed for the columns of X")
    
    return(list(X=X,ncomp=ncomp,X.names=X.names,ind.names=ind.names))
}




# --------------------------------------
# Check.entry.mixOmics.list.save
# --------------------------------------


Check.entry.mixOmics.list.save = function(X,
Y,
indY, #only use if Y not provided
ncomp,
keepX, #sparse
keepX.constraint, #hybrid
keepY, #sparse
keepY.constraint, #hybrid
study, #meta
design, #block
tau,# rgcca, number between 0,1 or "optimal"
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
    #need to give the default values of meta.block.spls to mixOmics
    
    if((missing(indY)& missing(Y)) & is.null(tau))
    stop("Either 'Y', 'indY' or 'tau' is needed")
    
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
    
    
    # construction of keepA and keepA.constraint
    keepA=list()
    keepA.constraint=list()
    if(missing(keepX.constraint))
    {
        if(missing(keepX))
        {
            #if both keepX.constraint and keepX are missing, pls-like: keepX=ncol(X)
            for(q in 1:length(X)) {keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))} #keepX
            names(keepA)=names(X)
        }else{
            
            for(q in 1:length(X))
            {
                if(q<=length(keepX))
                {keepA[[q]]=keepX[[q]]
                }else{
                    keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))##
                }
            }
            
        }
        
        for(q in 1:length(X)) {keepA.constraint[[q]]=list()} #keepX.constraint
        #print(keepA)
        
    }else{
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
            stop(paste0("you should have length(keepX.constraint[[",q,"]]) lower or equal to ncomp[",q,"]."))
        }
        
        if(missing(keepX))
        {
            #if not missing keepX.constraint but missing keepX, we complete keepX pls-like to have length(keepX.constraint)+length(keepX)=ncomp
            for(q in 1:length(A))
            {
                keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepX.constraint[[q]]))
                
            }
            names(keepA)=names(X)
        }else{ #complete keepA so that length(keepX.constraint[[q]])+length(keepX[[q]])=max(ncomp)
            
            for(q in 1:length(X))
            {
                if(q<=length(keepX))
                {keepA[[q]]=keepX[[q]]
                }else{
                    keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepX.constraint[[q]]))
                }
            }
            
        }
        
        keepA.constraint=keepX.constraint
        
    }
    # print("constraint")
    #print(keepA.constraint)
    #print("keepA")
    #print(keepA)
    
    #check that keepX is ok
    for(q in 1:length(X))
    {
        if(is.list(keepA[[q]]))
        stop(paste0("keepX[[",q,"]]' must be a vector"))
        if (any(keepA[[q]] > ncol(X[[q]])))
        stop(paste0("each component of 'keepX[[",q,"]]' must be lower or equal to ", ncol(X[[q]]), "."))
        if(length(keepA[[q]])>max(ncomp))
        stop(paste0("length of 'keepX[[",q,"]]' must be lower or equal to ", ncomp[q], "."))
    }
    
    
    
    if(is.null(tau))
    {
        
        # =====================================================
        # meta.block.spls_iteration algo
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
                # check that max(ncomp)>=length(keepY.constraint)
                if(length(keepY.constraint)>max(ncomp))
                stop(paste0("you should have length(keepY.constraint) lower or equal to ",max(ncomp),"."))
                
                if(missing(keepY)) {keepY=rep(ncol(Y),max(ncomp)-length(keepY.constraint))}  #if not missing keepY.constraint but missing keepY, we complete keepY pls-like to have length(keepY.constraint)+length(keepY)=ncomp
                
                #keepA[[length(X)+1]]
            }
            
            #check that keepY is ok
            if (any(keepY > ncol(Y)))
            stop(paste0("each component of 'keepY' must be lower or equal to ", ncol(Y), "."))
            if(length(keepY)>max(ncomp))
            stop("length of 'keepX' must be lower or equal to ", max(ncomp), ".")
            
            keepA[[length(X)+1]]=keepY
            
            # build the list A and indY
            A=X
            A$Y=Y
            indY=length(A)
            ncomp=c(ncomp,max(ncomp)) #adjust ncomp for Y
            
            
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
        
        
        if(missing(init)) init="svd"
        
        if(!init%in%c("svd","svd.single"))
        stop("init should be one of 'svd' or 'svd.single'")
        
        
        
        if (!abs(indY - round(indY) < 1e-25)) {stop ("indY must be an integer")}
        if (indY > length(A)) {stop ("indY must point to a block of A")}
        
        
    }else{
        # =====================================================
        #  RGCCA algo
        # =====================================================
        A=X #list provided as input
        indY=NULL #no indY for RGCCA
        
        
        if(!missing(study) && nlevels(study)!=1)
        stop("'study' cannot be used when 'tau' is provided.")
        
        study=factor(rep(1,nrow(A[[1]])))
        
        
        if(is.vector(tau))
        {
            if(length(tau)!=length(A)) stop(paste0("'tau' must be of length ",length(A),"."))
            tau = matrix(rep(tau, max(ncomp)), nrow = max(ncomp), ncol = length(tau), byrow = TRUE)
        }
        
        
        if(is.numeric(tau))
        {
            if(any(tau<0) | any(tau>1)) stop("'tau' contains either values between 0 and 1, or 'optimal'.")
            
        }else{
            if(tau!="optimal") stop("'tau' contains either values between 0 and 1, or 'optimal'.")
        }
        
        if(missing(init)) init="svd.single"
        
        if (init != "svd.single") stop("init should be 'svd.single'.")
        if (mode != "canonical") stop("Only canonical deflation can be done when 'tau' is provided. Try again with mode='canonical'")
        
    }
    
    # =====================================================
    # with or without tau (RGGCA or meta.block.spls algo)
    # =====================================================
    
    #check that length(keepA.constraint)+length(keepA)=ncomp for all blocks
    for(q in 1:length(A))
    {
        if (any(keepA[[q]] > ncol(A[[q]])))
        stop(paste0("each component of 'keepA[[",q,"]]' must be lower or equal to ", ncol(A[[q]]), "."))
        #check keepX and keepX.constraint
        if((length(keepA.constraint[[q]])+length(keepA[[q]]))!=max(ncomp)) stop(paste0("length (keepA.constraint[[",q,"]]) + length(keepA[[",q,"]]) should be max(ncomp)."))
    }
    
    
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
    
    return(list(A=A,ncomp=ncomp,study=study,keepA=keepA,keepA.constraint=keepA.constraint,
    indY=indY,design=design,init=init))
    
}



# --------------------------------------
# check keepA and keepA.constraint
# --------------------------------------

check.keepA.and.keepA.constraint=function(X,keepA,keepA.constraint,ncomp)
{
    # X:data
    # keepA
    # keepA.constraint
    
    if(missing(keepA.constraint))
    {
        if(missing(keepA))
        {
            keepA=list()
            for(q in 1:length(X)) {keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))}
            names(keepA)=names(X)
        }else{
            for(q in 1:length(X))
            {
                if(q>length(keepA))
                {keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))}
            }
            
        }
        for(q in 1:length(X)) {keepA.constraint[[q]]=list()} #keepA.constraint
    }else{
        # check that max(ncomp)>=length(keepA.constraint)
        if(length(keepA.constraint)>max(ncomp))
        stop(paste0("you should have length(keepA.constraint) lower or equal to ",max(ncomp),"."))
        
        if(length(keepA.constraint)>length(X))
        stop("length(keepA.constraint) is higher than the number of blocks in 'X'")
        
        if(length(keepA.constraint)<length(X))
        for(q in (length(keepA.constraint)+1):length(X))
        keepA.constraint[[q]]=list() #artificially creating an entry tht we won't use, so we can access all entries from 1 to length(X)
        
        #check names of keepA.constraint, gives the name of the blocks
        if(length(unique(names(keepA.constraint)))!=length(keepA.constraint) | sum(is.na(match(names(keepA.constraint),names(X)))))
        {names(keepA.constraint)=names(X)[1:length(keepA.constraint)]}
        
        
        # check that ncomp>=length(keepA.constraint)
        for(q in 1:length(X))
        {
            if(length(keepA.constraint[[q]])>ncomp[q])
            stop(paste0("you should have length(keepA.constraint[[",q,"]]) lower or equal to ncomp[",q,"]."))
        }
        
        if(missing(keepA))
        {
            #if not missing keepA.constraint but missing keepA, we complete keepA pls-like to have length(keepA.constraint)+length(keepA)=ncomp
            for(q in 1:length(X))
            {
                keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepA.constraint[[q]]))
                
            }
            names(keepA)=names(X)
        }else{ #complete keepA so that length(keepA.constraint[[q]])+length(keepA[[q]])=max(ncomp)
            
            for(q in 1:length(X))
            {
                if(q>length(keepA))
                {keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepA.constraint[[q]]))}
            }
            
        }
    }
    return(list(keepA=keepA,keepA.constraint=keepA.constraint))
}




# --------------------------------------
# Check.entry.wrapper.sparse.meta.block
# --------------------------------------


Check.entry.wrapper.sparse.meta.block = function(X,
Y,
indY, #only use if Y not provided
ncomp,
keepX, #sparse
keepX.constraint, #hybrid
keepY, #sparse
keepY.constraint, #hybrid
study, #meta
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
    #need to give the default values of meta.block.spls to mixOmics
    
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
    keepA=list()
    keepA.constraint=list()
    if(missing(keepX.constraint))
    {
        if(missing(keepX))
        {
            #if both keepX.constraint and keepX are missing, pls-like: keepX=ncol(X)
            for(q in 1:length(X)) {keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))} #keepX
            names(keepA)=names(X)
        }else{
            
            for(q in 1:length(X))
            {
                if(q<=length(keepX))
                {keepA[[q]]=keepX[[q]]
                    if(length(keepA[[q]])<max(ncomp)) {keepA[[q]]=c(keepA[[q]],rep(ncol(X[[q]]),max(ncomp)-length(keepA[[q]])))} #complete the keepX already provided
                }else{
                    keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))##
                }
            }
            
        }
        
        for(q in 1:length(X)) {keepA.constraint[[q]]=list()} #keepX.constraint
        #print(keepA)
        
    }else{
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
            stop(paste0("you should have length(keepX.constraint[[",q,"]]) lower or equal to ncomp[",q,"]."))
        }
        
        if(missing(keepX))
        {
            #if not missing keepX.constraint but missing keepX, we complete keepX pls-like to have length(keepX.constraint)+length(keepX)=ncomp
            for(q in 1:length(X))
            {
                keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepX.constraint[[q]]))
                
            }
            names(keepA)=names(X)
        }else{ #complete keepA so that length(keepX.constraint[[q]])+length(keepX[[q]])=max(ncomp)
            
            for(q in 1:length(X))
            {
                if(q<=length(keepX))
                {keepA[[q]]=keepX[[q]]
                    if(length(keepA[[q]])<max(ncomp)) {keepA[[q]]=c(keepA[[q]],rep(ncol(X[[q]]),max(ncomp)-length(keepA[[q]])))} #complete the keepX already provided
                }else{
                    keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepX.constraint[[q]]))
                }
            }
            
        }
        
        keepA.constraint=keepX.constraint
        
    }
    # print("constraint")
    #print(keepA.constraint)
    #print("keepA")
    #print(keepA)
    
    #check that keepX is ok
    for(q in 1:length(X))
    {
        if(is.list(keepA[[q]]))
        stop(paste0("keepX[[",q,"]]' must be a vector"))
        if (any(keepA[[q]] > ncol(X[[q]])))
        stop(paste0("each component of 'keepX[[",q,"]]' must be lower or equal to ", ncol(X[[q]]), "."))
        if(length(keepA[[q]])>max(ncomp))
        stop(paste0("length of 'keepX[[",q,"]]' must be lower or equal to ", ncomp[q], "."))
    }
    
    
    
    # =====================================================
    # meta.block.spls_iteration algo
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
            # check that max(ncomp)>=length(keepY.constraint)
            if(length(keepY.constraint)>max(ncomp))
            stop(paste0("you should have length(keepY.constraint) lower or equal to ",max(ncomp),"."))
            
            if(missing(keepY)) {keepY=rep(ncol(Y),max(ncomp)-length(keepY.constraint))}  #if not missing keepY.constraint but missing keepY, we complete keepY pls-like to have length(keepY.constraint)+length(keepY)=ncomp
            
            #keepA[[length(X)+1]]
        }
        
        #check that keepY is ok
        if (any(keepY > ncol(Y)))
        stop(paste0("each component of 'keepY' must be lower or equal to ", ncol(Y), "."))
        if(length(keepY)>max(ncomp))
        stop("length of 'keepX' must be lower or equal to ", max(ncomp), ".")
        
        keepA[[length(X)+1]]=keepY
        
        # build the list A and indY
        A=X
        A[[length(A)+1]]=Y
        names(A)[length(A)]="Y"
        indY=length(A)
        ncomp=c(ncomp,max(ncomp)) #adjust ncomp for Y
        
        
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
    
    
    if(missing(init)) init="svd"
    
    if(!init%in%c("svd","svd.single"))
    stop("init should be one of 'svd' or 'svd.single'")
    
    
    
    if (!abs(indY - round(indY) < 1e-25)) {stop ("indY must be an integer")}
    if (indY > length(A)) {stop ("indY must point to a block of A")}
    
    
    
    # =====================================================
    # with or without tau (RGGCA or meta.block.spls algo)
    # =====================================================
    
    #check that length(keepA.constraint)+length(keepA)=ncomp for all blocks
    for(q in 1:length(A))
    {
        if (any(keepA[[q]] > ncol(A[[q]])))
        stop(paste0("each component of 'keepA[[",q,"]]' must be lower or equal to ", ncol(A[[q]]), "."))
        #check keepX and keepX.constraint
        if((length(keepA.constraint[[q]])+length(keepA[[q]]))!=max(ncomp)) stop(paste0("length (keepA.constraint[[",q,"]]) + length(keepA[[",q,"]]) should be max(ncomp)."))
    }
    
    
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
    
    return(list(A=A,ncomp=ncomp,study=study,keepA=keepA,keepA.constraint=keepA.constraint,
    indY=indY,design=design,init=init))
    
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
    #need to give the default values of meta.block.spls to mixOmics
    
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
        
    if(is.vector(tau))
    {
        if(length(tau)!=length(A)) stop(paste0("'tau' must be of length ",length(A),"."))
        tau = matrix(rep(tau, max(ncomp)), nrow = max(ncomp), ncol = length(tau), byrow = TRUE)
    }
    
    
    if(is.numeric(tau))
    {
        if(any(tau<0) | any(tau>1)) stop("'tau' contains either values between 0 and 1, or 'optimal'.")
        
    }else{
        if(tau!="optimal") stop("'tau' contains either values between 0 and 1, or 'optimal'.")
    }
    
    if(missing(init)) init="svd.single"
    
    if (init != "svd.single") stop("init should be 'svd.single'.")
    if (mode != "canonical") stop("Only canonical deflation can be done when 'tau' is provided. Try again with mode='canonical'")
    
    
    
    # =====================================================
    # with or without tau (RGGCA or meta.block.spls algo)
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
    
    return(list(A=A,ncomp=ncomp,design=design,init=init,scheme=scheme,verbose=verbose,bias=bias,near.zero.var=near.zero.var))
    
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
if(missing(keepX.constraint))
{
    if(missing(keepX))
    {
        #if both keepX.constraint and keepX are missing, pls-like: keepX=ncol(X)
        for(q in 1:length(X)) {keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))} #keepX
        names(keepA)=names(X)
    }else{
        
        for(q in 1:length(X))
        {
            if(q<=length(keepX))
            {keepA[[q]]=keepX[[q]]
                if(length(keepA[[q]])<max(ncomp)) {keepA[[q]]=c(keepA[[q]],rep(ncol(X[[q]]),max(ncomp)-length(keepA[[q]])))} #complete the keepX already provided
            }else{
                keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))##
            }
        }
        
    }
    
    for(q in 1:length(X)) {keepA.constraint[[q]]=list()} #keepX.constraint
    #print(keepA)
    
}else{
    


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
        stop(paste0("you should have length(keepX.constraint[[",q,"]]) lower or equal to ncomp[",q,"]."))
    }
    
    if(missing(keepX))
    {
        #if not missing keepX.constraint but missing keepX, we complete keepX pls-like to have length(keepX.constraint)+length(keepX)=ncomp
        for(q in 1:length(X))
        {
            keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepX.constraint[[q]]))
            
        }
        names(keepA)=names(X)
    }else{ #complete keepA so that length(keepX.constraint[[q]])+length(keepX[[q]])=max(ncomp)
        
        for(q in 1:length(X))
        {
            if(q<=length(keepX))
            {keepA[[q]]=keepX[[q]]
                if(length(keepA[[q]])<max(ncomp)) {keepA[[q]]=c(keepA[[q]],rep(ncol(X[[q]]),max(ncomp)-length(keepA[[q]])))} #complete the keepX already provided

            }else{
                keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)-length(keepX.constraint[[q]]))
            }
        }
        
    }
    
    keepA.constraint=keepX.constraint
    
}

    #check that keepX is ok
    for(q in 1:length(X))
    {
        if(is.list(keepA[[q]]))
        stop(paste0("keepX[[",q,"]]' must be a vector"))
        if (any(keepA[[q]] > ncol(X[[q]])))
        stop(paste0("each component of 'keepX[[",q,"]]' must be lower or equal to ", ncol(X[[q]]), "."))
        if(length(keepA[[q]])>max(ncomp))
        stop(paste0("length of 'keepX[[",q,"]]' must be lower or equal to ", ncomp[q], "."))
    }
    

    #check that length(keepA.constraint)+length(keepA)=ncomp for all blocks
    for(q in 1:length(X))
    {
        if (any(keepA[[q]] > ncol(X[[q]])))
        stop(paste0("each component of 'keepA[[",q,"]]' must be lower or equal to ", ncol(X[[q]]), "."))
        #check keepX and keepX.constraint
        if((length(keepA.constraint[[q]])+length(keepA[[q]]))!=max(ncomp)) stop(paste0("length (keepX.constraint[[",q,"]]) + length(keepX[[",q,"]]) should be max(ncomp)."))
    }
    
    #print(keepA)
    #print("constr")
    #print(keepA.constraint)
    
    
    return(list(keepA=keepA,keepA.constraint=keepA.constraint))
}





