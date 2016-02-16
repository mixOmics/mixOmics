# Copyright (C) 2013
# Kim-Anh Le Cao, University of Queensland, Brisbane, Australia
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

wrapper.sgccda = function(
X,
Y,
design = NULL,
ncomp = rep(1, length(X)),
keepA,
keepA.constraint,
scheme = "centroid",
mode="canonical",
scale = TRUE,
init = "svd",
bias = TRUE,
max.iter = 500,
tol = 1e-06,
verbose = FALSE
){
    
    ### Start: check X and Y matrix
    if (!(is.matrix(X) || is.list(X)))
    stop ("X must be either a matrix or a list")
    
    if (is.matrix(X))
    X = list(X = X)
    
    X = lapply(X, function(x) {as.matrix(x, rownames.force = ifelse(is.null(row.names(x)), FALSE, TRUE))})
    
    if (any(sapply(X, function(x){length(dim(x)) != 2 || !is.numeric(x)})))
    stop("each 'X' matrix must be a numeric matrix.")
    
    if (length(unique(sapply(X, nrow))) != 1)
    stop ("unequal number of rows in 'X'")
    
    if (is.null(dim(Y))) {
        Y = as.factor(Y)
    } else {
        stop("'Y' should be a factor or a class vector.")
        
    }
    
    
    if (unique(sapply(X, nrow)) != length(Y))
    stop("unequal number of rows in 'X' and 'Y'.")
    ### End: check X and Y matrix
    
    ### Start check design matrix
    if (is.null(design)) {
        design = 1 - diag(length(X) + 1)
    } else if (ncol(design) != nrow(design) || ncol(design) < length(X) || ncol(design) > (length(X) + 1) || any(!design %in% c(0,1))){
        stop('invalid design matrix.')
    } else if (ncol(design) == length(X)){
        warning('Design matrix changed')
        design = rbind(cbind(design, 1), 1)
        diag(design) = 0
    }
    ### End check design matrix
    
    ### Start Check keepX, ncomp and penalty parameters
    
    #check length(ncomp)=length(X)
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
    if (!is.null(ncomp)) {
        if (length(ncomp) != length(X) || !is.numeric(ncomp) || any(ncomp <= 0) || any(ncomp > sapply(X, ncol))){stop("invalid number of components, 'ncomp'.")}
    } else {
        ncomp = sapply(X, ncol)
    }
    
    ### End check ncomp parameter
    
    #add names to the blocks if no names or not unique name for each block
    if(length(unique(names(X)))!=length(X))
    names(X)=paste0("block",1:length(X))
    
    
    
    check=check.keepA.and.keepA.constraint(X=X,keepX=keepA,keepX.constraint=keepA.constraint,ncomp=ncomp)
    keepA=check$keepA
    keepA.constraint=check$keepA.constraint
    

    #A = X; A$Y = unmap(Y)
    #keepA$Y = rep(nlevels(Y), max(ncomp))
    
    ### End Check keepX and penalty parameter
    
    Y.input=Y
    Y=unmap(Y)
    
    
    result <- wrapper.sparse.meta.block(X=X,Y=Y,ncomp=ncomp,keepX.constraint=keepA.constraint,
    keepX=keepA,scheme=scheme,mode=mode,scale=scale,
    bias=bias,init=init,tol=tol,verbose=verbose,max.iter=max.iter,near.zero.var=FALSE)
    
    
    
    
    
    cl = match.call()
    cl[[1]] = as.name("sgccda")
    result$ind.mat=result$Y[[1]]
    result$Y=Y.input    
    result$call = cl
    result$ncomp = ncomp
    result$names$Y = attr(result$Y[[1]], "levels")
    row.names(result$variates$Y) = row.names(X); row.names(result$loadings$Y) = paste0("Y", c(1 : nlevels(Y.input)))
    names(result)[names(result) == "keepA"] = "keepX"; result$keepX = keepA
    class(result) = "sgccda"
    return(invisible(result))
    
    
}


