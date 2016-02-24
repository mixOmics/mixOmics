# Copyright (C) 2015
# Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# created: 22-04-2015
# last modified: 18-02-2016

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


# ========================================================================================================
# block.spls: perform a horizontal PLS-DA on a combination of datasets, input as a list in X
# this function is a particular setting of sparse.mint.block, the formatting of the input is checked in wrapper.sparse.mint.block
# ========================================================================================================

# X: a list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets. \code{NA}s are not allowed.
# Y: outcome.
# indY: to supply if Y is missing, indicate the position of the outcome in the list X.
# ncomp: numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.
# keepX.constraint: A list of same length as X. Each entry keepX.constraint[[i]] is a list containing which variables of X[[i]] are to be kept on each of the first PLS-components
# keepY.constraint: Only used if Y is provided, otherwise extracted from keepX.constraint. A list containing which variables of Y are to be kept on each of the first PLS-components
# keepX: A vector of same length as X.  Each entry keepX[i] is the number of X[[i]]-variables kept in the model on the last components (once all keepX.constraint[[i]] are used).
# keepY: Only used if Y is provided. Each entry keepY[i] is the number of Y-variables kept in the model on the last components.
# design: the input design.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid". Default to "centroid"
# mode: input mode, one of "canonical", "classic", "invariant" or "regression". Default to "regression"
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# bias: boleean. A logical value for biaised or unbiaised estimator of the var/cov (defaults to FALSE).
# init: intialisation of the algorithm, one of "svd" or "svd.single". Default to "svd"
# tol: Convergence stopping value.
# verbose: if set to \code{TRUE}, reports progress on computing.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations


block.spls <- function(X,
Y,
indY,
ncomp=rep(2,length(X)),
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
    
    
    result <- wrapper.sparse.mint.block(X=X,Y=Y,indY=indY,ncomp=ncomp,keepX.constraint=keepX.constraint,
    keepY.constraint=keepY.constraint,keepX=keepX,keepY=keepY,design=design,scheme=scheme,mode=mode,scale=scale,
    bias=bias,init=init,tol=tol,verbose=verbose,max.iter=max.iter,near.zero.var=near.zero.var)
    
    
    cl = match.call()
    cl[[1]] = as.name("block.spls")
    
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
    
    
    out=list(call=cl,X=result$X,Y=result$Y[[1]],ncomp=result$ncomp,mode=result$mode,keepX=keepX,keepY=keepY,
    keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,
    variates=result$variates,loadings=result$loadings,names=result$names,
    tol=result$tol,iter=result$iter,nzv=result$nzv,scale=scale)

    if(!missing(ncomp))   out$ncomp=ncomp


    class(out) = "block.spls"
    return(invisible(out))
    
}



