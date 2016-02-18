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
# meta.spls: perform a vertical PLS on a combination of experiments, input as a matrix in X
# this function is a particular setting of meta.spls.hybrid, the formatting of the input is checked in wrapper.meta.spls.hybrid
# ========================================================================================================

# X: numeric matrix of predictors
# Y: numeric vector or matrix of responses
# ncomp: the number of components to include in the model. Default to 2.
# study: grouping factor indicating which samples are from the same study
# keepX.constraint: A list containing which variables of X are to be kept on each of the first PLS-components.
# keepY.constraint: A list containing which variables of Y are to be kept on each of the first PLS-components
# keepX: number of \eqn{X} variables kept in the model on the last components (once all keepX.constraint[[i]] are used).
# keepY: number of \eqn{Y} variables kept in the model on the last components.
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations


meta.spls <- function(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
keepX.constraint=NULL,
keepY.constraint=NULL,
keepX=rep(ncol(X), ncomp),
keepY=rep(ncol(Y), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 500,
near.zero.var = FALSE)
{

    result <- wrapper.meta.spls.hybrid(X=X,Y=Y,ncomp=ncomp,scale=scale,near.zero.var=near.zero.var,study=study,mode=mode,
    keepX=keepX,keepY=keepY,keepX.constraint=keepX.constraint,keepY.constraint=keepY.constraint,max.iter=max.iter,tol=tol)
    
    
    cl = match.call()
    cl[[1]] = as.name("meta.spls")
    
    out=list(call=cl,X=result$X[[1]],Y=result$Y[[1]],ncomp=result$ncomp,study=result$study,mode=result$mode,keepX=result$keepA[[1]],keepY=result$keepA[[2]],
    keepX.constraint=result$keepA.constraint[[1]],keepY.constraint=result$keepA.constraint[[2]],
    variates=result$variates,loadings=result$loadings,variates.partial=result$variates.partial,loadings.partial=result$loadings.partial,
    names=result$names,tol=result$tol,iter=result$iter,nzv=result$nzv,scale=scale)
    
    class(out) = "meta.spls"
    return(invisible(out))
 
    
    
    
    
}
