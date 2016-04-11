#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2015
# last modified: 01-03-2016
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


# ========================================================================================================
# pls: perform a PLS
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

# X: numeric matrix of predictors
# Y: numeric vector or matrix of responses
# ncomp: the number of components to include in the model. Default to 2.
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations


pls <- function(X,
Y,
ncomp = 2,
scale = TRUE,
mode = c("regression", "canonical", "invariant", "classic"),
tol = 1e-06,
max.iter = 500,
near.zero.var = FALSE,
logratio="none",   # one of "none", "CLR"
multilevel=NULL)    # multilevel is passed to multilevel(design=) in withinVariation. Y is ommited and shouldbe included in multilevel design
{
    
    
    #-- validation des arguments --#
    # most of the checks are done in the wrapper.mint.spls.hybrid function
    if(missing(ncomp)) ncomp=2

    result <- internal_wrapper.mint(X=X,Y=Y,ncomp=ncomp,scale=scale,near.zero.var=near.zero.var,mode=mode,
    max.iter=max.iter,tol=tol,logratio=logratio,multilevel=multilevel,DA=FALSE)

    

    cl = match.call()
    #cl[[1]] = as.name("pls")
    
    out=list(call=cl,X=result$X[-result$indY][[1]],Y=result$X[result$indY][[1]],ncomp=result$ncomp,mode=result$mode,variates=result$variates,loadings=result$loadings,
        names=result$names,tol=result$tol,iter=result$iter,max.iter=result$max.iter,nzv=result$nzv,scale=scale,explained_variance=result$explained_variance)
     

    class(out) = c("pls")
    
    if(!is.null(multilevel))
    {
        out$Xw=result$Xw
        out$multilevel=multilevel
        class(out)=c("mlpls",class(out))
    }
    
    #calcul explained variance
    #explX=explained_variance(out$X,out$variates$X,ncomp)
    #explY=explained_variance(out$Y,out$variates$Y,ncomp)
    #out$explained_variance=list(X=explX,Y=explY)
    
    return(invisible(out))


}
