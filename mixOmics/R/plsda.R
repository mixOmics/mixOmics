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


# ========================================================================================================
# plsda: perform a PLS-DA
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 2.
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations


plsda <- function(X,
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
    if(is.null(multilevel))
    {
        if (is.null(Y))
        stop("'Y' has to be something else than NULL.")
        
        if (is.null(dim(Y)))
        {
            Y = as.factor(Y)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        Y.mat=unmap(Y)
        colnames(Y.mat) = levels(Y)#paste0("Y", 1:ncol(Y.mat))
    }else{Y.mat=NULL}

    
    result <- internal_wrapper.mint(X=X,Y=Y.mat,ncomp=ncomp,scale=scale,near.zero.var=near.zero.var,mode=mode,
    max.iter=max.iter,tol=tol,logratio=logratio,multilevel=multilevel,DA=TRUE)

    
    cl = match.call()
    #cl[[1]] = as.name("plsda")
    
    out=list(call=cl,X=result$X[[1]],Y=if(is.null(multilevel)){Y}else{result$Y.factor},ind.mat=result$Y[[1]],ncomp=result$ncomp,
    mode=result$mode,variates=result$variates,loadings=result$loadings,
    names=result$names,tol=result$tol,iter=result$iter,nzv=result$nzv,scale=scale)
    out$names$Y = levels(out$Y)
    #row.names(out$variates$Y) = row.names(out$variates$X)
    #row.names(out$loadings$Y) = paste0("Y", c(1 : nlevels(out$Y)))
    
   
    class(out) = c("plsda","pls")

    if(!is.null(multilevel))
    {
        out$Xw=result$Xw
        out$multilevel=multilevel
        class(out)=c("mlplsda",class(out))
    }

    return(invisible(out))
    
    
    
    
}
