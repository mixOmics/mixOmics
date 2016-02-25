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



wrapper.sgcca = function(
X,
design = 1 - diag(length(X)),
penalty = NULL,
ncomp = rep(1, length(X)),
keepX,
keepX.constraint,
scheme = "centroid",
mode="canonical",
scale = TRUE,
init = "svd",
bias = TRUE,
tol = .Machine$double.eps,
verbose = FALSE,
max.iter=500,
near.zero.var = FALSE
){
  
  # call function
  #rgcca <- function(A, C = 1-diag(length(A)), tau = rep(1, length(A)), ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE , init="svd", bias = TRUE, tol = .Machine$double.eps, verbose=TRUE)
  
  
  check=Check.entry.sgcca(X=X, design=design ,ncomp=ncomp , scheme=scheme , scale=scale ,  bias=bias,
  init=init , tol=tol , verbose=verbose,mode=mode, max.iter=max.iter,near.zero.var=near.zero.var,keepX=keepX,keepX.constraint=keepX.constraint)
  

  A=check$A
  design=check$design
  ncomp=check$ncomp
  init=check$init
  scheme=check$scheme
  verbose=check$verbose
  bias=check$bias
  near.zero.var=check$near.zero.var
  keepA.constraint=check$keepA.constraint
  keepA=check$keepA
  nzv.A=check$nzv.A
  
  
  result.sgcca = internal_mint.block(A = A, design = design, tau = NULL,
                       ncomp = ncomp,
                       scheme = scheme, scale = scale,
                       init = init, bias = bias, tol = tol, verbose = verbose,
                       keepA.constraint=NULL,
                       keepA=keepA,
                       max.iter=max.iter,
                       study=factor(rep(1,nrow(A[[1]]))),#mint.rgcca not coded yet
                       mode=mode,penalty=penalty
                       )

  # outputs
#   out <- list(Y = shave.matlist(Y, ncomp),
#               a = shave.matlist(a, ncomp), 
#               astar = shave.matlist(astar, ncomp),
#               C = C, tau = tau_mat, scheme = scheme,
#               ncomp=ncomp, crit = crit,
#               mode = mode,
#               AVE=list(AVE_X=AVE_X,
#                        AVE_outer=AVE_outer,
#                        AVE_inner=AVE_inner),
#               #KA added names of rows and cols for plotIndiv and plotVar
#               names = list(indiv = rownames(A[[1]]))
#   )
#   class(out) <- "rgcca"
#   return(out)

  cl = match.call()
  cl[[1]] = as.name('sgcca')
  
  output = list(
    class = cl,
    X = X,
    variates = result.sgcca$variates,
    loadings = result.sgcca$loadings,
    loadings.star = result.sgcca$loadings.star,
    design = design,
    penalty=penalty,
    scheme = scheme,
    ncomp = ncomp, 
    crit = result.sgcca$crit,
    AVE = list(AVE.X = result.sgcca$AVE$AVE_X, result.sgcca$AVE$AVE_outer, result.sgcca$AVE$AVE_inner), #rename?
    names = list(indiv = rownames(X[[1]]), var = sapply(X, colnames)),
    nzv=result.sgcca$nzv
    
  )

  class(output) = 'sgcca'
  return(invisible(output))
  
}
