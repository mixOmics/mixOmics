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

wrapper.rgcca = function(
X,
design = 1 - diag(length(X)),
tau = rep(1, length(X)),
ncomp = rep(1, length(X)),
scheme = "centroid",
mode="canonical",
scale = TRUE,
init = "svd.single",
bias = FALSE,
tol = 1e-6,
verbose = FALSE,
max.iter=500,
near.zero.var = FALSE)
{
    
    # call function
    #rgcca <- function(A, C = 1-diag(length(A)), tau = rep(1, length(A)), ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE , init="svd", bias = TRUE, tol = .Machine$double.eps, verbose=TRUE)
    
    
    
    # check needed
    check=Check.entry.rgcca(X = X, design = design, tau = tau, ncomp = ncomp, scheme = scheme, scale = scale,
    init = init, bias = bias, tol = tol, verbose = verbose, max.iter=max.iter, mode=mode,near.zero.var=near.zero.var)
    X=check$A
    ncomp=check$ncomp
    design=check$design
    init=check$init
    scheme=check$scheme
    verbose=check$verbose
    bias=check$bias
    nzv.A=check$nzv.A

    
    result.rgcca = internal_mint.block(A = X, design = design, tau = tau,
    ncomp = ncomp,
    scheme = scheme, scale = scale,
    init = init, bias = bias, tol = tol, verbose = verbose,
    max.iter=max.iter,
    study=factor(rep(1,nrow(X[[1]]))),#mint.rgcca not coded yet
    mode="canonical"
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
    cl[[1]] = as.name('rgcca')
    
    output = list(
    class = cl,
    X = X,
    variates = result.rgcca$variates,
    loadings = result.rgcca$loadings,
    loadings.star = result.rgcca$loadings.star,
    design = design,
    tau = result.rgcca$tau,
    scheme = scheme,
    ncomp = ncomp,
    crit = result.rgcca$crit,
    AVE = list(AVE.X = result.rgcca$AVE$AVE_X, result.rgcca$AVE$AVE_outer, result.rgcca$AVE$AVE_inner), #rename?
    names = list(indiv = rownames(X[[1]]), var = sapply(X, colnames)),
    nzv=result.rgcca$nzv
    
    )
    
    class(output) = 'rgcca'
    return(invisible(output))
    
}
