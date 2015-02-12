# Copyright (C) 2009 
# S?bastien D?jean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonz?lez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh L? Cao, French National Institute for Agricultural Research and 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
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


plsda <-
  function(X,
           Y, 
           ncomp = 2,
           max.iter = 500,  	 
           tol = 1e-06,
           near.zero.var = TRUE,
           ...)
{
    
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- check that the user did not enter extra arguments
    arg.call = match.call()
    user.arg = names(arg.call)[-1]
    
    err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())), 
                   error = function(e) e)
    
    if ("simpleError" %in% class(err))
      stop(err[[1]], ".", call. = FALSE)
    
    default.arg = c("freqCut", "uniqueCut")
    function.arg = c(names(mget(names(formals()), sys.frame(sys.nframe()))),
                     default.arg)
    not.arg = !(user.arg %in% function.arg)
    
    if (any(not.arg)) {
      unused.arg = user.arg[not.arg]
      not.arg = which(not.arg) + 1
      output = rep("", length(not.arg))
      
      for (i in 1:length(not.arg)) {
        output[i] = paste0(unused.arg[i], " = ", arg.call[[not.arg[i]]])
      }
      
      output = paste0("(", paste(output, collapse = ", "), ").")
      msg = "unused argument "
      if (length(not.arg) > 1) msg = "unused arguments "  
      stop(msg, output, call. = FALSE)
    }
    
    #-- data set names --#
    data.names = c(deparse(substitute(X)), deparse(substitute(Y)))
    
    #-- X matrix
    if (is.data.frame(X)) X = as.matrix(X)
    
    if (!is.matrix(X)) {
      if (!is.vector(X) || is.list(X))
        stop("'X' must be a numeric matrix.", call. = FALSE)
    }
    
    X = as.matrix(X)
    
    if (is.character(X))
      stop("'X' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(X, 1, is.infinite))) 
      stop("infinite values in 'X'.", call. = FALSE)
    
    #-- Y factor
    if (is.null(dim(Y))) {
      Y = as.factor(Y)  
      ind.mat = unmap(as.numeric(Y))					
    }
    else {
      stop("'Y' should be a factor or a class vector.", 
           call. = FALSE)						
    }		
    
    #-- equal number of rows in X and Y
    if ((n = nrow(X)) != nrow(ind.mat)) 
      stop("unequal number of samples in 'X' and 'Y'.", call. = FALSE)
    
    #-- near.zero.var
    if (!is.logical(near.zero.var))
      stop("'near.zero.var' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)
    
    #-- ncomp
    if (is.null(ncomp) || !is.finite(ncomp) || ncomp <= 0)
      stop("invalid number of components, 'ncomp'.", call. = FALSE)
    
    if (near.zero.var == TRUE) { 
      nzv = nearZeroVar(X, ...)
      
      if (length(nzv$Position > 0)) {
        warning("zero- or near-zero variance predictors. \nReset predictors matrix to not near-zero variance predictors. \nSee $nzv for problematic predictors.", 
                call. = FALSE)
        
        X = X[, -nzv$Position]
      }
    }
    
    p = ncol(X)
    ncomp = round(ncomp)
    
    if (ncomp > p) {
      warning("reset maximum number of variates 'ncomp' to ncol(X) = ", p, ".", 
              call. = FALSE)
      
      ncomp = p
    }
    
    #-- max.iter
    if (is.null(max.iter) || max.iter < 1 || !is.finite(max.iter))
      stop("invalid value for 'max.iter'.", call. = FALSE)
    
    max.iter = round(max.iter)  
    
    #-- tol
    if (is.null(tol) || tol < 0 || !is.finite(tol))
      stop("invalid value for 'tol'.", call. = FALSE)
    
    #-- end checking --#
    #------------------#

    
    #-- call pls approach ------------------------------------------------------#
    #---------------------------------------------------------------------------#

    result = pls(X, ind.mat, ncomp = ncomp, mode = "regression", 
                 max.iter = max.iter, tol = tol, near.zero.var = FALSE, 
                 ...)
    
    
    #-- output -----------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    cl = match.call()
    cl[[1]] = as.name('plsda')
    result$call = cl
    
    result$ind.mat = ind.mat
    result$names$Y = levels(Y)
    result$tol = tol
    result$max.iter = max.iter
    
    if (near.zero.var == TRUE && length(nzv$Position > 0)) result$nzv = nzv
    
    class(result) = "plsda"
    return(invisible(result))
}
