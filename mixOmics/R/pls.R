# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research, Toulouse France and 
# The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD 
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


pls <-
  function(X,
           Y,
           ncomp = 2,
           mode = "regression",
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
    
    #-- Y matrix
    if (is.data.frame(Y)) Y = as.matrix(Y)
    
    if (!is.matrix(Y)) {
      if (!is.vector(Y) || is.list(Y))
        stop("'Y' must be a numeric matrix.", call. = FALSE)
    }
    
    Y = as.matrix(Y)
    
    if (is.character(Y))
      stop("'Y' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(Y, 1, is.infinite))) 
      stop("infinite values in 'Y'.", call. = FALSE)
    
    #-- equal number of rows in X and Y
    if ((n = nrow(X)) != nrow(Y)) 
      stop("unequal number of rows in 'X' and 'Y'.", call. = FALSE)
    
    p = ncol(X)
    q = ncol(Y)
    
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
    
    #-- mode
    choices = c("regression", "canonical", "invariant", "classic")
    mode = choices[pmatch(mode, choices)]
   
    if (is.na(mode)) 
      stop("'mode' should be one of 'regression', 'canonical', 'invariant' or 'classic'.", 
           call. = FALSE)
    
    #-- max.iter
    if (is.null(max.iter) || max.iter < 1 || !is.finite(max.iter))
      stop("invalid value for 'max.iter'.", call. = FALSE)
    
    max.iter = round(max.iter)  
    
    #-- tol
    if (is.null(tol) || tol < 0 || !is.finite(tol))
      stop("invalid value for 'tol'.", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    
    #-- put names to variables and samples --#
    X.names = colnames(X)
    if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")    
    
    if (dim(Y)[2] == 1) Y.names = data.names[2]
    else {
      Y.names = colnames(Y)
      if (is.null(Y.names)) Y.names = paste("Y", 1:q, sep = "")
    }
    
    ind.names = rownames(X)
    if (is.null(ind.names)) {
      ind.names = rownames(Y)
      rownames(X) = ind.names
    }
    
    if (is.null(ind.names)) {
      ind.names = 1:n
      rownames(X) = rownames(Y) = ind.names
    }

    
    #-- pls approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- scale the matrices --#
    X = scale(X, center = TRUE, scale = TRUE)
    Y = scale(Y, center = TRUE, scale = TRUE) 
    
    #-- initialization of matrices --#
    X.temp = X
    Y.temp = Y
    mat.t = mat.u = matrix(nrow = n, ncol = ncomp)
    mat.a = mat.c = matrix(nrow = p, ncol = ncomp)
    mat.b = mat.d = mat.e = matrix(nrow = q, ncol = ncomp)
    n.ones = rep(1, n)
    p.ones = rep(1, p)
    q.ones = rep(1, q)
    na.X = FALSE
    na.Y = FALSE
    is.na.X = is.na(X)
    is.na.Y = is.na(Y)
    if (any(is.na.X)) na.X = TRUE
    if (any(is.na.Y)) na.Y = TRUE
    
    #-- loop on h --#
    for (h in 1:ncomp) {
      
      #-- initialization --#
      u = Y.temp[, 1] 
      if (any(is.na(u))) u[is.na(u)] = 0
      a.old = 0
      iter = 1
      
      if (na.X) {
        X.aux = X.temp
        X.aux[is.na.X] = 0
      }
      
      if (na.Y) {
        Y.aux = Y.temp
        Y.aux[is.na.Y] = 0
      }
      
      repeat {
        #-- compute loading vectors and variates associated to X 
        if (na.X) {
          a = crossprod(X.aux, u)				
          U = drop(u) %o% p.ones
          U[is.na.X] = 0
          u.norm = crossprod(U)				
          a = a / diag(u.norm)
          a = a / drop(sqrt(crossprod(a)))
          tt = X.aux %*% a
          A = drop(a) %o% n.ones
          A[t(is.na.X)] = 0
          a.norm = crossprod(A)
          tt = tt / diag(a.norm)				
        }
        else {			
          a = crossprod(X.temp, u)
          a = a / drop(sqrt(crossprod(a)))
          tt = X.temp %*% a / drop(crossprod(a))
        }
        
        #-- compute loading vectors and variates associated to Y 
        if (na.Y) {
          b = crossprod(Y.aux, tt)
          TT = drop(tt) %o% q.ones
          TT[is.na.Y] = 0
          t.norm = crossprod(TT)				
          b = b / diag(t.norm)
          u = Y.aux %*% b
          B = drop(b) %o% n.ones
          B[t(is.na.Y)] = 0
          b.norm = crossprod(B)
          u = u / diag(b.norm)					
        }
        else {			
          b = crossprod(Y.temp, tt) / drop(crossprod(tt)) 
          u = Y.temp %*% b / drop(crossprod(b))
        }
        
        if (crossprod(a - a.old) < tol) break
        
        if (iter == max.iter) {
          warning(paste("maximum number of iterations reached for dimension", h),
                  call. = FALSE)
          break
        }
        
        a.old = a
        iter = iter + 1
      }
      
      #-- deflation of matrices --#
      if (na.X) {
        X.aux = X.temp
        X.aux[is.na.X] = 0
        c = crossprod(X.aux, tt)				
        TT = drop(tt) %o% p.ones
        TT[is.na.X] = 0
        t.norm = crossprod(TT)				
        c = c / diag(t.norm)
      }
      else {
        c = crossprod(X.temp, tt) / drop(crossprod(tt))
      }	
      
      X.temp = X.temp - tt %*% t(c)   
      
      #-- mode canonique --#
      if (mode == "canonical") {
        if (na.Y) {
          Y.aux = Y.temp
          Y.aux[is.na.Y] = 0
          e = crossprod(Y.aux, u)
          U = drop(u) %o% q.ones
          U[is.na.Y] = 0
          u.norm = crossprod(U)				
          e = e / diag(u.norm)					
        }
        else {
          e = crossprod(Y.temp, u) / drop(crossprod(u))
        }
        
        Y.temp = Y.temp - u %*% t(e)
      }
      
      #-- mode classic --#
      if(mode == "classic") Y.temp = Y.temp - tt %*% t(b)                 
      
      #-- mode regression --#
      if(mode == "regression") {
        if (na.Y) {
          Y.aux = Y.temp
          Y.aux[is.na.Y] = 0
          d = crossprod(Y.aux, tt)
          TT = drop(tt) %o% q.ones
          TT[is.na.Y] = 0
          t.norm = crossprod(TT)				
          d = d / diag(t.norm)
        }
        else {				
          d = crossprod(Y.temp, tt) / drop(crossprod(tt))
        }
        
        Y.temp = Y.temp - tt %*% t(d)
      }
      
      #-- mode invariant --#
      if (mode == "invariant") Y.temp = Y 
      
      mat.t[, h] = tt
      mat.u[, h] = u
      mat.a[, h] = a
      mat.b[, h] = b
      mat.c[, h] = c
      if (mode == "regression") mat.d[, h] = d
      if (mode == "canonical") mat.e[, h] = e
      
    } #-- end loop on h --#
    
    
    #-- output -----------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    rownames(mat.a) = rownames(mat.c) = X.names
    rownames(mat.b) = Y.names
    rownames(mat.t) = rownames(mat.u) = ind.names
    
    comp = paste("comp", 1:ncomp)
    colnames(mat.t) = colnames(mat.u) = comp
    colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = comp 
    
    cl = match.call()
    cl[[1]] = as.name('pls')
    
    result = list(call = cl,
                  X = X, 
                  Y = Y, 
                  ncomp = ncomp, 
                  mode = mode, 
                  mat.c = mat.c,
                  mat.d = mat.d,
                  mat.e = mat.e, 
                  variates = list(X = mat.t, Y = mat.u),
                  loadings = list(X = mat.a, Y = mat.b), 
                  names = list(X = X.names, Y = Y.names, sample = ind.names,
                               data = data.names),
                  tol = tol,
                  max.iter = max.iter)
    
    if (near.zero.var == TRUE && length(nzv$Position > 0)) result$nzv = nzv  
    
    class(result) = "pls"
    return(invisible(result))
}

