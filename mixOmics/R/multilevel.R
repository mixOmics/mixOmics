# Copyright (C) 2012 
# Benoit Liquet, Université de Bordeaux, France
# Kim-Anh Lê Cao, Queensland Facility for Advanced Bioinformatics, University of Queensland, Brisbane, Australia
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
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


# ---------------------------------------------
# Generic function for multilevel analysis 
# ---------------------------------------------

multilevel <- 
  function(X, Y = NULL,
           design = NULL,
           ncomp = 2,
           keepX = NULL,
           keepY = NULL,  
           method = NULL,
           mode = "regression",
           ...)
  {
    #-- checking general input arguments ---------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- check that the user did not enter extra arguments
    arg.call = match.call()
    user.arg = names(arg.call)[-1]
    
    err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())), 
                   error = function(e) e)
    
    if ("simpleError" %in% class(err))
      stop(err[[1]], ".", call. = FALSE)
    
    default.arg = c("max.iter", "tol", "near.zero.var",
                    "freqCut", "uniqueCut")
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
    
    #-- method
    choices = c("spls", "splsda")
    method = choices[pmatch(method, choices)]
    
    if (is.na(method)) 
      stop("'method' should be one of 'spls' or 'splsda'.", call. = FALSE)
    
    #-- mode
    if (method == "spls") {
      choices = c("regression", "canonic")
      mode = choices[pmatch(mode, choices)]
      
      if (is.na(mode)) 
        stop("'mode' should be one of 'regression' or 'canonic'.", call. = FALSE)
    }
        
    #-- design
    if(is.null(design)) 
      stop("the 'design' matrix is missing.", call. = FALSE)
    
    design = as.data.frame(design)
    
    if ((nrow(design) != nrow(X))) 
      stop("unequal number of rows in 'X' and 'design'.", call. = FALSE)
    
    if (ncol(design) < 2) 
      stop("'design' must be a matrix or data frame with at least 2 columns.",
           call. = FALSE)
    
    #-- ncomp
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
      stop("invalid number of variates, 'ncomp'.", call. = FALSE)
    
    p = ncol(X)
    ncomp = round(ncomp)
    
    if(ncomp > p) {
      warning("reset maximum number of variates 'ncomp' to ncol(X) = ", p, ".", 
              call. = FALSE)
      ncomp = p
    }
    
    #-- keepX
    if (is.null(keepX)) {
      keepX = rep(p, ncomp)
    }
    else {
      if (!is.numeric(keepX) || !is.vector(keepX))
        stop("'keepX' must be a numeric vector of length equal to ", ncomp, ".",
             call. = FALSE)
      
      if (length(keepX) != ncomp) 
        stop("length of 'keepX' must be equal to ", ncomp, ".", call. = FALSE)
      
      if (any(keepX > p)) 
        stop("each component of 'keepX' must be lower or equal than ", p, ".",
             call. = FALSE)
    }
    
    #-- Y matrix and keepY if method = 'spls'
    if (method == "spls") {
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
      
      q = ncol(Y)
      
      #-- keepY
      if (is.null(keepY)) {
        keepX = rep(q, ncomp)
      }
      else {
        if (!is.numeric(keepY) || !is.vector(keepY))
          stop("'keepY' must be a numeric vector of length equal to ", ncomp, ".",
               call. = FALSE)
        
        if (length(keepY) != ncomp) 
          stop("length of 'keepY' must be equal to ", ncomp, ".", call. = FALSE)
        
        if (any(keepY > q)) 
          stop("each component of 'keepY' must be lower or equal than ", q, ".",
               call. = FALSE)
      }
    }
    
    #-- end checking --#
    #------------------#
    
    
    #-- multilevel approach ----------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    desing = apply(design, 2, as.factor)
    rep.measures = design[, 1]
    factors = as.data.frame(design[, -1])
    
    #-- within-subject deviation matrix function for 1 factor
    within.dev <- function(X, rep.measures) {
      Xb = apply(X, 2, tapply, rep.measures, mean)
      Xb = Xb[rep.measures, ]
      Xw = X - Xb
    }
    
    if(method == 'splsda') { 
      #-- apply multi-factor analysis with splsda --#
      #---------------------------------------------#
      
      #-- within-subject deviation matrix for 1 or 2 factors
      if (ncol(factors) == 1) {
        Xw = within.dev(X, rep.measures)
      }
      else { #-- if 2 factors 
        
        ###### off set term
        Xm = colMeans(X)
        
        ###### compute the between subject variation 
        Xs = apply(X, 2, tapply, rep.measures, mean)
        Xs = Xs[rep.measures, ]
        
        xbfact1 = apply(X, 2, tapply, paste0(rep.measures, factors[, 1]), mean)
        xbfact1 = xbfact1[paste0(rep.measures, factors[, 1]), ]
        
        xbfact2 = apply(X, 2, tapply, paste0(rep.measures, factors[, 2]), mean)
        xbfact2 = xbfact2[paste0(rep.measures, factors[, 2]), ]
        
        ##### fixed effect
        Xmfact1 = apply(X, 2, tapply, factors[, 1], mean)
        Xmfact1 = Xmfact1[factors[, 1], ]
        
        Xmfact2 = apply(X, 2, tapply, factors[, 2], mean)
        Xmfact2 = Xmfact2[factors[, 2], ]
        
        Xw = X + Xs - xbfact1 - xbfact2 + Xmfact1 + Xmfact2
        Xw = sweep(Xw, 2, 2*Xm)
      }
      
      #-- apply sPLS-DA on the within deviation matrix --#
      Y = apply(factors, 1, paste, collapse = "")
      res = splsda(Xw, Y, ncomp = ncomp, keepX = keepX, ...)
      
      result = c(res, list(Xw = Xw, design = design))
      class(result) = c("mlsplsda", "plsda")
    }
    else {
      #-- apply multilevel analysis with spls --#
      #-----------------------------------------#
      
      #-- within-subject deviation matrices --#
      Y = as.matrix(Y)
      Xw = within.dev(X, rep.measures)
      Yw = within.dev(Y, rep.measures)
      
      #-- apply sPLS on the within deviation matrices --#
      if (is.null(keepY)) keepY = rep(ncol(Y), ncomp)
      
      res = spls(Xw, Yw, ncomp = ncomp, mode = mode, keepY = keepY, keepX = keepX, ...) 
      
      result = c(res, list(Xw = Xw, Yw = Yw, design = design))
      class(result) = c("mlspls", "pls")
    }
    
    result$names$data = data.names
    return(invisible(result))
  }


#-------------------------------------------------------
# within-subject deviation matrix function for 1 factor
#-------------------------------------------------------
Split.variation.one.level <- 
  function(X, rep.measures) {
    rownames(X) = as.character(rep.measures)
    Xb = apply(X, 2, tapply, rep.measures, mean)
    Xb = Xb[as.character(rep.measures), ]
    Xw = X - Xb
    return(invisible(list(Xw = Xw)))
}


#-----------------------------------------------
# within-subject deviation matrix for 2 factors
#-----------------------------------------------
Split.variation.two.level <- 
  function(X, factor1, factor2, rep.measures) {
    
    factors = data.frame(factor1, factor2)
    
    ###### off set term
    Xm = colMeans(X)
    
    ###### compute the between subject variation 
    Xs = apply(X, 2, tapply, rep.measures, mean)
    Xs = Xs[rep.measures, ]
    
    xbfact1 = apply(X, 2, tapply, paste0(rep.measures, factors[, 1]), mean)
    xbfact1 = xbfact1[paste0(rep.measures, factors[, 1]), ]
    
    xbfact2 = apply(X, 2, tapply, paste0(rep.measures, factors[, 2]), mean)
    xbfact2 = xbfact2[paste0(rep.measures, factors[, 2]), ]
    
    #### fixed effect
    Xmfact1 = apply(X, 2, tapply, factors[, 1], mean)
    Xmfact1 = Xmfact1[factors[, 1], ]
    
    Xmfact2 = apply(X, 2, tapply, factors[, 2], mean)
    Xmfact2 = Xmfact2[factors[, 2], ]
    
    Xw = X + Xs - xbfact1 - xbfact2 + Xmfact1 + Xmfact2
    Xw = sweep(Xw, 2, 2*Xm)
    return(invisible(list(Xw = Xw)))
}
