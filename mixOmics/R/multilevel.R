# Copyright (C) 2015 
# Benoit Liquet, Universite de Bordeaux, France
# Kim-Anh Le Cao, University of Queensland, Brisbane, Australia
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
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
  function(X, 
           Y = NULL,
           design,  # design = NULL #*KA* added a stop missing below, because if missing, the whole multilevel module has no point
           ncomp = 2,
           keepX = NULL,
           keepY = NULL,  
           method = c("spls", "splsda"),  # *KA*
           mode = c("regression", "canonical"), # *KA*
           max.iter = 500, # *KA*
           tol = 1e-06, # *KA*
           ...)
{
    #-- checking general input arguments ---------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- check that the user did not enter extra arguments
    user.arg = names(match.call())[-1]
    function.arg = c("X", "Y", "design", "ncomp",
                     "keepX", "keepY", "method", "mode",
                     "max.iter", "tol", "near.zero.var",
                     "freqCut", "uniqueCut")
    
    in.arg = !(user.arg %in% function.arg)
    
    if (any(in.arg)) {
      unused.arg = user.arg[in.arg]
      unused.arg = paste(unused.arg, collapse = ", ")
      unused.arg = paste0("(", unused.arg, ").")
      stop("unused argument(s) ", unused.arg, call. = FALSE)
    }
    
    #-- method
    choices = c("spls", "splsda")
    method = choices[pmatch(method, choices)]
    
    if (is.na(method) || length(method) != 1) # *BG* See Ignacio email 18/02/2015, another way to deal with "method" pb 
      stop("'method' should be one of 'spls' or 'splsda'.", call. = FALSE)
    
    #-- mode
    if (method == "spls") {
      choices = c("regression", "canonical")
      mode = choices[pmatch(mode, choices)]
      
      if (is.na(mode) || length(mode) != 1) # *BG* See Ignacio email 18/02/2015, another way to deal with "mode" pb 
        stop("'mode' should be one of 'regression' or 'canonical'.", call. = FALSE)
    }
        
    #-- design
    # *KA* added this condition
    if(missing(design)) stop('You forgot to input the design matrix.', call. = FALSE) # *BG* add ", call. = FALSE"
    
    if(is.null(design)) 
      stop("the 'design' matrix is missing.", call. = FALSE)
    design = as.data.frame(design)
    
    #-- X matrix
    # *BG* Start: simplify check X matrix
#     if (is.data.frame(X)) 
#     	X = as.matrix(X)
    
#     if (!is.matrix(X)) {
#       if (!is.vector(X) || is.list(X))
#         stop("'X' must be a numeric matrix.", call. = FALSE)
#     }
    
#     X = as.matrix(X)
    
    if (is.data.frame(X) || !is.vector(X)) {
     	X = as.matrix(X)
    } else {
    	stop("'X' must be a numeric matrix.", call. = FALSE)
    }
     # *BG* End: simplify check X matrix
    
    if (is.character(X))
      stop("'X' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(X, 1, is.infinite))) 
      stop("infinite values in 'X'.", call. = FALSE)
    
    #-- equal number of rows in X and design
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
    # if keepX is missing, we assume it is a non sparse model
    if (is.null(keepX)) {
      keepX = rep(p, ncomp)
    } else {
      if (!is.numeric(keepX) || !is.vector(keepX))
        stop("'keepX' must be a numeric vector of length equal to ", ncomp, ".",
             call. = FALSE)
      
      if(any(keepX <= 0) || any(!is.finite(keepX)))
        stop("'keepX' must be positive integer vector.", call. = FALSE)
    }
    
    #-- Y matrix and keepY
    if (method == "spls") {
      
      # *BG* Start simplify check Y matrix
#      if (is.data.frame(Y)) 
#      	Y = as.matrix(Y)
      
#      if (!is.matrix(Y)) {
#        if (!is.vector(Y) || is.list(Y))
#          stop("'Y' must be a numeric matrix.", call. = FALSE)
#      }
      
#      Y = as.matrix(Y)
      
      if (is.data.frame(Y) || is.vector(Y)){ 
      	Y = as.matrix(Y)
      } else {
      	stop("'Y' must be a numeric matrix.", call. = FALSE)
      }
      # *BG* End simplify check Y matrix
      
      if (is.character(Y))
        stop("'Y' must be a numeric matrix.", call. = FALSE)
      
      if (any(apply(Y, 1, is.infinite))) 
        stop("infinite values in 'Y'.", call. = FALSE)
      
      q = ncol(Y)
      
      # if keepY is missing, we assume it is a non sparse model
      if (is.null(keepY)) {
        # keepX = rep(q, ncomp) # *BG* see Ignacio email 19/02/2015
        keepY = rep(q, ncomp) # *BG* see Ignacio email 19/02/2015
      } else {
        if (!is.numeric(keepY) || !is.vector(keepY))
          stop("'keepY' must be a numeric vector of length equal to ", ncomp, ".",
               call. = FALSE)
        
        if(any(keepY <= 0) || any(!is.finite(keepY)))
          stop("'keepY' must be positive integer vector.", call. = FALSE)
      }
    }
    
    #-- end checking --#
    #------------------#
    
    
    #-- multilevel approach ----------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    # design = apply(design, 2, as.factor) # *BG* break factor if present in design 
    
    # the first column of the design should indicate the sample ID
    # rep.measures = design[, 1] # *BG* new withinVariation function with only one argument called "design"
    # factors = as.data.frame(design[, -1]) # *BG* new withinVariation function with only one argument called "design"
    
    if (method == 'splsda') { 
      #-- apply multi-factor analysis with splsda --#
      #---------------------------------------------#
      
      # *BG* Start: reshape withinVariation using design argument
      #-- within-subject deviation matrix for 1 or 2 factors
      # if (ncol(factors) == 1){  
      #  Xw = suppressMessages(
      #    withinVariation(X, factors = NULL, rep.measures = rep.measures))
      #}
      #else{ 
      #  Xw = suppressMessages(
      #    withinVariation(X, factors = factors, rep.measures = rep.measures))
      #} # end 1 or 2 factor
      
      #-- apply sPLS-DA on the within deviation matrix --#
      # Y = apply(factors, 1, paste, collapse = "")
      
      #-- within-subject deviation matrix for 1 or 2 factors
      # remove suppressMessage to let the user know which design is analysed
      Xw = withinVariation(X, design  = design)
      
      #-- Need to set Y variable for 1 or 2 factors
      Y = design[, -1]
      if (!is.null(dim(Y))) {
        Y = apply(Y, 1, paste, collapse = ".")  #  paste is to combine in the case we have 2 levels
      }
      Y = as.factor(Y) 
      # *BG* End: reshape withinVariation using design argument      
      
      res = splsda(Xw, Y, ncomp = ncomp, keepX = keepX, max.iter = max.iter, tol = tol,...) # *KA* added max.iter and tol
      
      result = c(res, list(Xw = Xw, design = design))
      
      # *KA* different classes depending on the number of factors. Added to fit pheatmap (not the best patching)
      # if(ncol(factors) == 1){ # *BG* no feature factors anymore
      if(ncol(design) == 2){ # *BG* no feature factors anymore
        class(result) = c("mlsplsda", "splsda", "plsda", "splsda1fact") # *KA* added multilevel for pheatmap
      }else{  # if more than 2 columns (i.e. 2 factors)
        class(result) = c("mlsplsda", "splsda", "plsda", "splsda2fact")
      }
      ## class(result) = c("mlsplsda", "splsda", "plsda") 
    } else {
      #-- apply multilevel analysis with spls --#
      #-----------------------------------------#
      
      #-- within-subject deviation matrices --#
      # Y = as.matrix(Y) # *BG* already done before (check part)
      # Xw = suppressMessages(withinVariation(X, rep.measures = rep.measures)) # *KA* # *BG* use design argument with new withinVariation function
      # Yw = suppressMessages(withinVariation(Y, rep.measures = rep.measures)) # *KA* # *BG* use design argument with new withinVariation function
      Xw = withinVariation(X, design = design) # *BG* use design argument with new withinVariation function
      Yw = withinVariation(Y, design = design) # *BG* use design argument with new withinVariation function
      
      #-- apply sPLS on the within deviation matrices --#
      # if keepY is missing, we assume it is a non sparse model
      # line below is probably not needed since it is in the check function? # *KA*
      ##if (is.null(keepY)) keepY = rep(ncol(Y), ncomp)
      
      res = spls(Xw, Yw, ncomp = ncomp, mode = mode, 
                 keepY = keepY, keepX = keepX, max.iter = max.iter, tol = tol, ...) # *KA*
      
      result = c(res, list(Xw = Xw, Yw = Yw, design = design))
      class(result) = c("mlspls", "spls", "pls")
    } # end if spls / splsda
    
    return(invisible(result))
  }

