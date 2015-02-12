# Copyright (C) 2014 
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Amrit Singh, University of British Columbia, Vancouver.
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
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


# ---------------------------------------------------
# perf for pls object
# ---------------------------------------------------
perf.pls <-
  function(object,
           criterion = c("all", "MSEP", "R2", "Q2"), 
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE,
           ...)
  {
    
    #-- validation des arguments --#
    # these are the centered and scaled matrices output from pls 
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    
    if (length(dim(X)) != 2) 
      stop("'X' must be a numeric matrix for validation.")
    
    if(object$mode == 'canonical') stop('PLS mode should be set to regression, invariant or classic')
    
    validation = match.arg(validation)
    
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    if (any(criterion == "all" | criterion == "Q2") & ncomp == 1)
      stop("'ncomp' must be > 1 for Q2 criterion.")
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    
    #-- define  M fold or loo cross validation --------------------#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }
      else {
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    #-- compute MSEP and/or R2 --#
    
    ## add Q2
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = array(0, c(n, q, ncomp+1))
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    
    #KA: all criteria included in the computation
    if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))){
      press.mat = Ypred = array(0, c(n, q, ncomp))
      MSEP = R2 = matrix(0, nrow = q, ncol = ncomp)
      
      # set up dimnames
      rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
      dimnames(press.mat)[[2]] = colnames(Y)
      
      # in case the test set only includes one sample, it is better to advise the user to perform loocv
      stop.user = FALSE
      if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
      
      for (i in 1:M) {
        if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
        
        omit = folds[[i]]
        # see below, we stop the user if there is only one sample drawn on the test set using MFold
        if(length(omit) == 1) stop.user = TRUE
        
        # the training set is NOT scaled
        X.train = X[-omit, ]
        Y.train = Y[-omit, ]
        # the test set is scaled either in the predict function directly (for X.test)
        # or below for Y.test
        X.test = matrix(X[omit, ], nrow = length(omit))
        Y.test = matrix(Y[omit, ], nrow = length(omit))
        
        
        #-- pls --#
        pls.res = pls(X = X.train, Y = Y.train, ncomp = ncomp, 
                      mode = mode, max.iter = max.iter, tol = tol)
        
        if (!is.null( pls.res$nzv$Position)) X.test = X.test[, - pls.res$nzv$Position]
        # in the predict function, X.test is already normalised w.r.t to training set X.train, so no need to do it here
        Y.hat = predict( pls.res, X.test)$predict
        
        for (h in 1:ncomp) {
          Ypred[omit, , h] = Y.hat[, , h]
          
          # compute the press and the RSS
          # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
          press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
          RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
          
        } # end h
      } #end i (cross validation)
      
      # warn the user that at least test set had a length of 1
      if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
      
      # these criteria are computed across all folds.
      for (h in 1:ncomp) { 
        MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[, , h], use = "pairwise")))^2
        
        # on en profite pour calculer le PRESS par composante et le Q2
        if(q>1){
          RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
          PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
          RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
          PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
      }
      
      colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
      
      if (q == 1){
        rownames(MSEP) = rownames(R2) = ""
      }
      
      if(ncomp>1){
        # compute Q2 total
        if(q>1){
          Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
        }else{ # for q == 1
          Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
        }
      }else{
        Q2.total = NA
      }
      
      names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
      
    } # end all, MSEP, R2, Q2
    
    
    if (progressBar == TRUE) cat('\n')
    
    # ----------  final outputs:
    if (any(criterion %in% c("all", "MSEP"))) res$MSEP = MSEP
    if (any(criterion %in% c("all", "R2"))) res$R2 = R2
    if (any(criterion %in% c("all", "Q2"))) res$Q2 = t(Q2.inside)
    if (any(criterion %in% c("all", "Q2"))) res$Q2.total = Q2.total
    
    res$press.mat=press.mat
    res$RSS.indiv=RSS.indiv
    res$PRESS.inside=PRESS.inside
    res$RSS=RSS
    
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
  }

