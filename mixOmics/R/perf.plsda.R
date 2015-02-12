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
# perf for plsda object
# ---------------------------------------------------
perf.plsda <-
  function(object,
           method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE, ...) 
  {
    
    #-- validation des arguments --#
    # these data are the centered and scaled X output or the unmapped(Y) scaled and centered
    X = object$X
    level.Y = object$names$Y  # to make sure the levels are ordered
    Y = object$ind.mat
    Y = map(Y)
    Y = factor(Y,labels = level.Y)
    ncomp = object$ncomp
    n = nrow(X)
    
    tol = object$tol
    max.iter = object$max.iter
    
    method.predict = match.arg(method.predict, several.ok = TRUE)
    if (any(method.predict == "all")) {
      nmthdd = 3
    }else{
      nmthdd = length(method.predict)
    }
    
    
    error.fun = function(x, y) {
      error.vec = sweep(x, 1, y, FUN = "-")
      error.vec = (error.vec != 0)
      error.vec = apply(error.vec, 2, sum) / length(y)
      return(error.vec)
    }
    
    #-- define the folds --#
    if (validation == "Mfold") {
      if (is.list(folds)) {
        if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds.")
        
        M = length(folds)
      }else{
        if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.")
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    }else{
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    
    error.mat = array(0, dim = c(ncomp, nmthdd, M))
    
    # in case the test set only includes one sample, it is better to advise the user to
    # perform loocv
    stop.user = FALSE
    # set up a progress bar
    if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
    
    for (i in 1:M) {
      if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
      
      #set up leave out samples.
      omit = folds[[i]]
      
      # see below, we stop the user if there is only one sample drawn on the test set using MFold
      if(length(omit) == 1) stop.user = TRUE
      
      # the training set is NOT scaled
      X.train = X[-omit, ]
      Y.train = Y[-omit]      
      X.test = matrix(X[omit, ], nrow = length(omit))
      
      
      plsda.res = plsda(X = X.train, Y = Y.train, ncomp = ncomp, 
                        max.iter = max.iter, tol = tol)
      
      if (!is.null(plsda.res$nzv$Position)) X.test = X.test[, -plsda.res$nzv$Position]
      Y.predict = predict(plsda.res, X.test, method = method.predict)$class
      error.mat[, , i] = sapply(Y.predict, error.fun, y = as.numeric(Y[omit]))
    }
    
    # warn the user that at least test set had a length of 1
    if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
    
    if (progressBar == TRUE) cat('\n')
    
    #-- compute the error --#
    error.rate = apply(error.mat, 1:2, mean)
    rownames(error.rate) = paste('ncomp', 1:ncomp, sep = " ")
    colnames(error.rate) = names(Y.predict)
    
    result = list()
    result$error.rate = error.rate
    
    method = "plsda.mthd"
    class(result) = c("perf", method)
    #updated outputs
    return(invisible(result))
  }

