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


# --------------------------
# declare the S3 function:
# -------------------------
perf <- function(object, ...) UseMethod("perf")


#------------------------------------------------------#
#-- Includes perf for PLS, sPLS, PLS-DA and sPLS-DA --#
#------------------------------------------------------#
# ---------------------------------------------------
# perf for pls object
# ---------------------------------------------------

perf.pls <-
  function(object,
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE,
           near.zero.var = FALSE,
           ...)
  {
    

    #------------------#
    #-- check entries --#
    
    #-- check pls mode
    if(object$mode == 'canonical')
      stop("PLS 'mode' should be set to 'regression', 'invariant' or 'classic'.", call. = FALSE)
    
    #-- validation
    choices = c("Mfold", "loo")
    validation = choices[pmatch(validation, choices)]
    
    if (any(is.na(validation)) || length(validation) > 1)
      stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- progressBar
    if (!is.logical(progressBar))
      stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    #-- cross-validation approach  ---------------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- initialising arguments --#
    # these are the centered and scaled matrices output from pls
    X = object$X
    Y = object$Y
    
    tol = object$tol
    max.iter = object$max.iter
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()
    
    if (any(is.na(X)) || any(is.na(Y))) 
      stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
    
    # -------------------------------------
    # added: first check for near zero var on the whole data set
    nzv = nearZeroVar(X)
    if (length(nzv$Position > 0))
    {
      warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
      X = X[, -nzv$Position, drop=FALSE]
      
      if(ncol(X)==0) {stop("No more predictors after Near Zero Var has been applied!")}
      
    }
    # and then we start from the X data set with the nzv removed
    

    #-- M fold or loo cross validation --#
    #-- define the folds --#
    if (validation == "Mfold") {
        if (is.list(folds)) {
            
            if (length(folds) < 2 || length(folds) > n)
              stop("Invalid number of folds.", call. = FALSE)
            
            if (length(unlist(folds)) != n)
              stop("Invalid folds. The total number of samples in folds must be equal to ", n, ".", call. = FALSE)
            
            if (length(unique(unlist(folds))) != n)
              stop("Invalid folds. Repeated samples in folds.", call. = FALSE)
            
            M = length(folds)
        } else {
            if (is.null(folds) || !is.finite(folds) || folds < 2 || folds > n) {
                stop("Invalid number of folds.", call. = FALSE)
            } else {
                M = round(folds)
                folds = split(sample(1:n), rep(1:M, length = n))
            }
        }
    } else {
        folds = split(1:n, rep(1:n, length = n))
        M = n
    }
    
    
    ## add Q2
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = q)})
    RSS.indiv[[1]] = X
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)

    
    #KA: all criteria included in the computation

    press.mat = Ypred = MSEP.mat = array(0, c(n, q, ncomp))
    MSEP = R2 = matrix(0, nrow = q, ncol = ncomp)
    
    # set up dimnames
    rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
    dimnames(press.mat)[[2]] = colnames(Y)
    
    # in case the test set only includes one sample, it is better to advise the user to perform loocv
    stop.user = FALSE
    if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
    

    
    for (i in 1:M)
    {
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
        mode = mode, max.iter = max.iter, tol = tol, near.zero.var = near.zero.var)
        
        if (!is.null( pls.res$nzv$Position)) X.test = X.test[, - pls.res$nzv$Position,drop=FALSE]
        
        # in the predict function, X.test is already normalised w.r.t to training set X.train, so no need to do it here
        Y.hat = predict( pls.res, X.test)$predict #used for Ypred and MSEP.mat
        
        
        for (h in 1:ncomp)
        {
            Ypred[omit, , h] = Y.hat[, , h]
            # compute the MSEP
            MSEP.mat[omit, , h] = (Y.test - Y.hat[, , h])^2

            # update 5.0-4, another way to calculate Y.hat is used for press.mat
            a.cv=pls.res$loadings$X[,h]
            d.cv = t(Y.train) %*% (X.train %*% a.cv) / norm((X.train %*% a.cv), type = "2")^2
            Y.hat.cv = (X.test %*% a.cv) %*% t(d.cv)
            press.mat[omit, , h] = (Y.test - Y.hat.cv)^2

        } # end h
    } #end i (cross validation)
    
    
    # these criteria are computed across all folds.
    for (h in 1:ncomp)
    {
        MSEP[, h] = apply(as.matrix(MSEP.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[, , h], use = "pairwise")))^2
        
        # RSS calculated on the whole data, nothing to do with the folds
        d = object$mat.d[, h]
        tt = object$variates$X[, h]
        RSS.indiv[[h + 1]] = Y - tt %*% t(d)
        RSS[h + 1, ] = colSums((Y - tt %*% t(d))^2)
        
        # on en profite pour calculer le PRESS par composante et le Q2
        if(q>1){
            PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
            PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
    }
    

    
    Q2.total = matrix(1 - rowSums(PRESS.inside) / rowSums(RSS[-(ncomp+1), , drop = FALSE]), nrow = 1, ncol = ncomp,
                dimnames = list("Q2.total", paste0(1:ncomp, " comp")))
    
    # set up dimnames
    colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
    rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = object$names$Y
    
    if (progressBar == TRUE) cat('\n')


    res = list()
    res$MSEP = MSEP
    res$R2 = R2
    res$Q2 = t(Q2.inside)
   	res$Q2.total =  t(Q2.total)
    res$RSS = RSS
    res$PRESS.inside = PRESS.inside
    res$press.mat = press.mat
    res$RSS.indiv = RSS.indiv
    res$nzvX = nzv$Position

    class(res) = c("perf", "pls.mthd")
    return(invisible(res))
  }

# ===========================================================================================
#---------------------------------------------------
# perf for spls object
perf.spls <-
  function(object,
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE,
           near.zero.var = FALSE,
           ...)
  {
    
    #------------------#
    #-- check entries --#
    
    #-- check spls mode
    if(object$mode == 'canonical')
    stop("SPLS 'mode' should be set to 'regression'.", call. = FALSE)
    
    #-- validation
    choices = c("Mfold", "loo")
    validation = choices[pmatch(validation, choices)]
    
    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    # initialize new objects:= to record feature stability
    featuresX  = featuresY =  list()
    for(k in 1:ncomp){
        featuresX[[k]] = featuresY[[k]] = NA
    }
    #-- cross-validation approach  ---------------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- initialising arguments --#
    # these are the centered and scaled matrices output from pls
    X = object$X
    Y = object$Y

    tol = object$tol
    max.iter = object$max.iter
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    res = list()

    if (any(is.na(X)) || any(is.na(Y)))
    stop("missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.", call. = FALSE)
    
    
    #-- tells which variables are selected in X and in Y --#
    keepX = object$keepX   
    keepY = object$keepY
        
    # -------------------------------------
    # added: first check for near zero var on the whole data set
    nzv = nearZeroVar(X)
    if (length(nzv$Position > 0))
    {
      warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
      X = X[, -nzv$Position, drop=TRUE]
      
      if(ncol(X)==0) {stop("No more predictors after Near Zero Var has been applied!")}
      
    }
    # and then we start from the X data set with the nzv removed
    
    
    
    #-- M fold or loo cross validation --#
    #- define the folds
    if (validation == "Mfold") {
      if (is.list(folds)) {
        
        if (length(folds) < 2 || length(folds) > n)
          stop("Invalid number of folds.", call. = FALSE)
        
        if (length(unlist(folds)) != n)
          stop("Invalid folds. The total number of samples in folds must be equal to ", 
               n, ".", call. = FALSE)
        
        if (length(unique(unlist(folds))) != n)
          stop("Invalid folds. Repeated samples in folds.", call. = FALSE)
        
        M = length(folds)
      } else {
        if (is.null(folds) || !is.finite(folds) || folds < 2 || folds > n)
          stop("Invalid number of folds.", call. = FALSE)
        else {
          M = round(folds)
          folds = split(sample(1:n), rep(1:M, length = n)) 
        }
      }
    } else { 
      folds = split(1:n, rep(1:n, length = n)) 
      M = n
    }
    

    ## add Q2
    RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
    RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = q)})
    RSS.indiv[[1]] = X
    PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
    
    
    #KA: all criteria included in the computation
    
    press.mat = Ypred = MSEP.mat = array(0, c(n, q, ncomp))
    MSEP = R2 = matrix(0, nrow = q, ncol = ncomp)
    
    # set up dimnames
    rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
    dimnames(press.mat)[[2]] = colnames(Y)
    
    # in case the test set only includes one sample, it is better to advise the user to perform loocv
    stop.user = FALSE
    if (progressBar == TRUE) pb <- txtProgressBar(style = 3)
    
    
    
    for (i in 1:M)
    {
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

        #-- spls --#
        # added the near.zero.var option
        spls.res = spls(X.train, Y.train, ncomp=ncomp, mode=mode,
            max.iter=max.iter, tol=tol, keepX=keepX, keepY=keepY, near.zero.var = near.zero.var)
        
        # added: record selected features in each set
        for(k in 1:ncomp)
        {
          featuresX[[k]] = c(unlist(featuresX[[k]]), select.var(spls.res, comp = k)$name.X)
          featuresY[[k]] = c(unlist(featuresY[[k]]), select.var(spls.res, comp = k)$name.Y)
        }
        
        if (!is.null( spls.res$nzv$Position)) X.test = X.test[, - spls.res$nzv$Position,drop=FALSE]
        
        # in the predict function, X.test is already normalised w.r.t to training set X.train, so no need to do it here
        Y.hat = predict( spls.res, X.test)$predict #used for Ypred and MSEP.mat
        
        
        for (h in 1:ncomp)
        {
            Ypred[omit, , h] = Y.hat[, , h]
            # compute the MSEP
            MSEP.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
            
            # update 5.0-4, another way to calculate Y.hat is used for press.mat
            a.cv=spls.res$loadings$X[,h]
            d.cv = t(Y.train) %*% (X.train %*% a.cv) / norm((X.train %*% a.cv), type = "2")^2
            Y.hat.cv = (X.test %*% a.cv) %*% t(d.cv)
            press.mat[omit, , h] = (Y.test - Y.hat.cv)^2
            
        } # end h
    } #end i (cross validation)
    
    
    # these criteria are computed across all folds.
    for (h in 1:ncomp)
    {
        MSEP[, h] = apply(as.matrix(MSEP.mat[, , h]), 2, mean, na.rm = TRUE)
        R2[, h] = (diag(cor(Y, Ypred[, , h], use = "pairwise")))^2
        
        # RSS calculated on the whole data, nothing to do with the folds
        d = object$mat.d[, h]
        tt = object$variates$X[, h]
        RSS.indiv[[h + 1]] = Y - tt %*% t(d)
        RSS[h + 1, ] = colSums((Y - tt %*% t(d))^2)
        
        # on en profite pour calculer le PRESS par composante et le Q2
        if(q>1){
            PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
        }else{
            PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
        }
        
        Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
        
    }
    
    Q2.total = matrix(1 - rowSums(PRESS.inside) / rowSums(RSS[-(ncomp+1), , drop = FALSE]), nrow = 1, ncol = ncomp,
    dimnames = list("Q2.total", paste0(1:ncomp, " comp")))
    
    # set up dimnames
    colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
    rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = object$names$Y
    
    if (progressBar == TRUE) cat('\n')


    #---- extract stability of features -----#
    list.features.X = features.final.X = list()
    list.features.Y = features.final.Y = list()
    
    for(k in 1:ncomp){            
        #remove the NA value that was added for initialisation
        remove.naX = which(is.na(featuresX[[k]]))
        remove.naY = which(is.na(featuresY[[k]]))
        # then summarise as a factor and output the percentage of appearance
        list.features.X[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
        list.features.Y[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
      
      #-- extract features selected from the full model
      features.final.X[[k]] = row.names(object$loadings$X)[object$loadings$X[, 1, drop = FALSE] != 0]
      features.final.Y[[k]] = row.names(object$loadings$Y)[object$loadings$Y[, 1, drop = FALSE] != 0]
    }
    
    names(features.final.X) = names(list.features.X) = names(features.final.Y) = names(list.features.Y) = paste('comp', 1:ncomp)
    

    
    res = list()
    res$MSEP = MSEP
    res$R2 = R2
    res$Q2 = t(Q2.inside)
   	res$Q2.total =  t(Q2.total)
    res$RSS = RSS
    res$PRESS.inside = PRESS.inside
    res$press.mat = press.mat
    res$RSS.indiv = RSS.indiv
    
    # features
    res$features$stable.X = list.features.X
    res$features$stable.Y = list.features.Y
    res$features$final.X = features.final.X
    res$features$final.Y = features.final.Y
    res$nzvX = nzv$Position
    
    class(res) = c("perf", "spls.mthd")
    return(invisible(res))
}

# ---------------------------------------------------
# perf for plsda object
# ---------------------------------------------------
perf.plsda <-
  function(object,
           method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
           validation = c("Mfold", "loo"),
           folds = 10,
           progressBar = TRUE, 
           near.zero.var = FALSE,...) 
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
    
    
    # -------------------------------------
    # added: first check for near zero var on the whole data set
    nzv = nearZeroVar(X)
    if (length(nzv$Position > 0))
    {
      warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
      X = X[, -nzv$Position, drop=TRUE]
      
      if(ncol(X)==0)
      {
        stop("No more predictors after Near Zero Var has been applied!")
      }
      
    }
    # and then we start from the X data set with the nzv removed
    
    
    
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
                        max.iter = max.iter, tol = tol, near.zero.var = near.zero.var)
      
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
    # added
    result$nzvX = nzv$Position
    
    
    method = "plsda.mthd"
    class(result) = c("perf", method)
    #updated outputs
    return(invisible(result))
  }

# ---------------------------------------------------
# perf for splsda object
# ---------------------------------------------------
perf.splsda <- function(object,                                               
                        method.predict = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),     
                        validation = c("Mfold", "loo"),                                          
                        folds = 10,                                                              
                        progressBar = TRUE, 
                        near.zero.var = FALSE,...)                                                        
{
  
  #-- initialising arguments --#
  # these data are the centered and scaled X output or the unmapped(Y) scaled and centered
  X = object$X
  level.Y = object$names$Y  #to make sure the levels are ordered
  Y = object$ind.mat
  Y = map(Y)
  Y = factor(Y,labels = level.Y)
  ncomp = object$ncomp
  n = nrow(X)
  keepX = object$keepX  
  
  tol = object$tol
  max.iter = object$max.iter
  
  # initialize new objects:
  features <- list()
  for(k in 1:ncomp){
    features[[k]] = NA
  }
  
  method.predict = match.arg(method.predict, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
  if (any(method.predict == "all")) nmthdd = 3 
  else nmthdd = length(method.predict)  
  
  
  # -------------------------------------
  # added: first check for near zero var on the whole data set
  nzv = nearZeroVar(X)
  if (length(nzv$Position > 0))
  {
    warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
    X = X[, -nzv$Position, drop=TRUE]
    
    if(ncol(X)==0)
    {
      stop("No more predictors after Near Zero Var has been applied!")
    }
    
  }
  # and then we start from the X data set with the nzv removed
  
  
  
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
    }
    else {
      if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
        stop("Invalid number of folds.")
      else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n)) 
      }
    }
  } 
  else { 
    folds = split(1:n, rep(1:n, length = n)) 
    M = n
  }
  
  
  error.mat = array(0, dim = c(ncomp, nmthdd, M))
  
  # in case the test set only includes one sample, it is better to advise the user to perform loocv
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
    
    spls.res = splsda(X.train, Y.train, ncomp, max.iter, tol, keepX=keepX, near.zero.var = near.zero.var)  
    # added: record selected features
    for(k in 1:ncomp){
      features[[k]] = c(unlist(features[[k]]), select.var(spls.res, comp = k)$name)
    }
    
    if (!is.null(spls.res$nzv$Position)) X.test = X.test[, -spls.res$nzv$Position]
    Y.predict = predict(spls.res, X.test, method = method.predict)$class
    error.mat[, , i] = sapply(Y.predict, error.fun, y = as.numeric(Y[omit]))
    
    
  } # end loop on i
  
  # warn the user that at least test set had a length of 1
  if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
  
  if (progressBar == TRUE) cat('\n')
  
  #-- compute the error --#
  res = apply(error.mat, 1:2, mean)
  
  rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
  colnames(res) = names(Y.predict)
  
  # ---- extract stability of features ----- # NEW
  list.features = list()
  for(k in 1:ncomp){
    #remove the NA value that was added for initialisation
    remove.na = which(is.na(features[[k]]))
    # then summarise as a factor and output the percentage of appearance
    list.features[[k]] = sort(summary(as.factor(features[[k]][-remove.na]))/M, decreasing = TRUE)
  }
  
  # extract features selected from the full model ---------
  features.final = list()
  for(k in 1:ncomp){
    features.final[[k]] = select.var(object, comp = k)$value
  }
  
  names(features.final)  = names(list.features) = paste('comp', 1:ncomp)
  
  result = list()
  result$error.rate = res
  result$features$stable = list.features
  result$features$final = features.final
  
  # added
  result$nzvX = nzv$Position
  
  
  method = "splsda.mthd"
  #result$meth = "splsda.mthd"
  class(result) = c("perf", method)
  #updated outputs
  return(invisible(result))
}


