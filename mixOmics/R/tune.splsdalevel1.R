tune.splsdalevel1 <- function (X, ncomp = 1, test.keepX, dist = NULL, # *BG* remove sample, cond,
                               design = NULL, # *BG* add design argument
                              already.tested.X = NULL, validation, folds) {
  
  if (validation == "Mfold") {
    cat(paste("For a one-factor analysis, the tuning criterion is based on", folds, "cross validation"), "\n")
  } else {
    cat(paste("For a one-factor analysis, the tuning criterion is based on leave-one-out cross validation"), "\n")
  }
  
  # *BG* Replace cond by design
#   if (!is.factor(cond)) {
#     cond = as.factor(cond)
#     warning("cond was set as a factor", call. = FALSE)
#   }
  if (!is.factor(design[, 2])) {
    design[, 2] = as.factor(design[, 2])
    warning("design[, 2] was set as a factor", call. = FALSE)
  }
  
  if (is.null(dist)) 
    stop("Input dist is missing")
  
  if (length(already.tested.X) != (ncomp - 1)) 
    stop("The number of already tested parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)
  
  if ((!is.null(already.tested.X)) && (!is.numeric(already.tested.X))) 
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)
  
  if (!is.null(already.tested.X)) 
    cat("Number of variables selected on the first ", ncomp - 1, "component(s) was ", already.tested.X)
  
  sample = design[, 1]# *BG* add sample feature
  n = length(unique(sample))
  
  if (validation == "Mfold") {
    if (is.list(folds)) {
      if (length(folds) < 2 | length(folds) > n) 
        stop("Invalid number of folds.")
      if (length(unique(unlist(folds))) != n) 
        stop("Invalid folds.")
      M = length(folds)
    } else {
      if (is.null(folds) || !is.numeric(folds) || folds < 
            2 || folds > n) 
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
  
  vect.error = vector(length = length(test.keepX))
  names(vect.error) = paste("var", test.keepX, sep = "")
  error.sw = matrix(nrow = M, ncol = length(test.keepX), data = 0)
  
  # calculate first the within variation matrix, which is ok because Xw is not dependent on each subject (so there is not bias here)
  Xw = suppressMessages(withinVariation(X = X, design = design)) # *BG* Estimation of the design matrix. Outside the loop because Xw is "subject specific"

  for (j in 1:M) {
    omit = which(sample %in% folds[[j]] == TRUE)
    # *BG* Start: use Xw matrix
    # X.train = X[-c(omit), ]
    # X.test = matrix(X[omit, ], nrow = length(omit))
    trainxw = Xw[-c(omit), ]
    xwtest = Xw[omit, , drop = FALSE]
    # *BG* End: use Xw matrix
    
    cond.train = design[-c(omit), 2] # *BG* define Y.train variable
    cond.test = design[c(omit), 2] # *BG* defin Y.test variable
    
    # *BG* Start: remove unnecessary part 
    # cond.train = cond[-c(omit)]
    # cond.test = cond[omit]
    # xwtest <- X.test - matrix(colMeans(X.test), ncol = ncol(X.test), nrow = length(omit), byrow = TRUE)
    # samplew <- sample[-omit]
    # X.train.decomp <- Split.variation.one.level(X.train, cond.train, sample = samplew)
    # trainxw <- X.train.decomp$Xw
    # *BG* End: remove part unnecessary
    
    remove.zero = nearZeroVar(trainxw)$Position
    if (length(remove.zero) != 0) {
      trainxw = trainxw[, -c(remove.zero)]
      xwtest = xwtest[, -c(remove.zero)]
    }
    
    for (i in 1:length(test.keepX)) {
      if (ncomp == 1) {
        result.sw <- splsda(trainxw, cond.train, keepX = test.keepX[i], ncomp = ncomp)
      }
      else {
        result.sw <- splsda(trainxw, cond.train, keepX = c(already.tested.X, test.keepX[i]), ncomp = ncomp)
      }
      test.predict.sw <- predict(result.sw, xwtest, dist = dist)
      # Prediction.sw <- levels(cond)[test.predict.sw$class[[dist]][, ncomp]] *BG* Replace cond by change
      Prediction.sw <- levels(design[, 2])[test.predict.sw$class[[dist]][, ncomp]] 
      error.sw[j, i] <- sum(as.character(cond.test) != Prediction.sw)
    }
  }
  # result <- apply(error.sw, 2, sum)/length(cond) # *BG* Replace cond by design
  result <- apply(error.sw, 2, sum)/length(design[, 2]) # *BG* Replace cond by design

  names(result) = paste("var", test.keepX, sep = "")
  return(list(error = result))
}