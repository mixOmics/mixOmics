tune.splsdalevel2 <- function (X, #cond, sample,  # *BG* add design argument
                               design,
                               ncomp = 1, test.keepX = c(5, 10, 15), already.tested.X = NULL) {
  
  cat("For a two-factor analysis, the tuning criterion is based on the maximisation of the correlation between the components on the whole data set", "\n")
  
  if (!is.null(already.tested.X)) 
    cat("Number of variables selected on the first ", ncomp - 1, "component(s) was ", already.tested.X)
  
  if (length(already.tested.X) != (ncomp - 1)) 
    stop("The number of already tested parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)
  
  if ((!is.null(already.tested.X)) && (!is.numeric(already.tested.X))) 
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)
  
  # *BG* Start: check condition design matrix instead of cond
#   if ((is.matrix(cond)) && ncol(cond) > 2) 
#     stop("'cond' must be a matrix with 2 columns for the 2 factor analysis.")
#   
#   if (ncol(cond) == 2) {
#     if (!is.factor(cond[, 1])) {
#       warning("First cond response was set as a factor", call. = FALSE)
#     }
#     if (!is.factor(cond[, 2])) {
#       warning("Second cond response was set as a factor", call. = FALSE)
#     }
#     cond1 = as.factor(cond[, 1])
#     cond2 = as.factor(cond[, 2])
#   }

  if (is.data.frame(design[, -1]) && ncol(design[, -1]) > 2) 
    stop("'cond' must be a matrix with 2 columns for the 2 factor analysis.")

  if (ncol(design[, -1]) == 2){
    if (!is.factor(design[, 2])) {
      design[, 2] = as.factor(design[, 2])
      warning("First cond response was set as a factor", call. = FALSE)
    }
    if (!is.factor(design[, 3])) {
      design[, 3] = as.factor(design[, 3])
      warning("Second cond response was set as a factor", call. = FALSE)
    }
  }    
  # *BG* End: check condition design matrix instead of cond
  
  # cond.fact = as.factor(paste(cond1, cond2, sep = ".")) # *BG* use design matrix
  cond.fact = as.factor(paste(design[, 2], design[, 3], sep = ".")) # *BG* use design matrix

  # Xw <- Split.variation.two.level(X, cond1, cond2, sample)$Xw # *BG* use withinVariation instead
  Xw <- suppressMessages(withinVariation(X = X, design = design)) # *BG* use withinVariation instead
  cor.value = vector(length = length(test.keepX))
  names(cor.value) = paste("var", test.keepX, sep = "")
  
  for (i in 1:length(test.keepX)) {
    if (ncomp == 1) {
      spls.train = splsda(Xw, cond.fact, ncomp = ncomp, keepX = test.keepX[i])
    } else {
      spls.train = splsda(Xw, cond.fact, ncomp = ncomp, keepX = c(already.tested.X, test.keepX[i]))
    }
    
    # Note: this is performed on the full data set
    # (could be done with resampling (bootstrap) (option 1) and/or prediction (option 2))
    
    # *BG* Start: Carry out the computation of the correlation on the centered/scaled data
#     for (h in 1:ncomp) {
#       if (h == 1) {
#         X.deflated = Xw
#         cond.deflated = unmap(as.numeric(cond.fact))
#       } else {
#         X.deflated = X.deflated - spls.train$variates$X[, h - 1] %*% t(spls.train$mat.c[, h - 1])
#         cond.deflated = cond.deflated - spls.train$variates$Y[, h - 1] %*% t(spls.train$mat.d[, h - 1])
#       }
#       cor.value[i] = cor(as.matrix(X.deflated) %*% spls.train$loadings$X[, h], as.matrix(cond.deflated) %*% spls.train$loadings$Y[, h])
#     }
    cor.value[i] = cor(spls.train$variates$X[, ncomp], spls.train$variates$Y[, ncomp])
    # *BG* End: Carry out the computation of the correlation on the centered/scaled data
  }
  return(list(cor.value = cor.value))
}