tune.splslevel <- function (X, Y, # cond = NULL, sample = NULL, # *BG* set up design matrix instead cond and sample
                            design, # *BG* set up design matrix instead cond and sample
                            ncomp = NULL, 
                            test.keepX = rep(ncol(X), ncomp), test.keepY = rep(ncol(Y), ncomp), 
                            already.tested.X = NULL, already.tested.Y = NULL) {
  
  cat("For a multilevel spls analysis, the tuning criterion is based on the maximisation of the correlation between the components from both data sets", "\n")
  
  Y = as.matrix(Y)
  if (length(dim(Y)) != 2 || !is.numeric(Y)) 
    stop("'Y' must be a numeric matrix.")
    
  if (!is.null(already.tested.X)) 
    cat("Number of X variables selected on the first ", ncomp - 1, "component(s) was ", already.tested.X, "\n")
  
  if (!is.null(already.tested.Y)) 
    cat("Number of Y variables selected on the first ", ncomp - 1, "component(s) was ", already.tested.Y, "\n")
  
  if ((!is.null(already.tested.X)) && is.null(already.tested.Y)) 
    stop("Input already.tested.Y is missing")
  
  if ((!is.null(already.tested.Y)) && is.null(already.tested.X)) 
    stop("Input already.tested.X is missing")
  
  if (length(already.tested.X) != (ncomp - 1)) 
    stop("The number of already.tested.X parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)
  
  if (length(already.tested.Y) != (ncomp - 1)) 
    stop("The number of already.tested.Y parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)
  
  if ((!is.null(already.tested.X)) && (!is.numeric(already.tested.X))) 
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)
  
  if ((!is.null(already.tested.Y)) && (!is.numeric(already.tested.Y))) 
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)
  
  # *BG* Start: replace Split.variation.one.level by withinVariation
  # Xw <- Split.variation.one.level(X, Y = cond, sample)$Xw 
  # Yw <- Split.variation.one.level(X = Y, Y = cond, sample)$Xw
  Xw <- suppressMessages(withinVariation(X = X, design = design))
  Yw <- suppressMessages(withinVariation(X = Y, design = design))
  # *BG* Start: replace Split.variation.one.level by withinVariation
  
  cor.value = matrix(nrow = length(test.keepX), ncol = length(test.keepY))
  rownames(cor.value) = paste("varX", test.keepX, sep = "")
  colnames(cor.value) = paste("varY", test.keepY, sep = "")
  mode = "canonical" # *BG* !!! Hard coding !!! should be an argument?
  
  for (i in 1:length(test.keepX)) {
    for (j in 1:length(test.keepY)) {
      if (ncomp == 1) {
        spls.train = spls(Xw, Yw, ncomp = ncomp, 
                          keepX = test.keepX[i], 
                          keepY = test.keepY[j],
                          mode = mode)
      } else {
        spls.train = spls(Xw, Yw, ncomp = ncomp, 
                          keepX = c(already.tested.X, test.keepX[i]),
                          keepY = c(already.tested.Y, test.keepY[j]), 
                          mode = mode)
      }
      
      # *BG* Start: Carry out the computation of the correlation on the centered/scaled data
#       for (h in 1:ncomp) {
#         if (h == 1) {
#           X.deflated = scale(Xw)
#           Y.deflated = scale(Yw)
#         } else {
#           X.deflated = X.deflated - spls.train$variates$X[, h - 1] %*% t(spls.train$mat.c[, h - 1])
#           Y.deflated = Y.deflated - spls.train$variates$Y[, h - 1] %*% t(spls.train$mat.e[, h - 1])
#         }
#         cor.value[i, j] = cor(as.matrix(X.deflated) %*% spls.train$loadings$X[, h], as.matrix(Y.deflated) %*% spls.train$loadings$Y[, h])
#       }
      cor.value[i, j] = cor(spls.train$variates$X[, ncomp], spls.train$variates$Y[, ncomp])
      # *BG* End: Carry out the computation of the correlation on the centered/scaled data
    }
  }
  return(list(cor.value = cor.value))
}