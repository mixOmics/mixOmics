# ---------------------
# this function will be included in mixOmics
# at the moment not fully tested
# date: 08/05/2015
# author: BG, KA
# --------------------

library(ellipse)

# note KA:
# I have added this argument:
#level  
#The confidence level of a pairwise confidence region. The default is 0.95, for a 95% region. This is used to control the size of the ellipse being plotted. A vector of levels may be used.

plotIndiv2 <-function (object, comp = 1:2, ind.names = TRUE, rep.space = "X-variate", 
                        X.label = NULL, Y.label = NULL, col = "black", cex = 1, pch = 1, 
                        abline.line = FALSE, plot.ellipse = TRUE, ellipse.ellipse.level = 0.95, ...) 
{
  if (length(comp) != 2) 
    stop("'comp' must be a numeric vector of length 2.")
  if (!is.numeric(comp)) 
    stop("invalid vector for 'comp'.")
  if (length(comp) == 1) 
    stop("Need at least 2 components to plot the graph")
  if (any(comp > object$ncomp)) 
    stop("the elements of 'comp' must be smaller or equal than ", 
         object$ncomp, ".")
  if (is.logical(ind.names)) {
    if (isTRUE(ind.names)) 
      ind.names = object$names$indiv
  }
  if (length(ind.names) > 1) {
    if (length(ind.names) != nrow(object$X)) 
      stop("'ind.names' must be a character vector of length ", 
           nrow(object$X), " or a boolean atomic vector.")
  }
  
  
  # KA added: need to set up the colors as a list of vectors
  # the first element indicates the colors of the samples, can be one value or several (length = # samples)
  # the second element should indicate the color of the ellipses and should be the length of the number of groups - indicated by ind.mat
  if(length(col) == 1){
    col = list(col, rep(col, ncol(object$ind.mat)))
  }
  
  comp1 = round(comp[1])
  comp2 = round(comp[2])
  rep.space = match.arg(rep.space, c("XY-variate", "X-variate", 
                                     "Y-variate"))
  match.user = names(match.call())
  match.function = c("object", "comp", "rep.space", "X.label", 
                     "Y.label", "col", "cex", "pch", "abline.line", "plot.ellipse")
  
#   if (length(setdiff(match.user[-1], match.function)) != 0) 
#     warning("Some of the input arguments do not match the function arguments, see ?plotIndiv")
  
  
  if (rep.space == "X-variate") {
    x = object$variates$X[, comp1]
    y = object$variates$X[, comp2]
    if (is.null(X.label)) 
      X.label = paste("X-variate", comp1)
    if (is.null(Y.label)) 
      Y.label = paste("X-variate", comp2)
  }
  if (rep.space == "Y-variate") {
    x = object$variates$Y[, comp1]
    y = object$variates$Y[, comp2]
    if (is.null(X.label)) 
      X.label = paste("Y-variate", comp1)
    if (is.null(Y.label)) 
      Y.label = paste("Y-variate", comp2)
  }
  if (rep.space == "XY-variate") {
    x = (object$variates$X[, comp1] + object$variates$Y[, comp1])/2
    y = (object$variates$X[, comp2] + object$variates$Y[, comp2])/2
    if (is.null(X.label)) 
      X.label = paste("X-variate", comp1)
    if (is.null(Y.label)) 
      Y.label = paste("Y-variate", comp2)
  }
  
  ind.gp = matrice = cdg = variance = list()
  ind.gp = lapply(1:ncol(object$ind.mat), function(x){which(object$ind.mat[, x]==1)})
  matrice = lapply(ind.gp, function(z){matrix(c(x[z], y[z]), ncol = 2)})
  cdg = lapply(matrice, colMeans)
  variance = lapply(matrice, var) 
  coord.ellipse = lapply(1:ncol(object$ind.mat), function(x){ellipse(variance[[x]], centre = cdg[[x]], ellipse.level = 0.95)})
  #variance = lapply(variance, function(z){z/length(x)})
  max.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, max)})
  min.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, min)})
  
  if (length(ind.names) > 1) {
    plot(x, y, type = "n", xlab = X.label, ylab = Y.label, 
         xlim = c(min(x, min.ellipse[1, ]), max(x, max.ellipse[1, ])), ylim = c(min(y, min.ellipse[2, ]), max(y, max.ellipse[2, ])))
    text(x, y, ind.names, col = col[[1]], cex = cex, ...)  
    
    if (plot.ellipse == TRUE){
      for (i in 1 : length(ind.gp)){
        points(coord.ellipse[[i]], type = "l", col = col[[2]][i], cex = 50)
      }
    }
    
    if (abline.line) 
      abline(v = 0, h = 0, lty = 2)
  }
  else {
    if (isTRUE(ind.names)) {    
      plot(x, y, type = "n", xlab = X.label, ylab = Y.label,
           xlim = c(min(x, min.ellipse[1, ]), max(x, max.ellipse[1, ])), ylim = c(min(y, min.ellipse[2, ]), max(y, max.ellipse[2, ])))
      text(x, y, ind.names, col = col[[1]], cex = cex, ...)
      
      if (plot.ellipse == TRUE){
        for (i in 1: length(ind.gp)){
          points(coord.ellipse[[i]], type = "l", col = col[[2]][i], cex = 50)
        }
      }
      
      if (abline.line) 
        abline(v = 0, h = 0, lty = 2)
    }
    else {      
      plot(x, y, xlab = X.label, ylab = Y.label, col = col[[1]], cex = cex, pch = pch,
           xlim = c(min(x, min.ellipse[1, ]), max(x, max.ellipse[1, ])), ylim = c(min(y, min.ellipse[2, ]), max(y, max.ellipse[2, ])))
      
      if (plot.ellipse == TRUE){
        for (i in 1: length(ind.gp)){
          points(coord.ellipse[[i]], type = "l", col = col[[2]][i], cex = 50)
        }
      }
      
      if (abline.line) 
        abline(v = 0, h = 0, lty = 2)
    }
  }
}
