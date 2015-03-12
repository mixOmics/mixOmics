tune.multilevel <- function (X, Y = NULL, # *BG* cond = NULL, sample, replace cond and sample arguments by a design matrix
                             design = NULL, #  *BG* add design
                             ncomp = 1, test.keepX = c(5, 10, 15), test.keepY = NULL, 
                             already.tested.X = NULL, already.tested.Y = NULL, 
                             method = NULL, dist, validation = c("Mfold", "loo"), folds = 10) {
  
  X = as.matrix(X)
  if (length(dim(X)) != 2 || !is.numeric(X)) 
    stop("'X' must be a numeric matrix.")
  
# *BG* Start: check design instead of cond and sample
#   if (is.null(cond)) 
#     stop("Vector cond is missing", call. = FALSE)
#   
#   if (!is.null(dim(sample))) 
#     stop("The sample  vector should indicate the repeated measurements")
# 
# if (length(sample) != nrow(X)) 
#   stop("X and the vector sample should have the same number of subjects")
# 
# if (is.factor(sample)) {
#   sample = as.numeric(sample)
#   warning("the vector sample was converted into a numeric vector", call. = FALSE)
# }
# 
# if (length(summary(as.factor(sample))) == nrow(X)) 
#   stop("Check that the vector sample reflects the repeated measurements")
# 
# if (!any(names(summary(as.factor(sample))) == "1")) {
#   cat("The vector sample includes the values: ", as.vector(names(summary(as.factor(sample)))), "\n")
#   stop("sample vector", call. = FALSE)
# }
  
  if(missing(design)) 
    stop('You forgot to input the design matrix.', call. = FALSE)
  
  if(is.null(design)) 
    stop("the 'design' matrix is missing.", call. = FALSE)

  if ((nrow(design) != nrow(X))) 
    stop("unequal number of rows in 'X' and 'design'.", call. = FALSE)
  
  if (ncol(design) < 2) 
    stop("'design' must be a matrix or data frame with at least 2 columns.", call. = FALSE)

  design = as.data.frame(design)

  if (length(design[, 1]) != nrow(X)) 
    stop("X and the vector sample should have the same number of subjects")
  
  if (is.factor(design[, 1])) {
    design[, 1] = as.numeric(design[, 1])
    warning("the vector sample was converted into a numeric vector", call. = FALSE)
  }

  if (length(summary(as.factor(design[, 1]))) == nrow(X)) 
    stop("Check that the vector sample reflects the repeated measurements")
  
  if (!any(names(summary(as.factor(design[, 1]))) == "1")) {
    cat("The vector sample includes the values: ", as.vector(names(summary(as.factor(sample)))), "\n")
    stop("sample vector", call. = FALSE)
  }

  # *BG* End: check design instead of cond and sample

  if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0) 
    stop("invalid number of variates, 'ncomp'.")
  
  if (is.null(method)) 
    stop("Input method missing, should be set to splsda or spls", call. = FALSE)
  
  if ((method == "splsda") && (!is.null(already.tested.Y))) 
    warning("Unecessary already.tested.Y parameter")
  
  if (method == "splsda"){
    # if (is.null(dim(cond))) {# *BG* Replace cond by design
    if (dim(design)[2] == 2) {
      result = tune.splsdalevel1(X = X, # *BG* remove cond = cond, sample, 
                                 design = design, #  *BG* add design
                                 ncomp = ncomp, test.keepX = test.keepX, dist = dist, 
                                 already.tested.X = already.tested.X, validation, folds)
    } else {
      result = tune.splsdalevel2(X = X, # *BG* remove cond = cond, sample, 
                                 design = design, #  *BG* add design
                                 ncomp = ncomp, test.keepX = test.keepX, already.tested.X = already.tested.X)
    }
  } else {
    result = tune.splslevel(X = X, Y = Y, # *BG* remove cond = cond, sample = sample, 
                            design = design, #  *BG* add design
                            ncomp = ncomp, test.keepX = test.keepX, test.keepY = test.keepY, 
                            already.tested.X = already.tested.X, already.tested.Y = already.tested.Y)
  }
  return(result)
}
