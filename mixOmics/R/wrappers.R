wrapper.sgcca <- function (blocks, design = NULL, penalty = NULL, scheme = "centroid", 
                           scale = TRUE, bias = FALSE, max.iter = 1000,
                           tol = .Machine$double.eps, verbose = FALSE, near.zero.var = FALSE,
                           keep.blocks = NULL, ncomp = rep(2, length(blocks))) {
    
  result.sgcca = srgcca(blocks = blocks, indY = NULL, design = design, tau = rep(1, length(blocks)), 
                        ncomp = ncomp, scheme = scheme, scale = scale, bias = bias, init = "svd.sgcca", 
                        tol = tol, verbose = verbose, mode = "canonical", max.iter = max.iter, 
                        keep.blocks = keep.blocks, near.zero.var = near.zero.var, penalty = penalty)

  result.sgcca$class = "sgcca"
  class(result.sgcca) = "sgcca"
  return(invisible(result.sgcca))
}

wrapper.rgcca <- function (blocks, design = NULL, penalty = NULL, scheme = "centroid", 
                           scale = TRUE, bias = FALSE, max.iter = 1000,
                           tol = .Machine$double.eps, verbose = FALSE, near.zero.var = FALSE,
                           ncomp = rep(2, length(blocks)), tau = "optimal") {
  
  result.sgcca = srgcca(blocks = blocks, indY = NULL, design = design, tau = tau, 
                        ncomp = ncomp, scheme = scheme, scale = scale, bias = bias, init = "svd.rgcca", 
                        tol = tol, verbose = verbose, mode = "canonical", max.iter = max.iter, 
                        keep.blocks = NULL, near.zero.var = near.zero.var, penalty = NULL)
  
  result.sgcca$class = "rgcca"
  class(result.sgcca) = "rgcca"
  return(invisible(result.sgcca))
}

wrapper.sgccda <- function (blocks, Y, design = NULL, scheme = "centroid", 
                           scale = TRUE, bias = FALSE, max.iter = 1000,
                           tol = .Machine$double.eps, verbose = FALSE, near.zero.var = FALSE,
                           keep.blocks = NULL, ncomp = rep(2, length(blocks))) {

  #-- Define blocks
  if (!is.list(blocks))
    stop("'blocks' must be a list containing the data set.", call. = FALSE)  
  
  if (is.null(dim(Y))) {
    Y = as.factor(Y)
    ind.mat = data.frame(unmap(as.numeric(Y)))
    names(ind.mat) = levels(Y)    
  } else {
    stop("'Y' should be a factor or a class vector.")
  }
  
  #-- ncomp
  if (!is.vector(ncomp, mode = "double"))
    stop("'ncomp' must be a vector")
  
  if (length(ncomp) == length(blocks)){
    ncomp = c(ncomp, max(ncomp))
    message("'ncomp' has changed and extended to Y")
  }
  
  #-- Design
  if (is.null(design)) {
    design = 1 - diag(length(blocks) + 1)
  } else if (ncol(design) != nrow(design) || ncol(design) < length(blocks) || ncol(design) > (length(blocks) + 1) || any(!design %in% c(0,1))){
    stop('Invalid design matrix.')
  } else if (ncol(design) == length(blocks)){ 
    message('Design matrix has changed and extended to Y)')
    design = rbind(cbind(design, 1), 1)
    diag(design) = 0
  }
  
  #-- keep.blocks
  if (is.list(keep.blocks) & (length(keep.blocks) == length(blocks))){
    keep.blocks[[length(keep.blocks) + 1]] = rep(nlevels(Y), ncomp[length(blocks) + 1])
    message("'keep.blocks' matrix has changed and extended to Y)")
  }
  
  blocks[[length(blocks) + 1]] = ind.mat
  names(blocks)[length(blocks)] = "Y"
    
  result.sgcca = srgcca(blocks = blocks, indY = length(blocks), design = design, tau = rep(1, length(blocks)), 
                        ncomp = ncomp, scheme = scheme, scale = scale, bias = bias, init = "svd.da", 
                        tol = tol, verbose = verbose, mode = "regression", max.iter = max.iter,
                        keep.blocks = keep.blocks, near.zero.var = near.zero.var, penalty = NULL)
  
  result.sgcca$Y = Y; result.sgcca$ind.mat = ind.mat;
  result.sgcca$class = "sgccda"
  class(result.sgcca) = "sgccda"
  return(invisible(result.sgcca))
}