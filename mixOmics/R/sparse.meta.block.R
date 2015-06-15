#############################################################################################################
#
# Date: July 20, 2014
# Author: Amrit Singh
# sparse generalized canonical correlation discriminant analysis (sgccda)
#
# sgccda functions; calcPenalty(),, sgccda(), wrapper.sGCCA(), valid.sgccda(),
#                   predict.sgccda(),
#                   wrapper.sgccda(), valid.sgcca(),
#                   sgccda.model(),  calcPenalty(),
#                   sens.spec(), plotVar.sgccda()
#
#
# Functions modified from RGCCA R-library
#   sgcca(), sgccak()
#
# Functions acquired from RGCCA R-library
#   cov2(), initsvd(), crossprod(),
#   defl.select(),
#
#
#############################################################################################################

#############################################################################################################
# Functions modified from RGCCA
#############################################################################################################
# ----------------------------------------------------------------------------------------------------------
# sgcca - Runs sgcca() modified from RGCCA
#   inputs: A - list of datasets each with the same number of rows (samples)
#           design - design matrix
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ----------------------------------------------------------------------------------------------------------
sparse.meta.block = function (A, indY = NULL,  design = 1 - diag(length(A)),tau=NULL,#rep(1, length(A)),
                            ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,  bias = FALSE,
                            init = "svd.single", tol = 1e-06, verbose = FALSE,
                            mode = "canonical", max.iter = 500,study = NULL, keepA = NULL,
                            keepA.constraint = NULL, near.zero.var = FALSE) { # meta.hybrid.spls
  

  # A: list of matrices
  # indY: integer, pointer to one of the matrices of A
  # design: design matrix, links between matrices. Diagonal must be 0
  # ncomp: vector of ncomp, per matrix
  # scheme: a function "g", refer to the article (thanks Benoit)
  # scale: do you want to scale ? mean is done by default and cannot be changed (so far)
  # bias: scale the data with n or n-1
  # init: one of "svd" or "random", initialisation of the algorithm
  # tol: nobody cares about this
  # verbose: show the progress of the algorithm
  # mode: canonical, classic, invariant, regression
  # max.iter: nobody cares about this
  # study: factor for each matrix of A, must be a vector
  # keepA: keepX of spls for each matrix of A. must be a list. Each entry must be of the same length (max ncomp)
  # keepA.constraint: keepX.constraint, which variables are kept on the first num.comp-1 components. It is a list of characters
  # near.zero.var: do you want to remove variables with very small variance
  
  
  # add check so that this function can run on its own (which shouldn't happen now, but don't know the future)
  
  
  #if(is.null(indY) & is.null(tau))
  #stop("Either 'indY' or 'tau' is needed")
  
  #check=Check.entry.meta.block.spls(A, indY, design ,ncomp , scheme , scale ,  bias,
  #                                  init , tol , verbose,mode, sparse , max.iter,study , keepA, keepA.constraint)
  #A=check$A
  #ncomp=check$ncomp
  #study=check$study
  
  
  # at this stage keepA.constraint need to be character, to remove easily variables with near zero variance
  ### near.zero.var, remove the variables with very small variances
  if(near.zero.var == TRUE)
  {
    nzv.A = lapply(A,nearZeroVar)
    for(q in 1:length(A))
    {
      if (length(nzv.A[[q]]$Position) > 0)
      {
        names.remove.X=colnames(A[[q]])[nzv.A[[q]]$Position]
        A[[q]] = A[[q]][, -nzv.A[[q]]$Position,drop=FALSE]
        if (verbose)
        warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
        if(ncol(A[[q]]) == 0) {stop("No more variables in X")}
        
        # at this stage, keepA.constraint need to be numbers
        if(length(keepA.constraint[[q]])>0)
        {
          #remove the variables from keepA.constraint if removed by near.zero.var
          keepA.constraint[[q]]=match.keepX.constraint(names.remove.X,keepA.constraint[[q]])
          # replace character by numbers
          keepA.constraint[[q]]= lapply(keepA.constraint[[q]],function(x){match(x,colnames(A[[q]]))})
        }
        #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
        if(any(keepA[[q]]>ncol(A[[q]])))
        {
            ind=which(keepA[[q]]>ncol(A[[q]]))
            keepA[[q]][ind]=ncol(A[[q]])
        }
      }

    }
  }
  
  # keepA is updated to be of length now, the first entries correspond to the keepA.constraint if it was provided
  
  for(q in 1:length(A))
    keepA[[q]]=c(unlist(lapply(keepA.constraint[[q]],length)),keepA[[q]]) #of length ncomp, can contains 0

#print(keepA)

  # center the data per study, per matrix of A, scale if scale=TRUE, option bias
  mean_centered = lapply(A, function(x){mean_centering_per_study(x, study, scale, bias)})
  A = lapply(mean_centered, function(x){as.matrix(x$concat.data)})
  ni=lapply(mean_centered[[1]]$data.list.study.scale,nrow) #number of samples per study


    ### Start: Initialization parameters
  pjs <- sapply(A, NCOL); nb_ind <- NROW(A[[1]]);
  J <- length(A); R <- A; N <- max(ncomp) # R: residuals matrices, will be a list of length ncomp
  
  AVE_inner <- AVE_outer <- rep(NA, max(ncomp))
  defl.matrix <- AVE_X <- crit <- loadings.partial.A <- variates.partial.A <- tau.rgcca <- list()
  P <- loadings.A <- loadings.Astar <- c <- t <- b <- variates.A <- vector("list",J)
  for (k in 1:J) t[[k]] <- variates.A[[k]] <- matrix(NA, nb_ind, N)
  for (k in 1:J) P[[k]] <- loadings.A[[k]] <- loadings.Astar[[k]]<- matrix(NA, pjs[[k]], N)
  for (k in 1:J)
  {
      loadings.partial.A[[k]]=variates.partial.A[[k]]=vector("list",length=nlevels(study))
      for(m in 1:nlevels(study))
      {
          loadings.partial.A[[k]][[m]]=matrix(nrow=NCOL(A[[k]]),ncol=N)
          variates.partial.A[[k]][[m]]=matrix(nrow=ni[[m]],ncol=N)
      }
      #variates.partial.A[[k]]=matrix(nrow=nb_ind,ncol=N)
  }
  
  defl.matrix[[1]] <- A; ndefl <- ncomp - 1;  J2 <- J-1;
  
  if (is.vector(tau)){
    tau = matrix(rep(tau, N), nrow = N, ncol = length(tau), byrow = TRUE)
  }
  ### End: Initialization parameters


  iter=NULL
  for (n in 1:N)
  {
    #if (verbose)
    #  cat(paste0("Computation of the SGCCA block components #", n, " is under progress... \n"))
    meta.block.result = NULL
    
    ### Start: Estimation ai
    if (is.null(tau))
    {
      meta.block.result <- sparse.meta.block_iteration(R, design,study = study,
                        keepA.constraint=if (!is.null(keepA.constraint)) {lapply(keepA.constraint, function(x){unlist(x[n])})} else {NULL} ,
                        keepA = if (!is.null(keepA)) {lapply(keepA, function(x){x[n]})} else {NULL},indY = indY,
                        scheme = scheme, init = init, max.iter = max.iter, tol = tol,   verbose = verbose)
    } else {
      meta.block.result <- sparse.rgcca_iteration(R, design, tau = if (is.matrix(tau)){tau[n, ]} else {"optimal"}, scheme = scheme, init = init, tol = tol,
                        verbose = verbose, max.iter = max.iter,
                        keepA.constraint=if (!is.null(keepA.constraint)) {lapply(keepA.constraint, function(x){unlist(x[n])})} else {NULL} ,
                        keepA = if (!is.null(keepA)) {lapply(keepA, function(x){x[n]})} else {NULL})
    }
    
    ### End: Estimation ai
    
    if(is.null(tau))
    {
        #recording loadings.partials, $Ai$study[,ncomp]
        # recording variates.partials, $Ai[,ncomp]
        for(k in 1:J)
        {
            for(m in 1:nlevels(study))
            {
                loadings.partial.A[[k]][[m]][,n]=matrix(meta.block.result$loadings.partial.A.comp[[k]][[m]],ncol=1)
                variates.partial.A[[k]][[m]][,n]=matrix(meta.block.result$variates.partial.A.comp[[k]][[m]],ncol=1)
            }
            #variates.partial.A[[k]][,n]=matrix(unlist(meta.block.result$variates.partial.A.comp[[k]]),ncol=1)
        }
    }
    
    AVE_inner[n] <- meta.block.result$AVE_inner
    crit[[n]] <- meta.block.result$crit
    tau.rgcca[[n]] <- meta.block.result$tau
    
    for (k in 1:J) variates.A[[k]][, n] <- meta.block.result$variates.A[, k]
    #for (q in which(n < ndefl)) if (sum(meta.block.result$loadings.A[[q]] != 0) <= 1)
    # warning(sprintf("Deflation failed because only one variable was selected for block #", q, "! \n"))
    
    
    # deflation if there are more than 1 component and if we haven't reach the max number of component(N)
    if (N != 1 & n != N) {
      defla.result <- defl.select(meta.block.result$variates.A, R, ndefl, n, nbloc = J, indY = indY, mode = mode, aa = meta.block.result$loadings.A)
      R <- defla.result$resdefl
      defl.matrix[[n + 1]] <- R
    }
    
    for (k in 1 : J) {
      if (N != 1) {
        P[[k]][, n - 1] <- defla.result$pdefl[[k]]
      }
      loadings.A[[k]][, n] <- meta.block.result$loadings.A[[k]]
    }
    
    if (n == 1) {
      for (k in 1 : J) loadings.Astar[[k]][, n] <- meta.block.result$loadings.A[[k]]
    } else {
      for (k in 1 : J) loadings.Astar[[k]][, n] <- meta.block.result$loadings.A[[k]] - loadings.Astar[[k]][, (1 : n - 1), drop = F] %*% drop(t(loadings.A[[k]][, n]) %*% P[[k]][, 1 : (n - 1), drop = F])
    }
    iter=c(iter,meta.block.result$iter)
  }
  
  if (verbose)
    cat(paste0("Computation of the SGCCA block components #", N , " is under progress...\n"))
  
  shave.matlist <- function(mat_list, nb_cols) mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY = FALSE)
  shave.veclist <- function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
  
  for (k in 1:J) {
    rownames(loadings.A[[k]]) = rownames(loadings.Astar[[k]])=colnames(A[[k]])
    rownames(variates.A[[k]]) = rownames(A[[k]]) #= rownames(variates.partial.A[[k]])
    colnames(variates.A[[k]]) = colnames(loadings.A[[k]]) = paste0("comp ", 1:max(ncomp))
    AVE_X[[k]] = apply(cor(A[[k]], variates.A[[k]])^2, 2, mean)
    if(is.null(tau))
    {
        names(loadings.partial.A[[k]])= names(variates.partial.A[[k]]) = levels(study)
        for(m in 1:nlevels(study))
        {
            rownames(loadings.partial.A[[k]][[m]])=colnames(A[[k]])
            colnames(loadings.partial.A[[k]][[m]])=paste0("comp ", 1:max(ncomp))
            rownames(variates.partial.A[[k]][[m]])=rownames(mean_centered[[1]]$data.list.study.scale[[m]])
            colnames(variates.partial.A[[k]][[m]])=paste0("comp ", 1:max(ncomp))
        }
    }
  }
  
  outer = matrix(unlist(AVE_X), nrow = max(ncomp))
  for (j in 1 : max(ncomp)) AVE_outer[j] <- sum(pjs * outer[j, ])/sum(pjs)
  variates.A = shave.matlist(variates.A, ncomp)
  AVE_X = shave.veclist(AVE_X, ncomp)
  AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)
  
  ### Start: output
  names(loadings.A)= names(variates.A) = names(A)
  
  if(is.null(tau))
  names(loadings.partial.A) = names(variates.partial.A) = names(A)
  
  names = lapply(1:J, function(x) {colnames(A[[x]])})
  names(names) = names(A)
  names[[length(names) + 1]] = row.names(A[[1]])
  names(names)[length(names)] = "indiv"
  
  #for (j in 1 : J) {
      #ifelse(is.null(names(A)[j]), names(c)[j] <- paste("X", j, sep = ""), names(c)[j] <- names(A)[j])# not happening as the check put names(A) by default
      #ifelse(is.null(names(A)[j]), names(b)[j] <- paste("X", j, sep = ""), names(b)[j] <- names(A)[j])# not happening as the check put names(A) by default
    
    #colnames(c[[j]]) = paste("comp", c(1:max(ncomp)));
    # colnames(loadings.A[[j]]) = paste("comp", c(1:max(ncomp))); #colnames(b[[j]]) = paste("comp", c(1:max(ncomp)))
    #colnames(variates.A[[j]]) = paste("comp", c(1:max(ncomp)));
    #row.names(c[[j]]) = colnames(A[[j]]);
    #row.names(loadings.A[[j]]) = colnames(A[[j]]); #row.names(b[[j]]) = colnames(A[[j]])
    #row.names(variates.A[[j]]) = rownames(A[[j]]);
    #}
  ### End: Output
  
  ### Start: Update names list with mixOmics package
  out <- list(X = if (is.null(indY)){A} else {A[-indY]}, Y = if (is.null(indY)){NULL} else {A[indY]}, ncomp = ncomp, mode = mode,
              keepA = keepA, keepA.constraint = keepA.constraint,
              variates = variates.A, loadings = shave.matlist(loadings.A, ncomp),
              variates.partial= if(is.null(tau)) {variates.partial.A} ,loadings.partial= if(is.null(tau)) {loadings.partial.A},
              loadings.star = shave.matlist(loadings.Astar, ncomp),
              names = names,tol = tol, iter=iter, nzv = if(near.zero.var) nzv.A,
              design = design,
              scheme = scheme,  crit = crit, AVE = AVE, defl.matrix = defl.matrix,
              init = init, bias = bias,
              scale = scale, tau = if(!is.null(tau)) tau.rgcca,study=study)
  ### End: Update names list with mixOmics package
  
  return(out)
}

# ----------------------------------------------------------------------------------------------------------
# sgccak - Runs sgccak() modified from RGCCA
#   inputs: A - list of datasets each with the same number of rows (samples)
#           design - design matrix
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ----------------------------------------------------------------------------------------------------------

sparse.meta.block_iteration <- function (A, design, study = NULL, keepA.constraint = NULL, keepA = NULL,  indY = NULL,
                                       scheme = "centroid", init = "svd", max.iter = 500, tol = 1e-06, verbose = TRUE, bias = FALSE)
{
  
  # keepA.constraint is a list of positions in A of the variables to keep on the component
  # keepA is a vector of numbers
  # study is a vector
  
  
  # == == == == no check needed as this function is only used in sgcca, in which the checks are conducted

  
  ### Check design matrix
  #stop("L1 constraints must be between 0 and 1.")


  ### Start: Initialization parameters
  J <- length(A); J2 <- J-1; pjs = sapply(A, NCOL)
  AVE_X <- rep(0, J);
  
  
  iter <- 1
  converg <- crit <- numeric()
  variates.A <- Z <- matrix(0, NROW(A[[1]]), J)
  misdata = any(sapply(A, function(x){any(is.na(x))})) # Detection of missing data
  g <- function(x) switch(scheme, horst = x, factorial = x^2, centroid = abs(x))
  

  # study split
  A_split=lapply(A, study_split, study)
  nlevels_study=nlevels(study)
  ### End: Initialization parameters

  ### Start: Initialisation "loadings.A" vector
  if (init == "svd") {
    ### Start: Change initialization of loadings.A
     if (misdata) {
       M = lapply(c(1:(J-1)), function(x){crossprod(replace(A[[x]], is.na(A[[x]]), 0), replace(A[[J]], is.na(A[[J]]), 0))})
     } else {
       M = lapply(c(1:(J-1)), function(x){crossprod(A[[x]], A[[J]])})
     }
     
     svd.M = lapply(M, function(x){svd(x, nu = 1, nv = 1)})
     loadings.A = lapply(c(1:(J-1)), function(x){svd.M[[x]]$u})
     loadings.A[[J]] = svd.M[[1]]$v
  } else if (init=="svd.single")
  {
      loadings.A <-  lapply(1 : J, function(y){initsvd(lapply(y, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])})
      
    ### End: Change initialization of a
    
    ### a <- lapply(lapply(c(1:J), function(x) {replace(A[[x]], is.na(A[[x]]), 0)}), function(x) return(svd(x, nu = 0, nv = 1)$v))
  } else {
    stop("init should be either 'svd' or 'svd.single'.")
  }
  ### End: Initialisation "a" vector
  

  variates.partial.A.comp = NULL
  loadings.partial.A.comp = list()
  
  for (q in 1:J)
  {
    variates.A[, q] <- A[[q]]%*%loadings.A[[q]]#apply(A[[q]], 1, crossprod, loadings.A[[q]])
    loadings.A[[q]] <- l2.norm(as.vector(loadings.A[[q]]))
    loadings.partial.A.comp[[q]] = list()
  }
  loadings.A_old <- loadings.A
  

  ### Start Algorithm 1 Sparse generalized canonical analysis (See Variable selection for generalized canonical correlation analysis (Tenenhaus))
  repeat {
    #variates.Aold can be used for the convergence of the algorithm, not used at the moment, so please disregard variates.Aold
    variates.Aold <- variates.A
    for (q in 1:J)
    {

      ### Start : !!! Impact of the diag of the design matrix !!! ###
      if (scheme == "horst") CbyCovq <- design[q, ]
      if (scheme == "factorial") CbyCovq <- design[q, ] * cov2(variates.A, variates.A[, q], bias = bias)
      if (scheme == "centroid") CbyCovq <- design[q, ] * sign(cov2(variates.A, variates.A[, q], bias = bias))
      ### End : !!! Impact of the diag of the design matrix !!! ###
      
      ### Step A start: Compute the inner components
      Z[, q] <- rowSums(mapply("*", CbyCovq, as.data.frame(variates.A)))
      Z_split=study_split(Z[,q,drop=FALSE],study)  # split Z by the study factor
      ### Step A end: Compute the inner components
      

      ### Step B start: Computer the outer weight ###
      # possibility of removing NA (replacing by 0) and use crossprod, further development
      temp=0
      for (m in 1:nlevels_study)
      {
        loadings.partial.A.comp[[q]][[m]] <-t(A_split[[q]][[m]])%*%Z_split[[m]] #apply(t(A_split[[q]][[m]]), 1, crossprod, Z_split[[m]])
        temp=temp+loadings.partial.A.comp[[q]][[m]]
      }
      loadings.A[[q]]=temp
      

      # check the following
      #loadings.partial.A.comp = lapply(1 : J, function(x){lapply(1 : nlevels_study, function(y){ A_split[[x]][[y]] %*% Z_split[[x]]})})
      
      
      loadings.A[[q]]=sparsity(loadings.A[[q]],keepA[[q]],keepA.constraint[[q]])
      loadings.A[[q]]=l2.norm(as.vector(loadings.A[[q]]))
      
      ### Step B end: Computer the outer weight ###
      variates.A[, q] <- A[[q]]%*%loadings.A[[q]] #apply(A[[q]], 1, crossprod, loadings.A[[q]])
    }
    
    crit[iter] <- sum(design * g(cov2(variates.A, bias = bias)))
    
    if (iter > max.iter)
      warning(cat("The SGCCA algorithm did not converge after", max.iter ,"iterations."))
    
    ### Start: Match algorithm with mixOmics algo (stopping point)
    ### if ((converg[iter] < tol & sum(stationnary_point) == J) | iter > max.iter)
    if (max(sapply(1:J, function(x){crossprod(loadings.A[[x]] - loadings.A_old[[x]])})) < tol | iter > max.iter)
    {break}
    ### End: Match algorithm with mixOmics algo (stopping point)
    
    
    loadings.A_old <- loadings.A
    iter <- iter + 1
  }
  ### End Algorithm 1 (See Variable selection for generalized canonical correlation analysis (Tenenhaus))
  
  
  #calculation variates.partial.A.comp
  variates.partial.A.comp = lapply(1 : J, function(x){lapply(1 : nlevels_study, function(y){ A_split[[x]][[y]] %*% loadings.A[[x]]})})  
  
  if (verbose)
    plot(crit, xlab = "iteration", ylab = "criteria")
  
  AVE_inner <- sum(design * cor(variates.A)^2/2)/(sum(design)/2)
  
  result <- list(variates.A = variates.A, loadings.A = loadings.A, crit = crit[which(crit != 0)],
                 AVE_inner = AVE_inner, loadings.partial.A.comp=loadings.partial.A.comp,variates.partial.A.comp=variates.partial.A.comp,iter=iter )
  return(result)
}

# ----------------------------------------------------------------------------------------------------------
# rgccak - Runs sgccak() modified from RGCCA
#   inputs: A - list of datasets each with the same number of rows (samples)
#           design - design matrix
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ----------------------------------------------------------------------------------------------------------

sparse.rgcca_iteration <- function (A, design, tau = "optimal", scheme = "centroid", scale = FALSE, max.iter = 500,
                    verbose = FALSE, init = "svd.single", bias = FALSE, tol = .Machine$double.eps, keepA = NULL,keepA.constraint = NULL) {
  ### Start: Initialisation parameters
  A <- lapply(A, as.matrix)
  J <- length(A)
  n <- NROW(A[[1]])
  pjs <- sapply(A, NCOL)
  variates.A <- matrix(0, n, J)
  ### End: Initialisation parameters
  
  if (!is.numeric(tau))
    tau = sapply(A, tau.estimate)
    
  loadings.A <- alpha <- M <- Minv <- K <- list()
  which.primal <- which((n >= pjs) == 1)
  which.dual <- which((n < pjs) == 1)
  
  if (init == "svd.single") {
    for (j in which.primal) {           
      loadings.A[[j]] <- initsvd(lapply(j, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])
    }
    for (j in which.dual) {
      alpha[[j]] <- initsvd(lapply(j, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])
      K[[j]] <- A[[j]] %*% t(A[[j]])
    }
  #} else if (init == "random") {
  # for (j in which.primal) {
  #    loadings.A[[j]] <- rnorm(pjs[j])
  #  }
  #  for (j in which.dual) {
  #    alpha[[j]] <- rnorm(n)
  #    K[[j]] <- A[[j]] %*% t(A[[j]])
  #  }
  } else {
      stop("init should be 'svd.single'.")#either random or by SVD.")
  }
  
  N = ifelse(bias, n, n - 1)
  for (j in which.primal) {      
    M[[j]] <- ginv(tau[j] * diag(pjs[j]) + (1 - tau[j]) * cov2(A[[j]], bias = bias))    
    loadings.A[[j]] <- drop(1/sqrt(t(loadings.A[[j]]) %*% M[[j]] %*% loadings.A[[j]])) * M[[j]] %*% loadings.A[[j]]
    variates.A[, j] <- A[[j]] %*% loadings.A[[j]]
  }
  
  for (j in which.dual) {
    M[[j]] = tau[j] * diag(n) + (1 - tau[j])/(N) * K[[j]]
    Minv[[j]] = ginv(M[[j]])
    alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% M[[j]] %*% K[[j]] %*% alpha[[j]])) * alpha[[j]]       
    loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
    variates.A[, j] = A[[j]] %*% loadings.A[[j]]
  }
  
  iter = 1
  crit = numeric()
  # converg = numeric()
  Z = matrix(0, NROW(A[[1]]), J)
  loadings.A_old = loadings.A
  g <- function(x) switch(scheme, horst = x, factorial = x^2, centroid = abs(x))
  
  repeat {
    variates.Aold <- variates.A
    
    for (j in which.primal) {
      if (scheme == "horst") CbyCovq <- design[j, ]
      if (scheme == "factorial") CbyCovq <- design[j, ] * cov2(variates.A, variates.A[, j], bias = bias)
      if (scheme == "centroid") CbyCovq <- design[j, ] * sign(cov2(variates.A, variates.A[, j], bias = bias))     
      
      # Compute the inner components
      Z[, j] = rowSums(mapply("*", CbyCovq, as.data.frame(variates.A)))
      
      # Computer the outer weight
      loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
      
      # sparse using keepA
      loadings.A[[j]]=sparsity(loadings.A[[j]], keepA[[j]], keepA.constraint = keepA.constraint[[j]])
      
      # Update variate
      variates.A[, j] = A[[j]] %*% loadings.A[[j]]
    }
    
    for (j in which.dual) {
      if (scheme == "horst") CbyCovq <- design[j, ]
      if (scheme == "factorial") CbyCovq <- design[j, ] * cov2(variates.A, variates.A[, j], bias = bias)
      if (scheme == "centroid") CbyCovq <- design[j, ] * sign(cov2(variates.A, variates.A[, j], bias = bias))   
      
      # Compute the inner components
      Z[, j] = rowSums(mapply("*", CbyCovq, as.data.frame(variates.A)))
      
      # Compute the outer weight
      alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
      loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
      
      # sparse using keepA
      loadings.A[[j]]=sparsity(loadings.A[[j]], keepA[[j]], keepA.constraint = keepA.constraint[[j]])
      
      # Update variate
      variates.A[, j] = A[[j]] %*% loadings.A[[j]]
    }
    
    crit[iter] <- sum(design * g(cov2(variates.A, bias = bias)))
    if (iter > max.iter)
      warning(cat("The RGCCA algorithm did not converge after", max.iter ,"iterations."))
    
    ### Start: Match algorithm with mixOmics algo (stopping point)
    ### if ((converg[iter] < tol & sum(stationnary_point) == J) | iter > max.iter)
    if (max(sapply(1:J, function(x){crossprod(loadings.A[[x]] - loadings.A_old[[x]])})) < tol | iter > max.iter)
      break
    ### End: Match algorithm with mixOmics algo (stopping point)
    
    loadings.A_old <- loadings.A
    iter <- iter + 1
  }
    
  if (verbose) 
    plot(crit, xlab = "iteration", ylab = "criteria")
  
  AVE_inner <- sum(design * cor(variates.A)^2/2)/(sum(design)/2)
  
  result <- list(variates.A = variates.A, loadings.A = loadings.A, crit = crit[which(crit != 0)], 
                 AVE_inner = AVE_inner, design = design, tau = tau, scheme = scheme,iter=iter,
                 keepA=keepA,keepA.constraint=keepA.constraint )
  return(result)
}

