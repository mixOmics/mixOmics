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
# Parallel functions
#   parallel.valid.sgccda(),
#
# Functions modified from RGCCA R-library
#   sgcca(), sgccak()
#
# Functions acquired from RGCCA R-library
#   cov2(), scale2(), initsvd(), crossprod(), soft.threshold(), BinarySearch(),
#   soft(), norm2(), defl.select(), deflation(),
#
# Functions acquired from gdata library
#   cbindX()
#
#############################################################################################################

#############################################################################################################
# Functions modified from RGCCA
#############################################################################################################
# ----------------------------------------------------------------------------------------------------------
# sgcca - Runs sgcca() modified from RGCCA
#   inputs: A - list of datasets each with the same number of rows (samples)
#           design - design matrix
#           lambda - vector of shirnkage penalties
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ----------------------------------------------------------------------------------------------------------
meta.block.spls = function (A, indY = NULL,  design = 1 - diag(length(A)), lambda = NULL,tau=NULL,#rep(1, length(A)),
ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE,  bias = FALSE,
init = "svd", tol = 1e-06, verbose = FALSE,
mode = "canonical", sparse = FALSE, max.iter = 500,study = NULL, keepA = NULL,
keepA.constraint = NULL, near.zero.var = FALSE) { # meta.hybrid.spls
    
    # A: list of matrices
    # indY: integer, pointer to one of the matrices of A
    # design: design matrix, links between matrices. Diagonal must be 0
    # lambda: Matrix of shrinkage penalties. Each row is for a component, each column for a matrix of A. Can be a vector, in that case same penalty to each component. Only used if keepA is missing
    # ncomp: vector of ncomp, per matrix
    # scheme: a function "g", refer to the article (thanks Benoit)
    # scale: do you want to scale ? mean is done by default and cannot be changed (so far)
    # bias: scale the data with n or n-1
    # init: one of "svd" or "random", initialisation of the algorithm
    # tol: nobody cares about this
    # verbose: show the progress of the algorithm
    # mode: canonical, classic, invariant, regression
    # sparse: use sgcca instead of rgcca when you use sparse function
    # max.iter: nobody cares about this
    # study: factor for each matrix of A, must be a vector
    # keepA: keepX of spls for each matrix of A. must be a list. Each entry must be of the same length (max ncomp)
    # keepA.constraint: keepX.constraint, which variables are kept on the first num.comp-1 components. It is a list of characters
    # near.zero.var: do you want to remove variables with very small variance
    
    
    # add check so that this function can run on its own (which shouldn't happen now, but don't know the future)
    
    
    check=Check.entry.meta.block.spls(A, indY, design , lambda,ncomp , scheme , scale ,  bias,
            init , tol , verbose,mode, sparse , max.iter,study , keepA, keepA.constraint)
    A=check$A
    ncomp=check$ncomp
    study=check$study
    

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
                    warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv$X for problematic predictors.")
                    if(ncol(A[[q]]) == 0) {stop("No more variables in X")}
                
                    # at this stage, keepA.constraint need to be numbers
                    if(!is.null(keepA.constraint[[q]]))
                    {
                        #remove the variables from keepA.constraint if removed by near.zero.var
                        keepA.constraint[[q]]=match.keepX.constraint(names.remove.X,keepA.constraint[[q]])
                        # replace character by numbers
                        keepA.constraint[[q]]= lapply(keepA.constraint[[q]],function(x){match(x,colnames(A[[q]]))})
                    }
                }
            }
    }
    
    # keepA is updated to be of length ncomp, the first entries correspond to the keepA.constraint if it was provided
    
    for(q in 1:length(A))
    keepA[[q]]=c(unlist(lapply(keepA.constraint[[q]],length)),keepA[[q]]) #of length ncomp, can contains 0

    if(sparse==FALSE) tau=lambda
    
   
   
    
    # center the data per study, per matrix of A, scale if scale=TRUE, option bias
    mean_centered = lapply(A, function(x){mean_centering_per_study(x, study, scale, bias)})
    A = lapply(mean_centered, function(x){x$concat.data})
    
    ### Start: Initialization parameters
    pjs <- sapply(A, NCOL); nb_ind <- NROW(A[[1]]);
    J <- length(A); R <- A; N <- max(ncomp) # R: residuals matrices, will be a list of length ncomp
    
    AVE_inner <- AVE_outer <- rep(NA, max(ncomp))
    defl.matrix <- AVE_X <- crit <- loadings.partial.A <- variates.partial.A <-list()
    P <- loadings.A <- loadings.Astar <- c <- t <- b <- variates.A <- NULL
    for (k in 1:J) t[[k]] <- variates.A[[k]] <- matrix(NA, nb_ind, N)
    for (k in 1:J) P[[k]] <- loadings.A[[k]] <- loadings.Astar[[k]] <- b[[k]] <- c[[k]] <- matrix(NA, pjs[[k]], N)
    
    defl.matrix[[1]] <- A; ndefl <- ncomp - 1;  J2 <- J-1;
    
    if(is.vector(lambda))
    {
        lambda=t(matrix(rep(lambda,N),ncol=N))
    }
    ### End: Initialization parameters
    
    iter=NULL
    for (n in 1:N)
    {
        if (verbose)
        cat(paste0("Computation of the SGCCA block components #", n, " is under progress... \n"))
        meta.block.result = NULL
        
        ### Start: Estimation ai
        if (sparse == TRUE)
        {
            meta.block.result <- meta.block.spls_iteration(R, design, lambda[n, ],study = study,  keepA.constraint=if (!is.null(keepA.constraint)) {lapply(keepA.constraint, function(x){unlist(x[n])})} else {NULL} ,
            keepA = if (!is.null(keepA)) {lapply(keepA, function(x){x[n]})} else {NULL},indY = indY,
            scheme = scheme, init = init, max.iter = max.iter, tol = tol,   verbose = verbose)
        } else {
            meta.block.result <- rgccak(R, design, lambda[n, ], scheme = scheme, init = init, tol = tol, verbose = verbose, max.iter = max.iter)
        }
        
        ### End: Estimation ai
        
        loadings.partial.A[[n]]=meta.block.result$loadings.partial.A.comp
        variates.partial.A[[n]]=meta.block.result$variates.partial.A.comp
        
        AVE_inner[n] <- meta.block.result$AVE_inner
        crit[[n]] <- meta.block.result$crit
        for (k in 1:J) variates.A[[k]][, n] <- meta.block.result$variates.A[, k]
        for (q in which(n < ndefl)) if (sum(meta.block.result$loadings.A[[q]] != 0) <= 1)
        warning(sprintf("Deflation failed because only one variable was selected for block #", q, "! \n"))
        
        
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
        rownames(loadings.A[[k]]) = rownames(loadings.Astar[[k]]) = colnames(A[[k]])
        rownames(variates.A[[k]]) = rownames(A[[k]])
        colnames(variates.A[[k]]) = paste0("comp ", 1:max(ncomp))
        AVE_X[[k]] = apply(cor(A[[k]], variates.A[[k]])^2, 2, mean)
    }
    
    outer = matrix(unlist(AVE_X), nrow = max(ncomp))
    for (j in 1 : max(ncomp)) AVE_outer[j] <- sum(pjs * outer[j, ])/sum(pjs)
    variates.A = shave.matlist(variates.A, ncomp)
    AVE_X = shave.veclist(AVE_X, ncomp)
    AVE <- list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)
    
    ### Start: output
    names(loadings.A) = names(A); names(variates.A) = names(A)
    
    names = lapply(1:J, function(x) {colnames(A[[x]])})
    names(names) = names(A)
    names[[length(names) + 1]] = row.names(A[[1]])
    names(names)[length(names)] = "indiv"
    
    for (j in 1 : J) {
        ifelse(is.null(names(A)[j]), names(c)[j] <- paste("X", j, sep = ""), names(c)[j] <- names(A)[j])
        ifelse(is.null(names(A)[j]), names(b)[j] <- paste("X", j, sep = ""), names(b)[j] <- names(A)[j])
        
        colnames(c[[j]]) = paste("comp", c(1:max(ncomp))); colnames(loadings.A[[j]]) = paste("comp", c(1:max(ncomp))); colnames(b[[j]]) = paste("comp", c(1:max(ncomp)))
        row.names(c[[j]]) = colnames(A[[j]]); row.names(loadings.A[[j]]) = colnames(A[[j]]); row.names(b[[j]]) = colnames(A[[j]])
    }
    ### End: Output
    
    ### Start: Update names list with mixOmics package
    out <- list(X = if (is.null(indY)){A} else {A[-indY]}, Y = if (is.null(indY)){NULL} else {A[indY]}, mode = mode, names = names,
    variates = variates.A, loadings = shave.matlist(loadings.A, ncomp),
    variates.partial=variates.partial.A,loadings.partial=loadings.partial.A,
    loadings.star = shave.matlist(loadings.Astar, ncomp), design = design, tau = if (sparse == FALSE) {lambda} else {NULL},
    penalty = if (sparse == FALSE) {NULL} else {lambda}, max.iter = max.iter, scheme = scheme, ncomp = ncomp, crit = crit, AVE = AVE,
    defl.matrix=defl.matrix, tol = tol, keepA = keepA, keepA.constraint = keepA.constraint, init = init, bias = bias, scale = scale,near.zero.var= if(near.zero.var) nzv.A,iter=iter)
    ### End: Update names list with mixOmics package
    
    return(out)
}

# ----------------------------------------------------------------------------------------------------------
# sgccak - Runs sgccak() modified from RGCCA
#   inputs: A - list of datasets each with the same number of rows (samples)
#           design - design matrix
#           lambda - vector of shirnkage penalties
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ----------------------------------------------------------------------------------------------------------

meta.block.spls_iteration <- function (A, design, lambda = NULL,#rep(1, length(A)),
study = NULL, keepA.constraint = NULL, keepA = NULL,  indY = NULL,
scheme = "centroid", init = "svd", max.iter = 500, tol = 1e-06, verbose = TRUE, bias = FALSE)
{

    # keepA.constraint is a list of positions in A of the variables to keep on the component
    # keepA is a vector of numbers
    # study is a vector
    
    
    # == == == == no check needed as this function is only used in sgcca, in which the checks are conducted
    # needs checks
    #if(!is.null(lambda) | !is.null(keepA) | !is.null(keepA))    {sparse=TRUE}
    #check are different whether sparse==TRUE
    
    
    ### Check design matrix
    #if (!is.null(lambda) & (any(lambda < 0 | lambda > 1)))
    #stop("L1 constraints must be between 0 and 1.")
    
    
    ### Start: Initialization parameters
    J <- length(A); J2 <- J-1; pjs = sapply(A, NCOL)
    AVE_X <- rep(0, J);
    if(!is.null(lambda))
    {
        lambda <- lambda * sqrt(pjs)
    }
    
    
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
        ### End: Change initialization of a
        
        ### a <- lapply(lapply(c(1:J), function(x) {replace(A[[x]], is.na(A[[x]]), 0)}), function(x) return(svd(x, nu = 0, nv = 1)$v))
    } else if (init == "random") {
        loadings.A <- lapply(pjs, rnorm)
    } else {
        stop("init should be either random or svd.")
    }
    ### End: Initialisation "a" vector
    

    variates.partial.A.comp = NULL
    loadings.partial.A.comp = list()
    
    for (q in 1:J)
    {
        variates.A[, q] <- A[[q]]%*%loadings.A[[q]]#apply(A[[q]], 1, crossprod, loadings.A[[q]])
        loadings.A[[q]] <- l2.norm(as.vector(loadings.A[[q]]))
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
                loadings.partial.A.comp[[m]] <-t(A_split[[q]][[m]])%*%Z_split[[m]] #apply(t(A_split[[q]][[m]]), 1, crossprod, Z_split[[m]])
                temp=temp+loadings.partial.A.comp[[m]]
            }
            loadings.A[[q]]=temp
            
            # check the following
            #loadings.partial.A.comp = lapply(1 : J, function(x){lapply(1 : nlevels_study, function(y){ A_split[[x]][[y]] %*% Z_split[[x]]})})


            loadings.A[[q]]=sparsity(loadings.A[[q]],keepA[[q]],keepA.constraint[[q]],lambda)
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
#           lambda - vector of shirnkage penalties
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ----------------------------------------------------------------------------------------------------------

rgccak <- function (A, design, tau = "optimal", scheme = "centroid",  max.iter = 500,
verbose = FALSE, init = "svd",  tol = .Machine$double.eps) {
    
    
    # add checks
    
    A <- lapply(A, as.matrix)
    J <- length(A)
    n <- NROW(A[[1]])
    pjs <- sapply(A, NCOL)
    variates.A <- matrix(0, n, J)
    
    
    
    if (!is.numeric(tau))
    tau = sapply(A, tau.estimate)
    
    loadings.A <- alpha <- M <- Minv <- K <- list()
    which.primal <- which((n >= pjs) == 1)
    which.dual <- which((n < pjs) == 1)
    
    if (init == "svd") {
        for (j in which.primal) {
            loadings.A[[j]] <- initsvd(lapply(j, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]]) ### Handle missing data
        }
        for (j in which.dual) {
            alpha[[j]] <- initsvd(lapply(j, function(x) {replace(A[[x]], is.na(A[[x]]), 0)})[[1]]) ### Handle missing data
            K[[j]] <- do.call(mapply, c(cbind, lapply(1:n, function(x){apply(A[[j]], 1, crossprod, A[[j]][x, ])}))) ### Handle missing data
            ### K[[j]] <- A[[j]] %*% t(A[[j]])
        }
    } else if (init == "random") {
        for (j in which.primal) {
            loadings.A[[j]] <- rnorm(pjs[j])
        }
        for (j in which.dual) {
            alpha[[j]] <- rnorm(n)
            K[[j]] <- do.call(mapply, c(cbind, lapply(1:n, function(x){apply(A[[j]], 1, crossprod, A[[j]][x, ])}))) ### Handle missing data
            ### K[[j]] <- A[[j]] %*% t(A[[j]])
        }
    } else {
        stop("init should be either random or by SVD.")
    }
    
    N = ifelse(bias, n, n - 1)
    for (j in which.primal) {
        ifelse(tau[j] == 1, {
            loadings.A[[j]] <- drop(1/sqrt(t(loadings.A[[j]]) %*% loadings.A[[j]])) * loadings.A[[j]]
            variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
            ### variates.A[, j] <- A[[j]] %*% loadings.A[[j]]
        }, {
            M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])/(N)) * (do.call(mapply, c(cbind, lapply(1:NROW(t(A[[j]])), function(x){apply(t(A[[j]]), 1, crossprod, A[[j]][, x])}))))) ### Handle missing data
            ### M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])/(N)) * (t(A[[j]]) %*% A[[j]]))
            loadings.A[[j]] <- drop(1/sqrt(t(loadings.A[[j]]) %*% M[[j]] %*% loadings.A[[j]])) * M[[j]] %*% loadings.A[[j]]
            variates.A[, j] <- apply(A[[j]], 1, crossprod, a[[j]]) ### Handle missing data
            ### variates.A[, j] <- A[[j]] %*% loadings.A[[j]]
        })
    }
    
    for (j in which.dual) {
        ifelse(tau[j] == 1, {
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% K[[j]] %*% alpha[[j]])) * alpha[[j]]
            loadings.A[[j]] = apply(t(A[[j]]), 1, crossprod, alpha[[j]]) ### Handle missing data
            ### a[[j]] = t(A[[j]]) %*% alpha[[j]]
            variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
            ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
        }, {
            M[[j]] = tau[j] * diag(n) + (1 - tau[j])/(N) * K[[j]]
            Minv[[j]] = ginv(M[[j]])
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% M[[j]] %*% K[[j]] %*% alpha[[j]])) * alpha[[j]]
            loadings.A[[j]] = apply(t(A[[j]]), 1, crossprod, alpha[[j]]) ### Handle missing data
            ### a[[j]] = t(A[[j]]) %*% alpha[[j]]
            variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
            ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
        })
    }
    
    iter = 1
    crit = numeric()
    converg = numeric()
    Z = matrix(0, NROW(A[[1]]), J)
    loadings.A_temp = loadings.A
    g <- function(x) switch(scheme, horst = x, factorial = x^2, centroid = abs(x))
    
    repeat {
        variates.Aold <- variates.A
        if (scheme == "horst") {
            for (j in which.primal) {
                ifelse(tau[j] == 1, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * variates.A)
                    loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% do.call(mapply, c(cbind, lapply(1:n, function(x){apply(A[[j]], 1, crossprod, A[[j]][x, ])}))) %*% Z[, j])) * t(A[[j]]) %*% Z[, j] ### Handle missing data
                    ### loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * t(A[[j]]) %*% Z[, j]
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                }, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * variates.A)
                    loadings.A[[j]] = drop(1/sqrt(t(apply(t(A[[j]]), 1, crossprod, Z[, j])) %*% M[[j]] %*% apply(t(A[[j]]), 1, crossprod, Z[, j]))) * M[[j]] %*% apply(t(A[[j]]), 1, crossprod, Z[, j]) ### Handle missing data
                    ### loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * M[[j]] %*% t(A[[j]]) %*% Z[, j]
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                })
            }
            for (j in which.dual) {
                ifelse(tau[j] == 1, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * variates.A)
                    alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
                    loadings.A[[j]] = apply(t(A[[j]]), 1, crossprod, alpha[[j]]) ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                }, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * variates.A)
                    alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                    loadings.A[[j]] = apply(t(A[[j]]), 1, crossprod, alpha[[j]]) ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                })
            }
        }
        
        if (scheme == "factorial") {
            for (j in which.primal) {
                ifelse(tau[j] == 1, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov(variates.A[, j],  variates.A), n), n, J, byrow = TRUE) * variates.A)
                    loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% do.call(mapply, c(cbind, lapply(1:n, function(x){apply(A[[j]], 1, crossprod, A[[j]][x, ])}))) %*% Z[, j])) * t(A[[j]]) %*% Z[, j] ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * (t(A[[j]]) %*% Z[, j])
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                }, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov(variates.A[, j], variates.A), n), n, J, byrow = TRUE) * variates.A)
                    loadings.A[[j]] = drop(1/sqrt(t(apply(t(A[[j]]), 1, crossprod, Z[, j])) %*% M[[j]] %*% apply(t(A[[j]]), 1, crossprod, Z[, j]))) * M[[j]] %*% apply(t(A[[j]]), 1, crossprod, Z[, j]) ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                })
            }
            for (j in which.dual) {
                ifelse(tau[j] == 1, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov(variates.A[, j], variates.A), n), n, J, byrow = TRUE) * variates.A)
                    alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
                    loadings.A[[j]] = apply(t(A[[j]]), 1, crossprod, alpha[[j]]) ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                }, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * matrix(rep(cov(variates.A[, j], variates.A), n), n, J, byrow = TRUE) * variates.A)
                    alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                    loadings.A[[j]] = apply(t(A[[j]]), 1, crossprod, alpha[[j]]) ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                })
            }
        }
        if (scheme == "centroid") {
            for (j in which.primal) {
                ifelse(tau[j] == 1, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov(variates.A[, j], variates.A), n), n, J, byrow = TRUE)) * variates.A)
                    loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% do.call(mapply, c(cbind, lapply(1:n, function(x){apply(A[[j]], 1, crossprod, A[[j]][x, ])}))) %*% Z[, j])) * t(A[[j]]) %*% Z[, j] ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% t(A[[j]]) %*% Z[, j])) * (t(A[[j]]) %*% Z[, j])
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                }, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov(variates.A[, j], variates.A), n), n, J, byrow = TRUE)) * variates.A)
                    loadings.A[[j]] = drop(1/sqrt(t(apply(t(A[[j]]), 1, crossprod, Z[, j])) %*% M[[j]] %*% apply(t(A[[j]]), 1, crossprod, Z[, j]))) * M[[j]] %*% apply(t(A[[j]]), 1, crossprod, Z[, j]) ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*% t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                })
            }
            for (j in which.dual) {
                ifelse(tau[j] == 1, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov(variates.A[, j], variates.A), n), n, J, byrow = TRUE)) * variates.A)
                    alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
                    loadings.A[[j]] = apply(t(A[[j]]), 1, crossprod, alpha[[j]]) ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                }, {
                    Z[, j] = rowSums(matrix(rep(design[j, ], n), n, J, byrow = TRUE) * sign(matrix(rep(cov(variates.A[, j], variates.A), n), n, J, byrow = TRUE)) * variates.A)
                    alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                    loadings.A[[j]] = apply(t(A[[j]]), 1, crossprod, alpha[[j]]) ### Handle missing data
                    variates.A[, j] <- apply(A[[j]], 1, crossprod, loadings.A[[j]]) ### Handle missing data
                    ### loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
                    ### variates.A[, j] = A[[j]] %*% loadings.A[[j]]
                })
            }
        }
        num_converg <- sum((rowSums(variates.Aold) - rowSums(variates.A))^2)
        den_converg <- sum(rowSums(variates.Aold)^2)
        ### Start: Introduction of if/else to avoid division by 0 (den_converg)
        if (den_converg != 0){
            converg[iter] <- num_converg/den_converg
        } else {converg[iter] = 0}
        ### End: Introduction of if/else to avoid division by 0 (den_converg)
        stationnary_point = rep(FALSE, length(A))
        for (j in 1:J) stationnary_point[j] = sum(round(abs(loadings.A_temp[[j]] - loadings.A[[j]]), 8) < tol) == NCOL(A[[j]])
        loadings.A_temp <- loadings.A
        crit[iter] <- sum(design * g(cov2(variates.A, bias = bias)))
        if (iter > max.iter)
        warning(cat("The RGCCA algorithm did not converge after", max.iter, "iterations."))
        if ((converg[iter] < tol & sum(stationnary_point) == J | iter > max.iter))
        break
        iter <- iter + 1
    }
    if (sum(stationnary_point) == J & verbose)
    cat("The RGCCA algorithm converged to a fixed point of the stationary equations after", iter - 1, "iterations \n")
    if (verbose)
    plot(crit, xlab = "iteration", ylab = "criteria")
    AVEinner <- sum(design * cor(variates.A)^2/2)/(sum(design)/2)
    
    result <- list(variates.A = variates.A, loadings.A = loadings.A, crit = crit[which(crit != 0)], converg = converg[which(converg != 0)],
    AVE_inner = AVEinner, design = design, tau = tau, scheme = scheme)
    return(result)
}
