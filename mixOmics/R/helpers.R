
# --------------------------------------
# match.keepX.constraint after removing variables with near.zer.var for instance
# --------------------------------------
match.keepX.constraint=function(names.remove,keepX.constraint)
{
    #matching the new X (after removing some variables) and keepX.constraint
    if(length(names.remove)>0)
    {
        
        ind.match= lapply(keepX.constraint,function(x){match(x,names.remove)})
        
        if(sum(!is.na(unlist(ind.match)))>0)
        {
            warnings("at least one variable was removed from keepX.constraint because of a null variance. Please check object$keepX.constraint to see which variables are used.")
            #remove from keepX.constraint
            keepX.constraint=lapply(keepX.constraint,function(x){temp=match(x,names.remove); if(sum(!is.na(temp))>0){x=x[-which(!is.na(temp))]}else{x=x}})
        }
        
        keepX=unlist(lapply(keepX.constraint,length))
        if(any(keepX==0))
        {
            ind.min=min(which(keepX==0))
            ncomp=ind.min-1
            stop(paste("keepX.constraint was reduced from", length(keepX.constraint)," components to",ncomp," components. Please change keepX.constraint or put near.zero.var=FALSE and restart the function"))
            #construction of the new keepX.constraint, using ncomp components
            keepX.constraint.temp=keepX.constraint
            for(i in (ncomp+1):length(keepX.constraint)) keepX.constraint.temp[[i]]=NULL
            keepX.constraint=keepX.constraint.temp
            
        }
        
    }
    
    out=keepX.constraint
}


# --------------------------------------
# nearZeroVar
# --------------------------------------

nearZeroVar <- function (x, freqCut = 95/5, uniqueCut = 10) {
    
    if (is.vector(x))
    x = matrix(x, ncol = 1)
    
    freqRatio = apply(x, 2, function(data) {
        data = na.omit(data)
        
        if (length(unique(data)) == length(data)){ # No duplicate
            return(1)
        } else if (length(unique(data)) == 1) { # Same value
            return(0)
        } else {
            t = table(data)
            return(max(t, na.rm = TRUE)/max(t[-which.max(t)], na.rm = TRUE))
        }
    })
    
    lunique = apply(x, 2, function(data) length(unique(data[!is.na(data)])))
    
    percentUnique = 100 * lunique/rep(nrow(x))
    zeroVar = (lunique == 1) | apply(x, 2, function(data) all(is.na(data)))
    
    out = list()
    out$Position = which((freqRatio > freqCut & percentUnique <= uniqueCut) | zeroVar)
    names(out$Position) = NULL
    out$Metrics = data.frame(freqRatio = freqRatio, percentUnique = percentUnique)
    out$Metrics = out$Metrics[out$Position, ]
    return(out)
}



# --------------------------------------
# study_split: is used for meata.pls and meta.spls
# --------------------------------------
study_split = function(data, study)
{
    
    data = as.matrix(data)
    
    M=length(levels(study))
    P=ncol(data)
    
    #---------------------- split data
    data.list.study = split(data,study)
    if(!is.null(rownames(data))) {study.name = split(rownames(data),study)}
    
    for(m in 1:M)
    {
        data.list.study[[m]]=matrix(data.list.study[[m]], ncol=P)
        if(!is.null(colnames(data))) {colnames(data.list.study[[m]]) = colnames(data)}
        if(!is.null(rownames(data))) {rownames(data.list.study[[m]]) = study.name[[m]]}
    }
    
    
    result = data.list.study
    
    return(invisible(result))
    
}


# --------------------------------------
# soft_thresholding
# --------------------------------------
soft_thresholding_L1=function(x,nx)
{
    #if (nx != 0) {
    #    x = ifelse(abs(x) > abs(x[order(abs(x))][nx]),
    #    (abs(x) - abs(x[order(abs(x))][nx])) * sign(x), 0)
    #}
    
    #selection on a (loadings.X). modified on 19/02/15 to make sure that a!=0
    if(nx!=0)
    {
        absa = abs(x)
        if(any(rank(absa, ties.method = "max") <= nx)) {
            x = ifelse(rank(absa, ties.method = "max") <= nx, 0,
            sign(x) * (absa - max(absa[rank(absa, ties.method = "max") <= nx])))
        }
    }
    
    x
}

# ----------------------------------------------------------------------------------------------------------
# soft.threshold() - soft-thresholds a vector such that the L1-norm constraint is satisfied.
# ----------------------------------------------------------------------------------------------------------
soft.threshold = function (x, sumabs = 1)
return(soft(x, BinarySearch(x, sumabs)))

BinarySearch = function(argu,sumabs){
    if(norm2(argu)==0 || sum(abs(argu/norm2(argu)))<=sumabs) return(0)
    lam1 <- 0
    lam2 <- max(abs(argu))-1e-5
    iter <- 1
    while(iter < 500){
        su <- soft(argu,(lam1+lam2)/2)
        if(sum(abs(su/norm2(su)))<sumabs){
            lam2 <- (lam1+lam2)/2
        } else {
            lam1 <- (lam1+lam2)/2
        }
        if((lam2-lam1)<1e-10) return((lam1+lam2)/2)
        iter <- iter+1
    }
    warning("Didn't quite converge")
    return((lam1+lam2)/2)
}

soft = function(x,d) return(sign(x)*pmax(0, abs(x)-d))

norm2 = function(vec){
    a <- sqrt(sum(vec^2))
    if(a==0) a <- .05
    return(a)
}


# --------------------------------------
# sparsity function
# --------------------------------------
sparsity=function(loadings.A,keepA,keepA.constraint)
{
    
    if(!is.null(keepA.constraint))
    {
        loadings.A[-keepA.constraint]=0
    }else if(!is.null(keepA))
    {
        nx=length(loadings.A) - keepA
        loadings.A=soft_thresholding_L1(loadings.A,nx)
    }
    return(loadings.A)
}


# --------------------------------------
# Mean centering/scaling per study
# --------------------------------------
mean_centering_per_study=function(data,study,scale,bias=FALSE) {
    
    M = length(levels(study))   # number of groups
    # split the data
    data.list.study = study_split(data, study)
    
    # center and scale data per group, and concatene the data
    concat.data = NULL
    data.list.study.scale= list()
    for (m in 1:M) {
        data.list.study.scale[[m]] = scale(data.list.study[[m]], center = TRUE, scale = scale)
        if(bias)
        data.list.study.scale[[m]]=data.list.study.scale[[m]]*(sqrt((nrow(data.list.study.scale[[m]])-1)/nrow(data.list.study.scale[[m]])))
        
        if(sum(is.na(data.list.study.scale[[m]]))>0)
        data.list.study.scale[[m]][is.na(data.list.study.scale[[m]])]=0
        
        concat.data = rbind(concat.data, unlist(data.list.study.scale[[m]]))
    }
    #rename rows and cols of concatenated centered (and/or scaled) data
    colnames(concat.data) = colnames(data)
    
    #sort the samples as in the original X
    indice.match=match(rownames(data),rownames(concat.data))
    concat.data=concat.data[indice.match,]
    
    if (M > 1 )
    {
        for(m in 1:M)
        {
            attr(concat.data,paste0("means:", levels(study)[m]))=attr(data.list.study.scale[[m]],"scaled:center")
            attr(concat.data,paste0("sigma:", levels(study)[m]))=attr(data.list.study.scale[[m]],"scaled:scale")
        }
    } else {
        attr(concat.data,"scaled:center")=attr(data.list.study.scale[[m]],"scaled:center")
        attr(concat.data,"scaled:scale")=attr(data.list.study.scale[[m]],"scaled:scale")
    }
    
    return(list(concat.data=concat.data,data.list.study.scale=data.list.study.scale))
}


# --------------------------------------
# l2.norm
# --------------------------------------
l2.norm=function(x)
{
    if(!is.vector(x)) stop("x has to be a vector")
    out=x/drop(sqrt(crossprod(x)))
}



#############################################################################################################
# Functions acquired from RGCCA R-library
#############################################################################################################
# ----------------------------------------------------------------------------------------------------------
# cov2() - Compute biased and unbiased covariance and variance estimates
# ----------------------------------------------------------------------------------------------------------
cov2 = function (x, y = NULL, bias = TRUE) {
    n = NROW(x)
    if (is.null(y)) {
        x = as.matrix(x)
        if (bias) {
            C = ((n - 1)/n) * cov(x, use = "pairwise.complete.obs")
        } else {
            C = cov(x, use = "pairwise.complete.obs")
        }
    } else {
        if (bias) {
            C = ((n - 1)/n) * cov(x, y, use = "pairwise.complete.obs")
        } else {
            C = cov(x, y, use = "pairwise.complete.obs")
        }
    }
    return(C)
}

# ----------------------------------------------------------------------------------------------------------
# scale2() - scales columns (variables) of each datasets in list A
# ----------------------------------------------------------------------------------------------------------
scale2 = function (A, center = TRUE, scale = TRUE, bias = TRUE) {
    ### Start: add row.names to the output
    rownames = NULL
    if (!is.null(row.names(A)))
    rownames = row.names(A)
    ### End: add row.names to the output
    
    if (center == TRUE & scale == TRUE) {
        A = scale(A, center = TRUE, scale = FALSE)
        std = sqrt(apply(A, 2, function(x) cov2(x, bias = bias)))
        #if (any(std == 0)) {
        #    sprintf("there were %d constant variables", sum(std ==
        #    0))
        #    std[std == 0] = 1
        #}
        A = A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
        attr(A, "scaled:scale") = std
        ### return(A) Wait end script
    }
    if (center == TRUE & scale == FALSE) {
        A = scale(A, center = TRUE, scale = FALSE)
        ### return(A) Wait end script
    }
    if (center == FALSE & scale == TRUE) {
        std = apply(A, 2, function(x) cov2(x, bias = bias))
        A = A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
        attr(A, "scaled:scale") = std
        ### return(A) Wait end script
    }
    ### Start: add row.names to the output
    if (!is.null(rownames))
    row.names(A) = rownames
    return(A)
    ### End: add row.names to the output
}

# ----------------------------------------------------------------------------------------------------------
# initsvd() - performs SVD on matrix X
# ----------------------------------------------------------------------------------------------------------
initsvd = function (X) {
    n = NROW(X)
    p = NCOL(X)
    ifelse(n >= p, return(svd(X, nu = 0, nv = 1)$v), return(svd(X, nu = 1, nv = 0)$u))
}

# ----------------------------------------------------------------------------------------------------------
# miscrossprod() - Compute cross-product between vectors x and y
# ----------------------------------------------------------------------------------------------------------
miscrossprod = function (x, y) {
    d.p = sum(drop(x) * drop(y), na.rm = TRUE)
    #d.p = as.vector(d.p)/norm2(d.p)     ## change made
    return(d.p)
}


# ----------------------------------------------------------------------------------------------------------
# defl.select() - computes residual matrices
# ----------------------------------------------------------------------------------------------------------
defl.select = function(yy, rr, nncomp, nn, nbloc, indY = NULL, mode = "canonical", aa = NULL) { ### Start: Add new parameter for estimation classic mode
    resdefl <- NULL
    pdefl <- NULL
    for (q in 1 : nbloc) {
        ### Start: insertion of new deflations (See La regression PLS Theorie et pratique (page 139))
        if ( nn <= nncomp[q] ) {
            if ((mode == "canonical") || (q != indY)) {
                defltmp <- deflation(rr[[q]], yy[ , q])
                resdefl[[q]] <- defltmp$R
                pdefl[[q]]   <- defltmp$p
            } else if (mode == "classic") {
                resdefl[[q]] <- Reduce("+", lapply(c(1:nbloc)[-q], function(x) {rr[[q]] - yy[ ,x]%*%t(aa[[q]])}))/(nbloc-1)
                pdefl[[q]]   <-  rep(0,NCOL(rr[[q]]))
            } else if (mode == "invariant") {
                resdefl[[q]] <- rr[[q]]
                pdefl[[q]]   <-  rep(0,NCOL(rr[[q]]))
            } else if (mode == "regression") {
                resdefl[[q]] <- Reduce("+", lapply(c(1:nbloc)[-q], function(x) {deflation(rr[[q]],yy[, x])$R}))/(nbloc-1)
                pdefl[[q]]   <-  rep(0,NCOL(rr[[q]]))
            }
            ### End: insertion of new deflations (See La regression PLS Theorie et pratique (page 139))
        } else {
            resdefl[[q]] <- rr[[q]]
            pdefl[[q]]   <-  rep(0,NCOL(rr[[q]]))
        }
    }
    return(list(resdefl=resdefl,pdefl=pdefl))
}


deflation <- function(X, y){
    # Computation of the residual matrix R
    # Computation of the vector p.
    if (any(is.na(X)) || any(is.na(y))) {
        p <- apply(t(X),1,miscrossprod,y)/as.vector(crossprod(y))
    } else {
        p <- t(X)%*%y/as.vector(crossprod(y))
    }
    
    R <- X - y%*%t(p)
    return(list(p=p,R=R))
}

#############################################################################################################
# Functions acquired from gdata R-library
#############################################################################################################
# ---------------------------------------------------
# Column-bind objects with different number of rows
# ---------------------------------------------------
cbindX=function (...)
{
    x <- list(...)
    test <- sapply(x, function(z) is.matrix(z) | is.data.frame(z))
    if (any(!test))
    stop("only matrices and data.frames can be used")
    tmp <- sapply(x, nrow)
    maxi <- which.max(tmp)
    test <- tmp < tmp[maxi]
    for (i in 1:length(tmp)) {
        if (test[i]) {
            add <- matrix(nrow = tmp[maxi] - tmp[i], ncol = ncol(x[[i]]))
            if (is.data.frame(x[[i]])) {
                add <- as.data.frame(add)
            }
            colnames(add) <- colnames(x[[i]])
            x[[i]] <- rbind(x[[i]], add)
        }
    }
    ret <- x[[1]]
    for (i in 2:length(tmp)) {
        ret <- cbind(ret, x[[i]])
    }
    ret
}

#############################################################################################################
# Functions acquired from mixOmics R-library
#############################################################################################################
# ---------------------------------------------------
# unmap variates.A variable for (s)plsda
# ---------------------------------------------------

unmap <- function (classification, groups = NULL, noise = NULL, ...) {
    n <- length(classification)
    u <- sort(unique(classification))
    levels =  levels(classification)### Add levels
    
    if (is.null(groups)) {
        groups <- u
    } else {
        if (any(match(u, groups, nomatch = 0) == 0))
        stop("groups incompatible with classification")
        miss <- match(groups, u, nomatch = 0) == 0
    }
    
    cgroups <- as.character(groups)
    if (!is.null(noise)) {
        noiz <- match(noise, groups, nomatch = 0)
        if (any(noiz == 0))
        stop("noise incompatible with classification")
        groups <- c(groups[groups != noise], groups[groups ==
        noise])
        noise <- as.numeric(factor(as.character(noise), levels = unique(groups)))
    }
    
    groups <- as.numeric(factor(cgroups, levels = unique(cgroups)))
    classification <- as.numeric(factor(as.character(classification), levels = unique(cgroups)))
    k <- length(groups) - length(noise)
    nam <- levels(groups)
    
    if (!is.null(noise)) {
        k <- k + 1
        nam <- nam[1:k]
        nam[k] <- "noise"
    }
    
    z <- matrix(0, n, k, dimnames = c(names(classification), nam))
    for (j in 1:k) z[classification == groups[j], j] <- 1
    attr(z, "levels") = levels
    z
}

map <- function (Y, ...) {
    nrowY <- nrow(Y)
    cl <- numeric(nrowY)
    I <- 1:nrowY
    J <- 1:ncol(Y)
    for (i in I) {
        cl[i] <- (J[Y[i, ] == max(Y[i, ])])[1]
    }
    return(cl)
}

### Estimation of tau accoring to Strimmer formula
tau.estimate <- function (x)
{
    if (is.matrix(x) == TRUE && is.numeric(x) == FALSE)
    stop("The data matrix must be numeric!")
    p <- NCOL(x)
    n <- NROW(x)
    #covm <- cov(x)
    corm <- cor(x)
    xs <- scale(x, center = TRUE, scale = TRUE)
    xs2=xs^2
    v <- (n/((n - 1)^3)) * (crossprod(xs2) - 1/n * (crossprod(xs))^2)
    diag(v) <- 0
    m <- matrix(rep(apply(xs2, 2, mean), p), p, p)
    I <- diag(NCOL(x))
    d <- (corm - I)^2
    tau <- (sum(v))/sum(d)
    tau <- max(min(tau, 1), 0)
    return(tau)
}