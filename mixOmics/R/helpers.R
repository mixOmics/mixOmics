# Author : F.Rohart
# created: pre 04-08-2014
# last modification: 04-08-2014

# -------------------------------------------------------------------------------
# helpers for meta.pls and meta.spls
# -------------------------------------------------------------------------------

# --------------------------------------
# Check.entry.meta.pls
# --------------------------------------

Check.entry.pls = function(X, Y, ncomp, keepX, keepY)
{
    
    
    #  if(length(levels(study)) == 1)  # Aida
    #  stop("\nstudys must have more than one level")      #WHY?
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")
    
    X = as.matrix(X)
    Y = as.matrix(Y)
    
    if (!is.numeric(X) || !is.numeric(Y))
    stop("'X' and/or 'Y' must be a numeric matrix.")
    
    N = nrow(X)
    Q = ncol(Y)
    P= ncol(X)

    if ((N != nrow(Y)))
    stop("unequal number of rows in 'X' and 'Y'.")
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    if (length(keepX) != ncomp)
    stop("length of 'keepX' must be equal to ", ncomp, ".")
    if (any(keepX > ncol(X)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    
    if (length(keepY) != ncomp)
    stop("length of 'keepX' must be equal to ", ncomp, ".")
    if (any(keepY > ncol(Y)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    
    
    
    # if(near.zero.var == TRUE){
    #  nzv = nearZeroVar(X, ...)
    #  if (length(nzv$Position > 0)) {
    #   warning("Zero- or near-zero variance predictors.
    #           Reset predictors matrix to not near-zero variance predictors.
    #           See $nzv for problematic predictors.")
    #   X = X[, -nzv$Position]
    # }
    #}
    
    ncomp = round(ncomp)
    if(ncomp > P)
    {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", P, ".")
        ncomp = P
    }
    
    
    # stop("each component of 'keepY' must be lower or equal than ", q, ".")
    
    #mode = match.arg(mode)
    
    #-- initialisation des matrices --#
    X.names = dimnames(X)[[2]]
    if (is.null(X.names)) X.names = paste("X", 1:P, sep = "")
    
    if (dim(Y)[2] == 1) Y.names = "Y"
    if (dim(Y)[2] > 1)
    {
        Y.names = dimnames(Y)[[2]]
        if (is.null(Y.names)) Y.names = paste("Y", 1:Q, sep = "")
    }
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = dimnames(Y)[[1]]
        rownames(X) = ind.names
    }
    
    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X) = rownames(Y) = ind.names
    }
    
    
}



# --------------------------------------
# study_split
# --------------------------------------
study_split = function(data, study)
{
    
    data = as.matrix(data)
    
    M=length(levels(study))
    P=ncol(data)
    
    #---------------------- split data
    data.list.study = split(data,study)
    study.name = split(rownames(data),study)
    
    for(m in 1:M)
    {
        data.list.study[[m]]=matrix(data.list.study[[m]], ncol=P)
        colnames(data.list.study[[m]]) = colnames(data)
        rownames(data.list.study[[m]]) = study.name[[m]]
    }


    result = list(
    data.list.study = data.list.study)
    
    return(invisible(result))
    
}


# --------------------------------------
# l2.norm
# --------------------------------------
l2.norm=function(x)
{
    if(!is.vector(x)) stop("x has to be a vector")
    out=x/drop(sqrt(crossprod(x)))
}



# --------------------------------------
# soft_thresholding
# --------------------------------------
soft_thresholding=function(x,nx)
{
    if (nx != 0) {
        x = ifelse(abs(x) > abs(x[order(abs(x))][nx]),
        (abs(x) - abs(x[order(abs(x))][nx])) * sign(x), 0)
    }
    x
}


# --------------------------------------
# mean_centering_per_study
# --------------------------------------
mean_centering_per_study=function(data,study,scale)
{
    # initialise some parameters
    M=length(levels(study))   # number of groups
    
    # split the data
    study.data = study_split(data, study)
    data.list.study = study.data$data.list.study
    
    # center and scale data per group, and concatene the data
    concat.data = NULL
    data.list.study.scale= list()
    means.data=sigma.data=matrix(0,nrow=M,ncol=ncol(data))
    
    for (m in 1:M)
    {
        
        data.list.study.scale[[m]] = scale(data.list.study[[m]], center = TRUE, scale = scale)
        means.data[m,]=attr(data.list.study.scale[[m]],"scaled:center")
        if(scale) sigma.data[m,]=attr(data.list.study.scale[[m]],"scaled:scale") else sigma.data[m,]=NA
        
        concat.data = rbind(concat.data, unlist(data.list.study.scale[[m]]))
    }
    
    
    #rename rows and cols of concatenated centered (and/or scaled) data
    colnames(concat.data) = colnames(data)
    
    #sort the samples as in the original X
    indice.match=match(rownames(data),rownames(concat.data))
    concat.data=concat.data[indice.match,]
    
    
    return(list(concat.data=concat.data,means.data=means.data,sigma.data=sigma.data,data.list.study.scale=data.list.study.scale))
    
}

