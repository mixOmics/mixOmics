#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD

# created: 2011
# last modified: 05-10-2017
#
# Copyright (C) 2011
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


# ========================================================================================================
# splsda: perform a sPLS-DA
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# ncomp: the number of components to include in the model. Default to 2.
# keepX: number of \eqn{X} variables kept in the model on the last components.
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# logratio: one of "none", "CLR"
# multilevel: repeated measurement. `multilevel' is passed to multilevel(design = ) in withinVariation. Y is ommited and shouldbe included in `multilevel'
# all.outputs: calculation of non-essential outputs (e.g. explained variance, loadings.Astar, etc)

Check.entry.single = function(X,  ncomp, q)
{
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop(paste0("'X[[", q, "]]' must be a numeric matrix."))
    
    X = as.matrix(X)
    
    if (!is.numeric(X))
    stop(paste0("'X[[", q, "]]'  must be a numeric matrix."))
    
    N = nrow(X)
    P = ncol(X)
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop(paste0("invalid number of variates 'ncomp' for matrix 'X[[", q, "]]'."))
    
    ncomp = round(ncomp)
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]] = X.names
    }
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X)  = ind.names
    }
    
    if (length(unique(rownames(X))) != nrow(X))
    stop("samples should have a unique identifier/rowname")
    if (length(unique(X.names)) != P)
    stop("Unique indentifier is needed for the columns of X")
    
    return(list(X=X, ncomp=ncomp, X.names=X.names, ind.names=ind.names))
}


explained_variance = function(data, variates, ncomp)
{
    #check input data
    check = Check.entry.single(data, ncomp)
    data = check$X
    ncomp = check$ncomp
    
    isna = is.na(data)
    if (sum(isna > 0))
    {
        warning("NA values put to zero, results will differ from PCA methods used with NIPALS")
        data[isna] = 0
    }
    nor2x <- sum((data)^2) # total variance in the data
    
    exp.varX = NULL
    for (h in 1:ncomp)
    {
        a <- t(variates[, h, drop=FALSE]) %*% data
        ta = t(a)
        exp_var_new <- a%*%ta /crossprod(variates[, h],variates[, h])/nor2x
        
        
        exp.varX = append(exp.varX, exp_var_new)
        
    }
    names(exp.varX) = paste("comp", 1:ncomp)
    
    # result: vector of length ncomp with the explained variance per component
    exp.varX
}



unmap = function (classification, groups = NULL, noise = NULL)
{
    n = length(classification)
    u = sort(unique(classification))
    levels =  levels(classification)### Add levels
    
    if (is.null(groups))
    {
        groups = u
    } else {
        if (any(match(u, groups, nomatch = 0) == 0))
        stop("groups incompatible with classification")
        miss = match(groups, u, nomatch = 0) == 0
    }
    
    cgroups = as.character(groups)
    if (!is.null(noise))
    {
        noiz = match(noise, groups, nomatch = 0)
        if (any(noiz == 0))
        stop("noise incompatible with classification")
        
        groups = c(groups[groups != noise], groups[groups == noise])
        noise = as.numeric(factor(as.character(noise), levels = unique(groups)))
    }
    
    groups = as.numeric(factor(cgroups, levels = unique(cgroups)))
    classification = as.numeric(factor(as.character(classification), levels = unique(cgroups)))
    k = length(groups) - length(noise)
    nam = levels(groups)
    
    if (!is.null(noise))
    {
        k = k + 1
        nam = nam[1:k]
        nam[k] = "noise"
    }
    
    z = matrix(0, n, k, dimnames = c(names(classification), nam))
    for (j in 1:k) z[classification == groups[j], j] = 1
    attr(z, "levels") = levels
    z
}

Check.entry.pls = function(X, Y, ncomp, keepX, keepY, test.keepX, test.keepY, mode, scale,
near.zero.var, max.iter, tol, logratio, DA, multilevel)
{
    
    if (missing(mode))
    mode = "regression"
    
    if (length(mode)>1)
    mode = mode[1]
    
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
    stop("Choose one of the four following modes: canonical, invariant, classic or regression")
    
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix.")
    
    X = as.matrix(X)
    
    if (!(logratio %in% c("none", "CLR")))
    stop("Choose one of the two following logratio transformation: none or CLR")
    
    if(!is.null(multilevel))
    {
        #multilevel analysis: withinVariation and then pls-like
        # if it's DA analysis, Y and 'multilevel' are combined
        if(DA)
        {
            Y = multilevel
        }else{
            if ((nrow(X) != nrow(multilevel)))
            stop("unequal number of rows in 'X' and 'multilevel'.")
            
            Y = as.matrix(Y)
            if (!is.numeric(X) || !is.numeric(Y))
            stop("'X' and/or 'Y' must be a numeric matrix.")
        }
    }else{
        Y = as.matrix(Y)
        if (!is.numeric(X) || !is.numeric(Y))
        stop("'X' and/or 'Y' must be a numeric matrix.")
    }
    N = nrow(X)
    Q = ncol(Y)
    P= ncol(X)
    
    if ((N != nrow(Y)))
    stop("Unequal number of rows in 'X' and 'Y'.")
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || length(ncomp)>1)
    stop("invalid number of variates, 'ncomp'.")
    
    if(mode == "canonical" & ncomp>ncol(Y))
    stop("For `canonical mode', 'ncomp' needs to be lower than ncol(Y)= ",ncol(Y))
    
    
    ncomp = round(ncomp)
    if(ncomp > P)
    {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", P, ".")
        ncomp = P
    }
    
    if (!is.numeric(tol) | tol<=0)
    stop("tol must be non negative")
    
    if (!is.numeric(max.iter) | max.iter<=0)
    stop("max.iter must be non negative")
    
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]] = X.names
    }
    
    
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = dimnames(Y)[[1]]
        #rownames(X) = ind.names
    }
    
    if (is.null(ind.names))
    {
        ind.names = 1:N
        #rownames(X) = rownames(Y) = ind.names
    }
    
    rownames(X) = rownames(Y) = ind.names
    
    
    #if (dim(Y)[2] == 1) Y.names = "Y"
    Y.names = dimnames(Y)[[2]]
    if (is.null(Y.names))
    {
        if (dim(Y)[2] == 1)
        {
            Y.names = "Y"
        } else {
            Y.names = paste("Y", 1:Q, sep = "")
        }
        
        dimnames(Y)[[2]]=Y.names
    }
    
    if (length(unique(X.names)) != P)
    stop("Unique indentifier is needed for the columns of X")
    
    if (length(unique(Y.names)) != Q)
    stop("Unique indentifier is needed for the columns of Y")
    
    
    # check keepX
    if (missing(keepX))
    {
        keepX = rep(P, ncomp)
    } else {
        if (length(keepX)<ncomp)
        keepX = c(keepX, rep(P, ncomp - length(keepX))) #complete (with ncomp) the keepX already provided
    }
    
    # check keepY
    if (missing(keepY))
    {
        keepY = rep(Q, ncomp)
    } else {
        if (length(keepY) < ncomp)
        keepY = c(keepY, rep(Q, ncomp - length(keepY))) #complete the keepY already provided
    }
    
    
    
    
    if (any(keepX<0))
    stop("each component of 'keepX' must be non negative ")
    if (any(keepY<0))
    stop("each component of 'keepY' must be non negative ")
    
    if (any(keepX > ncol(X)))
    stop("each component of 'keepX' must be lower or equal than ", P, ".")
    if (any(keepY > ncol(Y)))
    stop("each component of 'keepY' must be lower or equal than ", Q, ".")
    
    
    if (!is.logical(scale))
    stop("'scale' must be either TRUE or FALSE")
    
    if (!is.logical(near.zero.var))
    stop("'near.zero.var' must be either TRUE or FALSE")
    
    
    ### near.zero.var, remove the variables with very small variances
    if (near.zero.var == TRUE)
    {
        nzv.A = nearZeroVar(X)
        
        if (length(nzv.A$Position) > 0)
        {
            names.remove.X = colnames(X)[nzv.A$Position]
            X = X[, -nzv.A$Position, drop=FALSE]
            warning("Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors.")
            if (ncol(X) == 0)
            stop("No more variables in X")
            
            #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
            if (any(keepX > ncol(X)))
            {
                ind = which(keepX > ncol(X))
                keepX[ind] = ncol(X)
            }
        }
        
    }else{nzv.A=NULL}
    
    
    return(list(X=X, Y=Y, ncomp=ncomp, X.names=X.names, Y.names=Y.names, ind.names=ind.names, mode=mode,
    keepX=keepX, keepY=keepY, nzv.A=nzv.A))
}


splsda_allin = function(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
keepX,
scale = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
logratio = "none",   # one of "none", "CLR"
multilevel = NULL,
all.outputs = TRUE)    
{
    
    
    #-- validation des arguments --#
    # most of the checks are done in the wrapper.mint.spls.hybrid function
    if(is.null(multilevel))
    {
        if (is.null(Y))
        stop("'Y' has to be something else than NULL.")
        
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
        stop("'Y' should be a factor with more than one level")

        Y.mat = unmap(Y)
        colnames(Y.mat) = levels(Y)#paste0("Y", 1:ncol(Y.mat))
    }else{
        # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
        multilevel = data.frame(multilevel)
        
        if ((nrow(X) != nrow(multilevel)))
        stop("unequal number of rows in 'X' and 'multilevel'.")
        
        if (ncol(multilevel) != 1)
        stop("'multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.")
        
        if (!is.null(ncol(Y)) && !ncol(Y) %in% c(0,1,2))# multilevel 1 or 2 factors
        stop("'Y' should either be a factor, a single column data.frame containing a factor, or a 2-columns data.frame containing 2 factors.")
        
        multilevel = data.frame(multilevel, Y)
        multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements
        
        Y.mat = NULL
    }
    
    # call to 'internal_wrapper.mint'
DA=TRUE

#-- validation des arguments --#
test.keepX=NULL
test.keepY=NULL

check = Check.entry.pls(X, Y.mat, ncomp, keepX, mode=mode, scale=scale,
near.zero.var=near.zero.var, max.iter=max.iter ,tol=tol ,logratio=logratio ,DA=DA, multilevel=multilevel)
X = check$X
input.X = X # save the checked X, before logratio/multileve/scale
Y = check$Y
ncomp = check$ncomp
mode = check$mode
keepX = check$keepX
keepY = check$keepY
nzv.A = check$nzv.A

rm(check) # free memory

#test.keepX and test.keepY must be checked before (in tune)

#set the default study factor
   study = factor(rep(1,nrow(X)))

if (length(study) != nrow(X))
stop(paste0("'study' must be a factor of length ",nrow(X),"."))

if (any(table(study) <= 1))
stop("At least one study has only one sample, please consider removing before calling the function again")
if (any(table(study) < 5))
warning("At least one study has less than 5 samples, mean centering might not do as expected")

design = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)


#-----------------------------#
#-- logratio transformation --#

#X = logratio.transfo(X=X, logratio=logratio)

#as X may have changed
if (ncomp > min(ncol(X), nrow(X)))
stop("'ncomp' should be smaller than ", min(ncol(X), nrow(X)), call. = FALSE)

#-- logratio transformation --#
#-----------------------------#


#---------------------------------------------------------------------------#
#-- multilevel approach ----------------------------------------------------#

if (!is.null(multilevel))
{
    if (!DA)
    {
        Xw = withinVariation(X, design = multilevel)
        Yw = withinVariation(Y, design = multilevel)
        X = Xw
        Y = Yw
    } else {
        Xw = withinVariation(X, design = multilevel)
        X = Xw
        
        #-- Need to set Y variable for 1 or 2 factors
        Y = multilevel[, -1,drop=FALSE]
        if (ncol(Y)>0)
        Y = apply(Y, 1, paste, collapse = ".")  #  paste is to combine in the case we have 2 levels
        
        Y = as.factor(Y)
        Y.factor = Y
        Y = unmap(Y)
        colnames(Y) = levels(Y)
        rownames(Y) = rownames(X)
        # if DA keepY should be all the levels (which is not happening in the check because of multilevel
        keepY = rep(ncol(Y),ncomp)
    }
}
#-- multilevel approach ----------------------------------------------------#
#---------------------------------------------------------------------------#


#---------------------------------------------------------------------------#
#-- keepA ----------------------------------------------------#

# shaping keepA, contains all the keepX/keepY models to be constructed

if(!is.null(test.keepX) & !is.null(test.keepY))
{
    test.keepA = lapply(list(X=test.keepX, Y=test.keepY),sort) #sort test.keepX so as to be sure to chose the smallest in case of several minimum
} else {test.keepA=NULL}

keepA = vector("list", length = ncomp) # one keepA per comp
names(keepA) = paste0("comp",1:ncomp)
for(comp in 1:length(keepX)) # keepA[[block]] [1:ncomp]
keepA[[comp]] = lapply(list(X=keepX, Y=keepY), function(x) x[comp])

if(!is.null(test.keepA))
keepA[[ncomp]] = test.keepA

keepA = lapply(keepA, expand.grid)

# keepA[[comp]] is a matrix where each row is all the keepX the test over the block (each block is a column)

#-- keepA ----------------------------------------------------#
#---------------------------------------------------------------------------#

#print(keepA)

A = list(X = X, Y = Y)
rm(X)
rm(Y)
indY = 2
mode = mode
ncomp = c(ncomp, ncomp)
scheme = "horst"
init="svd"
misdata = NULL
is.na.A = NULL
ind.NA = NULL
ind.NA.col = NULL
tau=NULL
penalty=NULL

names(ncomp) = names(A)

time = FALSE

if(time) time1=proc.time()

# center the data per study, per matrix of A, scale if scale=TRUE, option
mean_centered = lapply(A, function(x){mean_centering_per_study(x, study, scale)})
if(time) time1bis=proc.time()
if(time) print("scaling part1")
if(time) print(time1bis-time1)

A = lapply(mean_centered, function(x){as.matrix(x$concat.data)})

#save rownames study
mean_centered.rownames.study = vector("list", nlevels(study))
for (m in 1:nlevels(study))
mean_centered.rownames.study[[m]] = mean_centered[[1]]$rownames.study[[m]]

rm(mean_centered) #free memory

ni = table(study) #number of samples per study

### Start: Initialization parameters
pjs = sapply(A, NCOL)
nb_ind = NROW(A[[1]])
J = length(A)
R = A # R: residuals matrices, will be a list of length ncomp
N = max(ncomp)
AVE_inner = AVE_outer = rep(NA, max(ncomp))


# keepA[[comp]] is a matrix where each row is all the keepX the test over the block (each block is a column)

#number of models to be tested: either a keepA per component, or multiple (e.g. in tune functions)
number.models.per.comp = sapply(keepA, nrow)
one.model = !any( number.models.per.comp !=1)

if(time) print("one.model")
if(time) print(one.model)

AVE_X = crit = loadings.partial.A = variates.partial.A = tau.rgcca = list()
P = loadings.A = loadings.Astar = variates.A =  vector("list", J)


if(one.model) # more outputs that what is needed for tune functions
{
    for (k in 1:J)
    variates.A[[k]] = matrix(NA_real_, nb_ind, N)
    
    for (k in 1:J)
    {
        loadings.A[[k]] = matrix(NA_real_, pjs[[k]], N)
        if(all.outputs)
        P[[k]] = loadings.Astar[[k]]= matrix(NA_real_, pjs[[k]], N)
    }
    
    for (k in 1:J)
    {
        loadings.partial.A[[k]] = variates.partial.A[[k]] = vector("list", length = nlevels(study))
        for(m in 1:nlevels(study))
        {
            loadings.partial.A[[k]][[m]] = matrix(NA_real_,nrow = NCOL(A[[k]]), ncol = N)
            variates.partial.A[[k]][[m]] = matrix(NA_real_,nrow = ni[m], ncol = N)
        }
    }
} else {
    for (k in 1:J)
    {
        variates.A[[k]] = matrix(NA_real_, nb_ind, sum(number.models.per.comp))
        loadings.A[[k]] = matrix(NA_real_, pjs[[k]], sum(number.models.per.comp))
    }
    loadings.partial.A = variates.partial.A = NULL # not needed for tune functions
}


ndefl = ncomp - 1
J2 = J-1

if (is.vector(tau))
tau = matrix(rep(tau, N), nrow = N, ncol = length(tau), byrow = TRUE)

#save(list=ls(),file="temp.Rdata")

# if missing values are not given as input (only when direct call to a (mint).(block).(s)pls(da)), we search for them here (takes time)
if(is.null(misdata) &  is.null(is.na.A) & is.null(ind.NA) & is.null(ind.NA.col))
{
    misdata = sapply(A, anyNA) # Detection of missing data per block
    misdata.all = any(misdata) # is there any missing data overall
    
    if(time) print(misdata.all)
    
    if (misdata.all)
    {
        is.na.A = lapply(A, is.na) # size n*p, which entry is na. might be none, but at least one in all the block will be a TRUE
        
        ind.NA = ind.NA.col = list()
        for(q in 1:J)
        {
            ind.NA[[q]] = which(apply(is.na.A[[q]], 1, sum) >0) # indice of the row that have missing values. used in the calculation of loadings
            ind.NA.col[[q]] = which(apply(is.na.A[[q]], 2, sum) >0) # indice of the col that have missing values. used in the deflation
        }
    }else {
        is.na.A = NULL
        ind.NA = ind.NA.col = NULL
    }
} else{
    misdata.all = any(misdata)
}

if(time) print(misdata)

if(all.outputs & J==2 & nlevels(study) == 1 & one.model) #(s)pls(da) models, we calculate mat.c
{
    if(misdata.all)
    {
        p.ones = rep(1, ncol(A[[1]]))
        is.na.X = is.na.A[[1]]
    }
    mat.c = matrix(0, nrow = ncol(A[[1]]), ncol = N, dimnames = list(colnames(A[[1]],  paste0("comp ", 1:N))))
} else {mat.c = NULL}

### End: Initialization parameters

iter=NULL
compteur = 0
for (comp in 1 : N)
{
    if(time) cat("====================================================================================================\n")
    if(time) print(paste0("component ", comp))
    if(time) cat("====================================================================================================\n")
    
    if(misdata.all)# replace NA in A[[q]] by 0
    R = lapply(1:J, function(q){replace(R[[q]], is.na.A[[q]], 0)}) # if missing data, R is the one replace by 0 where NA are supposed to be
    
    # initialisation_by_svd, get the loadings.A
    loadings.A.init = initialisation_by_svd(R, indY, misdata, is.na.A, init = init)
    
    # loop on keepA[[comp]]: multiple values per block and we go through them. Need to have the same number of values per block.
    # we assume keepA[[comp]] is a grid here: columns are the blocks, rows are the different keepX
    for(ijk.keepA in 1:nrow(keepA[[comp]]))
    {
        compteur = compteur +1
        if(time) cat("---------------------------------------------------------\n")
        if(time) print(paste0("keepA ", ijk.keepA))
        if(time) cat("---------------------------------------------------------\n")
        
        keepA.ijk = keepA[[comp]][ijk.keepA,]
        
        ### start - repeat/convergence
        if (is.null(tau))
        {
            loadings.R = loadings.A.init
            
            
            
            
            time=FALSE
            
            ### Start: Initialization parameters
            J = length(R)
            J2 = J-1
            pjs = sapply(R, NCOL)
            AVE_X = rep(0, J)
            if (!is.null(penalty))
            penalty = penalty * sqrt(pjs)
            
            iter = 1
            converg = crit = numeric()
            variates.R = Z = matrix(0, NROW(R[[1]]), J)
            
            g = function(x) switch(scheme, horst = x, factorial = x^2, centroid = abs(x))
            
            # study split
            R_split = lapply(R, study_split, study)
            
            n = lapply(R_split, function(x){lapply(x,nrow)})
            p = lapply(R,ncol)
            
            nlevels_study = nlevels(study)
            ### End: Initialization parameters
            
            if(time) time2 = proc.time()
            
            ### End: Initialisation "a" vector
            variates.partial.R.comp = NULL
            loadings.partial.R.comp = list()
            for (q in 1:J)
            {
                if(misdata[q])
                {
                    loadings.temp = loadings.R[[q]]
                    variates.R.temp = R[[q]] %*% loadings.temp
                    
                    # we only want the diagonal, which is the norm of each column of temp
                    # loadings.R.norm = crossprod(temp)
                    # variates.R[, q] = variates.R.temp / diag(loadings.R.norm)
                    #only calculating the ones where there's a NR
                    d.variates.R.norm = rep(crossprod(loadings.temp), length(variates.R.temp))
                    
                    if(length(ind.NA[[q]])>0) # should always be true
                    {
                        temp = drop(loadings.temp) %o% rep(1, length(ind.NA[[q]])) #p*n -> p * where there are NA
                        temp[t(is.na.A[[q]][ind.NA[[q]],,drop=FALSE])] = 0
                        d.variates.R.norm[ind.NA[[q]]] = apply(temp,2, crossprod)
                    }
                    
                    variates.R[, q] = variates.R.temp / d.variates.R.norm
                    
                    # we can have 0/0, so we put 0
                    a = is.na(variates.R[, q])
                    if (any(a))
                    variates.R[a, q] = 0
                    
                }else{
                    variates.R[, q] = R[[q]]%*%loadings.R[[q]]
                }
                loadings.R[[q]] = l2.norm(as.vector(loadings.R[[q]]))
                loadings.partial.R.comp[[q]] = list()
            }
            loadings.R_old = loadings.R
            
            if(time) time3 = proc.time()
            if(time) print("loadings")
            if(time) print(time3-time2)
            if(time) print("repeat")
            
            ### Start Algorithm 1 Sparse generalized canonical analysis (See Variable selection for generalized canonical correlation analysis (Tenenhaus))
            repeat {
                #variates.Aold can be used for the convergence of the algorithm, not used at the moment, so please disregard variates.Aold
                variates.Rold = variates.R
                for (q in 1:J)
                {
                    if(time) print(paste("repeat, block",q))
                    
                    ### Start : !!! Impact of the diag of the design matrix !!! ###
                    if (scheme == "horst")
                    CbyCovq = design[q, ]
                    
                    if (scheme == "factorial")
                    CbyCovq = design[q, ] * cov2(variates.R, variates.R[, q])
                    
                    if (scheme == "centroid")
                    CbyCovq = design[q, ] * sign(cov2(variates.R, variates.R[, q]))
                    ### End : !!! Impact of the diag of the design matrix !!! ###
                    
                    ### Step R start: Compute the inner components
                    Z[, q] = rowSums(mapply("*", CbyCovq, as.data.frame(variates.R)))
                    Z_split = study_split(Z[,q,drop=FALSE],study)  # split Z by the study factor
                    ### Step A end: Compute the inner components
                    
                    
                    if(time) time5 = proc.time()
                    
                    ### Step B start: Computer the outer weight ###
                    # possibility of removing NA (replacing by 0) and use crossprod, further development
                    #time5 = proc.time()
                    temp=0
                    for (m in 1:nlevels_study)
                    {
                        loadings.partial.R.comp[[q]][[m]] = crossprod(R_split[[q]][[m]],Z_split[[m]])
                        temp=temp+loadings.partial.R.comp[[q]][[m]]
                    }
                    loadings.R[[q]] = temp
                    
                    if(time) time6 = proc.time()
                    if(time) print("loadings")
                    if(time) print(time6-time5)
                    if(time) time6 = proc.time()
                    
                    
                    # sparse using keepA / penalty
                    if (!is.null(penalty))
                    {
                        loadings.R[[q]] = sparsity(loadings.R[[q]], keepA = NULL, penalty = penalty[q])
                    }else{
                        loadings.R[[q]] = sparsity(loadings.R[[q]], keepA.ijk[[q]], penalty = NULL)
                    }
                    
                    loadings.R[[q]]=l2.norm(as.vector(loadings.R[[q]]))
                    
                    if(time) time7 = proc.time()
                    
                    ### Step B end: Computer the outer weight ###
                    if(misdata[q])
                    {
                        variates.R.temp = R[[q]] %*% loadings.R[[q]]
                        d.variates.R.norm = rep(crossprod(loadings.R[[q]]), length(variates.R.temp))
                        
                        if(length(ind.NA[[q]])>0)
                        {
                            temp = drop(loadings.R[[q]]) %o% rep(1, length(ind.NA[[q]]))
                            temp[t(is.na.A[[q]][ind.NA[[q]],,drop=FALSE])] = 0
                            d.variates.R.norm[ind.NA[[q]]] = apply(temp,2, crossprod)
                        }
                        variates.R[, q] = variates.R.temp / d.variates.R.norm
                        
                        # we can have 0/0, so we put 0
                        a = is.na(variates.R[, q])
                        if (any(a))
                        variates.R[a, q] = 0
                        
                    }else{
                        variates.R[, q] =  R[[q]]%*%loadings.R[[q]]
                    }
                    
                    if(time) time8 = proc.time()
                    if(time) print("variates")
                    if(time) print(time8-time7)
                    
                }
                
                crit[iter] = sum(design * g(cov2(variates.R)))
                
                if (iter > max.iter)
                warning("The SGCCA algorithm did not converge", call. = FALSE)#cat("The SGCCA algorithm did not converge after", max.iter ,"iterations."))
                
                ### Start: Match algorithm with mixOmics algo (stopping point)
                ### if ((converg[iter] < tol & sum(stationnary_point) == J) | iter > max.iter)
                if (max(sapply(1:J, function(x){crossprod(loadings.R[[x]] - loadings.R_old[[x]])})) < tol | iter > max.iter)
                break
                ### End: Match algorithm with mixOmics algo (stopping point)
                
                loadings.R_old = loadings.R
                iter = iter + 1
            }
            ### End Algorithm 1 (See Variable selection for generalized canonical correlation analysis (Tenenhaus))
            
            #calculation variates.partial.A.comp
            variates.partial.R.comp = apply(variates.R, 2, study_split, study)
            
            if(all.outputs){
                AVE_inner = sum(design * cor(variates.R)^2/2)/(sum(design)/2)
            } else{
                AVE_inner = NULL
            }
            
            names(loadings.R) = colnames(variates.R) = names(variates.partial.R.comp) = names(R)
            
            mint.block.result = list(variates.A = variates.R, loadings.A = loadings.R, crit = crit[which(crit != 0)],
            AVE_inner = AVE_inner, loadings.partial.A.comp = loadings.partial.R.comp, variates.partial.A.comp = variates.partial.R.comp, iter = iter)
            
        } else {
            mint.block.result = sparse.rgcca_iteration(R, design, tau = if (is.matrix(tau)){tau[comp, ]} else {"optimal"},
            scheme = scheme, init = init, tol = tol,
            max.iter = max.iter, penalty = penalty,
            keepA = keepA.ijk, all.outputs = all.outputs)
        }
        ### end - repeat/convergence
        
        if(time) time3 = proc.time()
        if(time) time4 = proc.time()
        if(time) print("one comp")
        if(time) print(time4-time3)
        if(time) time4 = proc.time()
        
        if(one.model)
        {
            # reshape outputs
            for (k in 1 : J)
            {
                loadings.A[[k]][, comp] = mint.block.result$loadings.A[[k]]
                variates.A[[k]][, comp] = mint.block.result$variates.A[, k]
                
                if(is.null(tau))
                {
                    # recording loadings.partials, $Ai$study[,ncomp]
                    # recording variates.partials, $Ai[,ncomp]
                    for(k in 1:J)
                    {
                        for(m in 1:nlevels(study))
                        {
                            loadings.partial.A[[k]][[m]][, comp] = matrix(mint.block.result$loadings.partial.A.comp[[k]][[m]], ncol=1)
                            variates.partial.A[[k]][[m]][, comp] = matrix(mint.block.result$variates.partial.A.comp[[k]][[m]], ncol=1)
                        }
                    }
                }
            }
        } else {
            # no record of partial component for multilple models, for gain of memory
            for (k in 1 : J)
            {
                loadings.A[[k]][, compteur] = mint.block.result$loadings.A[[k]]
                variates.A[[k]][, compteur] = mint.block.result$variates.A[, k]
            }
        }
        
        crit[[comp]] = mint.block.result$crit
        tau.rgcca[[comp]] = mint.block.result$tau
        if(all.outputs)
        AVE_inner[comp] = mint.block.result$AVE_inner
        
        if(all.outputs & J==2 & nlevels(study) == 1 & one.model)# mat.c, (s)pls(da)
        {
            if(misdata.all)
            {
                R.temp = R[[1]]
                R.temp[is.na.X] = 0
                c = crossprod(R.temp, variates.A[[1]][,comp])
                rm(R.temp) #free memory
                T = drop(variates.A[[1]][,comp]) %o% p.ones
                T[is.na.X] = 0
                t.norm = crossprod(T)
                c = c / diag(t.norm)
                mat.c[,comp] = c
            } else {
                mat.c[,comp] <- t(crossprod(variates.A[[1]][,comp], R[[1]])) / drop(crossprod (variates.A[[1]][,comp]))
            }
        } else {
            mat.c = NULL
        }
        
        # deflation if there are more than 1 component and if we haven't reach the max number of component (N)
        if (N != 1 & comp != N)
        {
            if(time) time4 = proc.time()
            
            defla.result = defl.select(yy=mint.block.result$variates.A, rr=R, nncomp=ndefl, nn=comp, nbloc = J, indY = indY, mode = mode, aa = mint.block.result$loadings.A,
            misdata=misdata, is.na.A=is.na.A, ind.NA = ind.NA.col)
            
            R = defla.result$resdefl
            
            if(!(all.outputs & one.model)) #rm only if not used in the loop below
            defla.result$resdefl=NULL #free memory
            
            if(time) time5 = proc.time()
            if(time) print("deflation")
            if(time) print(time5-time4)
        }
        
        
        if(all.outputs & one.model) #loadings.Astar
        {
            for (k in 1 : J)
            {
                if (N != 1)
                P[[k]][, comp - 1] = defla.result$pdefl[[k]]
            }
            
            if (comp == 1)
            {
                for (k in 1 : J)
                loadings.Astar[[k]][, comp] = mint.block.result$loadings.A[[k]]
            } else {
                for (k in 1 : J)
                loadings.Astar[[k]][, comp] = mint.block.result$loadings.A[[k]] - loadings.Astar[[k]][, (1 : comp - 1), drop = F] %*% drop(t(loadings.A[[k]][, comp]) %*% P[[k]][, 1 : (comp - 1), drop = F])
            }
        } else {
            loadings.Astar = NULL
        }
        iter = c(iter, mint.block.result$iter)
        
    } ### End loop on keepA
} ### End loop on ncomp


#### any model
# loadings.A[[block]][1:p, all.keepA.tested]
# variates.A[[block]][1:n, all.keepA.tested]

#### a single model
# loadings.partial.A[[block]][[study]][, 1:ncomp]
# variates.partial.A[[block]][[study]][, 1:ncomp]
# loadings.Astar[[block]][, 1:ncomp]

if(one.model)
{
    # only one model
    shave.matlist = function(mat_list, nb_cols) mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY = FALSE)
    shave.veclist = function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
    
    
    for (k in 1:J)
    {
        rownames(loadings.A[[k]]) = colnames(A[[k]])
        
        if(all.outputs)
        rownames(loadings.Astar[[k]]) = colnames(A[[k]])
        
        rownames(variates.A[[k]]) = rownames(A[[k]])
        colnames(variates.A[[k]]) = colnames(loadings.A[[k]]) = paste0("comp ", 1:max(ncomp))
        if(all.outputs)
        #AVE_X[[k]] = apply(cor(A[[k]], variates.A[[k]])^2, 2, mean)
        
        if (is.null(tau))
        {
            names(loadings.partial.A[[k]]) = names(variates.partial.A[[k]]) = levels(study)
            
            for (m in 1:nlevels(study))
            {
                rownames(loadings.partial.A[[k]][[m]]) = colnames(A[[k]])
                colnames(loadings.partial.A[[k]][[m]]) = paste0("comp ", 1:max(ncomp))
                rownames(variates.partial.A[[k]][[m]]) = mean_centered.rownames.study[[m]]
                colnames(variates.partial.A[[k]][[m]]) = paste0("comp ", 1:max(ncomp))
            }
        }
    }
    
    variates.A = shave.matlist(variates.A, ncomp)
    
    if(all.outputs)
    {
        # AVE
        outer = matrix(unlist(AVE_X), nrow = max(ncomp))
        #for (j in 1 : max(ncomp))
        #AVE_outer[j] = sum(pjs * outer[j, ])/sum(pjs)
        #AVE_X = shave.veclist(AVE_X, ncomp)
        #AVE = list(AVE_X = AVE_X, AVE_outer = AVE_outer, AVE_inner = AVE_inner)
        #names(AVE$AVE_X) = names(A)
        
        loadings.Astar = shave.matlist(loadings.Astar, ncomp)
        
        #calcul explained variance
        A_split=lapply(A, study_split, study) #split the data per study
        expl.A=lapply(1:length(A),function(x){
            if (nlevels(study) == 1)
            {
                temp = suppressWarnings(explained_variance(A[[x]], variates = variates.A[[x]], ncomp = ncomp[[x]]))
            }else{
                temp = lapply(1:nlevels(study), function(y){
                    suppressWarnings(explained_variance(A_split[[x]][[y]], variates = variates.partial.A[[x]][[y]], ncomp = ncomp[[x]]))})
                temp[[length(temp)+1]] = explained_variance(A[[x]], variates = variates.A[[x]], ncomp = ncomp[[x]])
                names(temp) = c(levels(study), "all data")
            }
            temp
        })
        names(expl.A) = names(A)
    } else {
        expl.A = NULL
        AVE = NULL
    }
    ### Start: output
    names(loadings.A) = names(variates.A) = names(A)
    
    if (is.null(tau))
    names(loadings.partial.A) = names(variates.partial.A) = names(A)
    
    names = lapply(1:J, function(x) {colnames(A[[x]])})
    names(names) = names(A)
    names[[length(names) + 1]] = row.names(A[[1]])
    names(names)[length(names)] = "indiv"
} else {
    # multiple models (tune)
    
    #### any model
    # loadings.A[[block]][1:p, all.keepA.tested]
    # variates.A[[block]][1:n, all.keepA.tested]
    
    #### a single model
    # loadings.partial.A[[block]][[study]][, 1:ncomp]
    # variates.partial.A[[block]][[study]][, 1:ncomp]
    # loadings.Astar[[block]][, 1:ncomp]
    
    
    keepA.names = unlist(lapply(1:N, function(x){
        paste(paste0("comp",x),apply(keepA[[x]],1,function(x) paste(x,collapse="_")), sep=":")
        
    }))
    
    for(k in 1:J)
    colnames(loadings.A[[k]]) = colnames(variates.A[[k]]) = keepA.names
    
    names(loadings.A) =  names(variates.A) = names(iter) = names(A)
    
    expl.A = NULL
    AVE = NULL
    
}

result = list(A = A, indY = indY, ncomp = ncomp, mode = mode,
keepA = keepA,
variates = variates.A, loadings = loadings.A,#shave.matlist(loadings.A, ncomp),
variates.partial= if(is.null(tau)) {variates.partial.A} ,loadings.partial= if(is.null(tau)) {loadings.partial.A},
loadings.star = loadings.Astar,
names = list(sample = row.names(A[[1]]), colnames = lapply(A, colnames), blocks = names(A)),
tol = tol, iter=iter, max.iter=max.iter,
design = design,
scheme = scheme,  crit = crit, mat.c = mat.c, #defl.matrix = defl.matrix,
init = init,
scale = scale, tau = if(!is.null(tau)) tau.rgcca, study = study,
explained_variance = expl.A)
### End: Output


# choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = if (is.null(multilevel))
            {
                Y
            } else {
                result$Y.factor
            },
        ind.mat = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        keepA=result$keepA,
        keepX = result$keepX,
        keepY = result$keepY,
        variates = result$variates,
        loadings = result$loadings,
        loadings.star = result$loadings.star,
        names = result$names,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = scale,
        logratio = logratio,
        explained_variance = result$explained_variance,#[-result$indY],
        input.X = result$input.X,
        mat.c = result$mat.c#,
        )
    
    class(out) = c("splsda","spls","DA")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        out$multilevel = multilevel
        class(out) = c("mlsplsda",class(out))
    }

    return(invisible(out))
}

