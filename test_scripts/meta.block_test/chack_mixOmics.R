###############################
############# PLS #############
###############################
#### Difference because of not the same starting point between the 2 algo "pls" and "spls" in mixOmics package when tol = 1e-06
rm(list=ls())
setwd("/Users/florian/Work/git/package-mixOmics/")
load("test_scripts/multigroup_test/Fibro-ESC-iPSC.6exp.167samples.light.Rdata") #load data, type.id, exp among others


library(mixOmics)


source("meta.block.spls/check_entry.R")
source("meta.block.spls/helpers.R")
source("meta.block.spls/meta.block.spls.R")
source("meta.block.spls/meta.spls.hybrid.R")
source("meta.block.spls/mixOmics.R")
source("meta.block.spls/pls.R")
source("meta.block.spls/plsda.R")
source("meta.block.spls/spls.R")
source("meta.block.spls/splsda.R")



#res.sgcca=meta.spls.hybrid.sgcca(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,5),ncomp=3,near.zero.var=FALSE,tol=1e-25)
Y.mat=unmap(type.id)
rownames(Y.mat)=rownames(data)

#wraper.meta.spls.hybrid
res=wrapper.meta.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,tol=1e-25)

## ======================================================================
## =========================        data for one study
## ======================================================================
ind=which(exp=="Briggs")
study.light=exp[ind]
type.id.light=type.id[ind]
data.light=data[ind,]
Y.mat.light=unmap(type.id.light)
rownames(Y.mat.light)=rownames(data.light)

#wraper.pls
res=wrapper.pls(X=data.light,Y=Y.mat.light,ncomp=3,near.zero.var=FALSE)


#wraper.plsda
res=wrapper.plsda(X=data.light,Y=type.id.light,ncomp=3,near.zero.var=FALSE)


#wraper.spls
res=wrapper.spls(X=data.light,Y=Y.mat.light,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15))


#wraper.splsda
res=wrapper.splsda(X=data.light,Y=type.id.light,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15))





res=mixOmics(data.light,type.id.light,ncomp=3) #plsda
res=mixOmics(data.light,type.id.light,ncomp=3,keepX=c(10,5,15))#splsda
res=mixOmics(data.light,Y.mat.light,ncomp=3) #pls
res=mixOmics(data.light,Y.mat.light,ncomp=3,keepX=c(10,5,15)) #spls


res=mixOmics(data,type.id,ncomp=3,study=exp) #pls
res=mixOmics(X=data,Y=type.id,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,15),ncomp=3)

## ======================================================================
## =========================        checks for wrong inputs
## ======================================================================

#Y too long
res=mixOmics(data,type.id.light,ncomp=3,study=exp) #pls
# bad keepX.constraint
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576",100)),keepX=c(10,15),ncomp=3)
# bad ncomp
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576",100)),keepX=c(10,15),ncomp=2)


#with gene names in keepX.constraint
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)
which(res.spls.hybrid$loadings$X[,1]!=0)
which(res.spls.hybrid$loadings$X[,2]!=0)

#with numbers in keepX.constraint
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),keepX=c(100,50),ncomp=3)
which(res.spls.hybrid$loadings$X[,1]!=0)
which(res.spls.hybrid$loadings$X[,2]!=0)




system.time(wrapper.meta.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE,tol=1e-25))
system.time(mixOmics(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE,tol=1e-25))


head(res$loadings.global$Y)
head(res.sgcca$loadings$Y)
head(res.sgcca$loadings$X)
res$iter
res.sgcca$iter


quartz()
plot(res.sgcca$variates$X[,1],res.sgcca$variates$X[,2],col=as.numeric(factor(type.id)),pch=as.numeric(exp))


head(res$loadings.global$Y)
head(res.sgcca$loadings$Y)

cov(res$variates.global$X[,1],res$variates.global$Y[,1])
cov(res.sgcca$variates$X[,1],res.sgcca$variates$Y[,1])

cov(res$variates.global$X[,2],res$variates.global$Y[,2])
cov(res.sgcca$variates$X[,2],res.sgcca$variates$Y[,2])


cov(res$variates.global$X[,3],res$variates.global$Y[,3])
cov(res.sgcca$variates$X[,3],res.sgcca$variates$Y[,3])



#with only two classes, same results
if(FALSE)
{
    source("mixOmics/R/helpers.meta.R")
    source("mixOmics/R/meta.pls.R")
    source("mixOmics/R/meta.plsda.R")
    source("mixOmics/R/meta.spls.R")
    source("mixOmics/R/meta.splsda.R")
    source("mixOmics/R/meta.spls.hybrid.R")
    source("mixOmics/R/predict.meta.R")
    source("mixOmics/R/plotIndiv.meta.R")
    
    type.id.light=type.id[-which(type.id=="hESC")]
    data.light=data[-which(type.id=="hESC"),]
    study.light=exp[-which(type.id=="hESC")]
    Y.mat=unmap(type.id.light)
    rownames(Y.mat)=rownames(data.light)
    res=meta.spls.hybrid(X=data.light,Y=Y.mat,study=study.light,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,5),ncomp=3,near.zero.var=FALSE,tol=1e-25)
    
    res=meta.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,tol=1e-25)
    class(res)="meta.spls"
    plotIndiv(res,col=as.numeric(factor(type.id)),pch=as.numeric(exp),ind.names=FALSE)
    
    #meta.sgcca
    source("sgcca/meta.spls.hybrid.R")
    source("sgcca/sgccdaFunctions.R")
    source("sgcca/helpers.R")
    source("sgcca/pls.sgcca.R")
    source("sgcca/spls.sgcca.R")
    res.sgcca=meta.spls.hybrid.sgcca(X=data.light,Y=Y.mat,study=study.light,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,5),ncomp=3,near.zero.var=FALSE,tol=1e-25)
    
    X=data.light
    Y=Y.mat
    study=study.light
    keepX.constraint=list(c("ENSG00000006576","ENSG00000008226"))
    keepX=c(10,5)
    ncomp=3
    near.zero.var=FALSE
    scale=FALSE
    max.iter = 500
    tol = 1e-06
    
}


X=data
Y=Y.mat
study=exp
keepX.constraint=list(c("ENSG00000006576","ENSG00000008226"),c("ENSG00000001461","ENSG00000000457"),c("ENSG00000008226"))
keepX=c(10)
ncomp=4
near.zero.var=FALSE
scale=FALSE
max.iter = 500
tol = 1e-06


#-- validation des arguments --#
if (length(dim(X)) != 2)
stop("'X' must be a numeric matrix.")

X = as.matrix(X)
Y = as.matrix(Y)

if (!is.numeric(X) || !is.numeric(Y))
stop("'X' and/or 'Y' must be a numeric matrix.")

keepY=rep(ncol(Y),ncomp)
keepY.constraint=list()



if((length(keepX.constraint)+length(keepX))!=ncomp) stop("length (keepX.constraint) + length(keepX) should be ncomp")
if((length(keepY.constraint)+length(keepY))!=ncomp) stop("length (keepY.constraint) + length(keepY) should be ncomp")

keepX.temp=c(unlist(lapply(keepX.constraint,length)),keepX) #of length ncomp
keepY.temp=c(unlist(lapply(keepY.constraint,length)),keepY) #of length ncomp


check=Check.entry.pls(X,Y,ncomp,keepX.temp,keepY.temp)
X=check$X
Y=check$Y
ncomp=check$ncomp
X.names=check$X.names
Y.names=check$Y.names
ind.names=check$ind.names

if(length(study)!=nrow(X)) stop("unequal number of observations in 'X' and 'study'")
if(length(unique(rownames(X)))!=nrow(X)) stop("samples should have a unique identifier/rowname")

study=factor(study)

design = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)

# match keepX.constraint and the colnames of X in order for keepX.constraint to be a list of character
# safety if keepX.constraint contains a mixed of character/numeric. It should one or the other, not a mix
if(length(keepX.constraint)>0)
{
    if(!is.numeric(unlist(keepX.constraint)))
    {
        ind=match(unlist(keepX.constraint),colnames(X))
        if(sum(is.na(ind))>0) stop("'keepX.constraint' must contains a subset of colnames(X) or the position of the X-variables you wish to keep.")
    }
    X.indice=X[,unlist(keepX.constraint),drop=FALSE]
    keepX.constraint=relist(colnames(X.indice),skeleton=keepX.constraint)
}

# same for keepY.constraint
if(length(keepY.constraint)>0)
{
    if(!is.numeric(unlist(keepY.constraint)))
    {
        ind=match(unlist(keepY.constraint),colnames(Y))
        if(sum(is.na(ind))>0) stop("'keepY.constraint' must contains a subset of colnames(Y) or the position of the Y-variables you wish to keep.")
    }
    Y.indice=Y[,unlist(keepY.constraint),drop=FALSE]
    keepY.constraint=relist(colnames(Y.indice),skeleton=keepY.constraint)
}

# we need numbers in keepX.constraint from now on
#keepX.constraint= lapply(keepX.constraint,function(x){match(x,colnames(X))})
#keepY.constraint= lapply(keepY.constraint,function(x){match(x,colnames(Y))})

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
# keepA.constraint: keepX.constraint, which variables are kept on the first num.comp-1 components. It is a list
# near.zero.var: do you want to remove variables with very small variance

result <- sgcca(A = list(X = X, Y = Y), indY = 2, mode = "regression", ncomp = c(ncomp, ncomp), tol = tol, max.iter = max.iter,
design = design, keepA = list(keepX,keepY),keepA.constraint = list(keepX.constraint,keepY.constraint), sparse = TRUE,
lambda = NULL, scale = scale,
scheme = "centroid", study = study,near.zero.var=near.zero.var)



A = list(X = X, Y = Y)
indY = 2
mode = "regression"
ncomp = c(ncomp, ncomp)
tol = tol
max.iter = max.iter
design = design
keepA = list(keepX,keepY)
keepA.constraint = list(keepX.constraint,keepY.constraint)
sparse = TRUE
lambda = NULL#rep(1, length(list(X = X, Y = Y)))
scale = scale
scheme = "centroid"
study = study
verbose = FALSE
init = "svd"
bias = FALSE
tol=1e-25
near.zero.var=FALSE



check=Check.entry.sgcca(A, indY, design , lambda,ncomp , scheme , scale ,  bias,
init , tol , verbose,mode, sparse , max.iter,study , keepA, keepA.constraint,near.zero.var)
A=check$A
ncomp=check$ncomp
study=check$study

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
        }
    }
}


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
#for (n in 1:N)
#    {
 n=1
#if (verbose)
#        cat(paste0("Computation of the SGCCA block components #", n, " is under progress... \n"))
#        sgcca.result = NULL
#
#        ### Start: Estimation ai
#        if (sparse == TRUE)
#        {
            sgcca.result <- sgccak(R, design, lambda[n, ],study = study, keepA.constraint=if (!is.null(keepA.constraint)) {lapply(keepA.constraint, function(x){unlist(x[n])})} else {NULL} ,
           keepA = if (!is.null(keepA)) {lapply(keepA, function(x){x[n]})} else {NULL},indY = indY,
            scheme = scheme, init = init, max.iter = max.iter, tol = tol,   verbose = verbose)
#        } else {
#            sgcca.result <- rgccak(R, design, lambda[n, ], scheme = scheme, init = init, tol = tol, verbose = verbose, max.iter = max.iter)
#        }

# ==============================    SGCCAK

A=R
lambda=lambda[n, ]
keepA.constraint=if (!is.null(keepA.constraint)) {lapply(keepA.constraint, function(x){unlist(x[n])})} else {NULL}
keepA = if (!is.null(keepA)) {lapply(keepA, function(x){x[n]})} else {NULL}

    
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
    variates.A[, q] <- apply(A[[q]], 1, crossprod, loadings.A[[q]])
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
        time1=proc.time()
        Z_split=study_split(Z[,q,drop=FALSE],study)  # split Z by the study factor
        time2=proc.time()
        temp=Z[,q,drop=FALSE]
        rownames(temp)=rownames(A[[q]])
        Z_split=study_split_names(temp,study)
        time3=proc.time()
        ### Step A end: Compute the inner components
        
        
        ### Step B start: Computer the outer weight ###
        # possibility of removing NA (replacing by 0) and use crossprod, further development
        temp=0
        for (m in 1:nlevels_study)
        {
            loadings.partial.A.comp[[m]] <- t(A_split[[q]][[m]])%*%Z_split[[m]]#apply(t(A_split[[q]][[m]]), 1, crossprod, Z_split[[m]])
            temp=temp+loadings.partial.A.comp[[m]]
        }
        loadings.A[[q]]=temp
        
        # check the following
        #loadings.partial.A.comp = lapply(1 : J, function(x){lapply(1 : nlevels_study, function(y){ A_split[[x]][[y]] %*% Z_split[[x]]})})
        
        
        loadings.A[[q]]=sparsity(loadings.A[[q]],keepA[[q]],keepA.constraint[[q]],lambda)
        loadings.A[[q]]=l2.norm(as.vector(loadings.A[[q]]))
        
        ### Step B end: Computer the outer weight ###
        variates.A[, q] <- A[[q]]%*%loadings.A[[q]]#apply(A[[q]], 1, crossprod, loadings.A[[q]])
    }
  
    print(iter)
    print(head(loadings.A_old[[1]]))
    print(head(loadings.A[[1]]))
    print(head(loadings.A[[2]]))
  
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




#----------------------------
# check convergence on covariance
#compute sum of covariance (criterion to maximize)






head(loadings.A[[1]])
head(loadings.A[[2]])
head(variates.A)

loadings.A_old.2[[1]][which(loadings.A_old.2[[1]]!=0)]
 loadings.A_old[[1]][which(loadings.A_old[[1]]!=0)]
loadings.A[[1]][which(loadings.A[[1]]!=0)]































        
        
        
        
        
        
        
        
        
        
        ### End: Estimation ai
        
        loadings.partial.A[[n]]=sgcca.result$loadings.partial.A.comp
        variates.partial.A[[n]]=sgcca.result$variates.partial.A.comp
        
        AVE_inner[n] <- sgcca.result$AVE_inner
        crit[[n]] <- sgcca.result$crit
        for (k in 1:J) variates.A[[k]][, n] <- sgcca.result$variates.A[, k]
        for (q in which(n < ndefl)) if (sum(sgcca.result$loadings.A[[q]] != 0) <= 1)
        warning(sprintf("Deflation failed because only one variable was selected for block #", q, "! \n"))
        
        
        ### Start: estimations of vector b, c and t
        # might be remove in the near future as the prediction fonction won't use mat.c mat.b etc
        t_temp <- c_temp <- b_temp <- NULL
        t_temp = lapply(1:J, function(x) {apply(R[[x]], 1, crossprod, sgcca.result$loadings.A[[x]])})
        
        if (is.null(indY)) {
            c_temp = lapply(1:J, function(x){as.matrix(apply(t(R[[x]]), 1, crossprod, t_temp[[x]])/(norm2(t_temp[[x]]))^2)})
        } else if (mode == "classic" || mode == "invariant"){
            c_temp = lapply(c(1:J), function(x){as.matrix(apply(t(R[[x]]), 1, crossprod, t_temp[[x]])/(norm2(t_temp[[x]]))^2)})
        } else if (mode == "canonical"){
            c_temp = lapply(c(1:J), function(x){as.matrix(apply(t(R[[x]]), 1, crossprod, t_temp[[x]])/(norm2(t_temp[[x]]))^2)})
            b_temp = lapply(indY, function(x){as.matrix(apply(t(R[[x]]), 1, crossprod, t_temp[[x]])/(norm2(t_temp[[x]]))^2)})
        } else if (mode == "regression"){
            c_temp = lapply(c(1:J), function(x){as.matrix(apply(t(R[[x]]), 1, crossprod, t_temp[[x]])/(norm2(t_temp[[x]]))^2)})
            b_temp = lapply(c(1:J)[-indY], function(x){as.matrix(apply(t(R[[indY]]), 1, crossprod, t_temp[[x]])/(norm2(t_temp[[x]]))^2)})
        }
        
        if (is.null(indY)){
            for (k in 1:J) {
                if (!is.null(c_temp)) c[[k]][, n] <- c_temp[[k]]
                if (!is.null(t_temp)) t[[k]][, n] <- t_temp[[k]]
            }
        } else {
            for (k in 1:J) {
                if (!is.null(b_temp) & k == indY) b[[k]][, n] <- b_temp[[1]]
                if (!is.null(c_temp) & k != indY) c[[k]][, n] <- c_temp[[k]]
                if (!is.null(t_temp)) t[[k]][, n] <- t_temp[[k]]
            }
        }
        ### End: Estimations of vectors b, c and t
        
        # deflation if there are more than 1 component and if we haven't reach the max number of component(N)
        if (N != 1 & n != N) {
            defla.result <- defl.select(sgcca.result$variates.A, R, ndefl, n, nbloc = J, indY = indY, mode = mode, aa = sgcca.result$loadings.A)
            R <- defla.result$resdefl
            defl.matrix[[n + 1]] <- R
        }
        
        for (k in 1 : J) {
            if (N != 1) {
                P[[k]][, n - 1] <- defla.result$pdefl[[k]]
            }
            loadings.A[[k]][, n] <- sgcca.result$loadings.A[[k]]
        }
        
        if (n == 1) {
            for (k in 1 : J) loadings.Astar[[k]][, n] <- sgcca.result$loadings.A[[k]]
        } else {
            for (k in 1 : J) loadings.Astar[[k]][, n] <- sgcca.result$loadings.A[[k]] - loadings.Astar[[k]][, (1 : n - 1), drop = F] %*% drop(t(loadings.A[[k]][, n]) %*% P[[k]][, 1 : (n - 1), drop = F])
        }
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









checklist = names(classic.linn.pls)
checklist = checklist[!checklist %in% c("iter")]

classic.linn.pls.abs = abslist(classic.linn.pls, c("variates", "loadings", "mat.c", "mat.d"))
classic.linn.pls.sgcca.abs = abslist(classic.linn.pls.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(classic.linn.pls.abs[names(classic.linn.pls.abs) == x],
                                        classic.linn.pls.sgcca.abs[names(classic.linn.pls.sgcca.abs) == x])})

head(classic.linn.pls.abs$loadings$X); head(classic.linn.pls.sgcca.abs$loadings$X);

indiv1 <- c(200, 40, 60); indiv2 <- c(190, 45, 45)
newdata <- rbind(indiv1, indiv2)
colnames(newdata) <- colnames(X)
newdata

pred.classic.linn.pls <- predict(classic.linn.pls, newdata)
pred.classic.linn.pls.sgcca <- unlist(predict.sgcca(classic.linn.pls.sgcca, newdata), recursive = FALSE)

names(pred.classic.linn.pls)

max(abs(pred.classic.linn.pls$predict) - abs(pred.classic.linn.pls.sgcca$predict))
max(abs(pred.classic.linn.pls$variates) - abs(pred.classic.linn.pls.sgcca$variates))
max(abs(pred.classic.linn.pls$B.hat) - abs(pred.classic.linn.pls.sgcca$B.hat))

system.time(perf(classic.linn.pls, method.predict = "all", validation = "loo", progressBar = FALSE))
system.time(perf.sgcca(classic.linn.pls.sgcca, validation = "loo"))

classic.linn.pls.perf = perf(classic.linn.pls, method.predict = "all", validation = "loo")
classic.linn.pls.perf.sgcca = perf.sgcca(classic.linn.pls.sgcca, validation = "loo")

names(classic.linn.pls.perf)
all.equal(classic.linn.pls.perf$MSEP, classic.linn.pls.perf.sgcca$MSEP$X.1)
all.equal(classic.linn.pls.perf$R2, classic.linn.pls.perf.sgcca$R2$X.1)
all.equal(classic.linn.pls.perf$Q2, classic.linn.pls.perf.sgcca$Q2$X.1) # Difference
all.equal(classic.linn.pls.perf$Q2.total, classic.linn.pls.perf.sgcca$Q2.total$X.1) # Difference
max(abs(classic.linn.pls.perf$press.mat[, , 1] - classic.linn.pls.perf.sgcca$press.mat$X.1[[1]]))
max(abs(classic.linn.pls.perf$press.mat[, , 2] - classic.linn.pls.perf.sgcca$press.mat$X.1[[2]]))
max(abs(classic.linn.pls.perf$RSS.indiv[, , 2] - classic.linn.pls.perf.sgcca$press.mat$X.1[[1]]))
max(abs(classic.linn.pls.perf$RSS.indiv[, , 3] - classic.linn.pls.perf.sgcca$press.mat$X.1[[2]]))
max(abs(classic.linn.pls.perf$PRESS.inside - classic.linn.pls.perf.sgcca$PRESS.inside$X.1))
max(abs(classic.linn.pls.perf$RSS - classic.linn.pls.perf.sgcca$RSS$X.1)) # Difference

classic.linn.pls.perf$Q2

classic.linn.pls.perf.sgcca$Q2

### Invariant
invariant.linn.pls <- pls(X, Y, mode = "invariant", tol = 1e-25)
invariant.linn.pls.sgcca <- pls.sgcca(X, Y, mode = "invariant", tol = 1e-25)

invariant.linn.pls.abs = abslist(invariant.linn.pls, c("variates", "loadings", "mat.c"))
invariant.linn.pls.sgcca.abs = abslist(invariant.linn.pls.sgcca, c("variates", "loadings", "mat.c"))

lapply(checklist, function(x){all.equal(invariant.linn.pls.abs[names(invariant.linn.pls.abs) == x], invariant.linn.pls.sgcca.abs[names(invariant.linn.pls.sgcca.abs) == x])})

head(invariant.linn.pls$loadings$X); head(invariant.linn.pls.sgcca$loadings$X);

pred.invariant.linn.pls <- predict(invariant.linn.pls, newdata)
pred.invariant.linn.pls.sgcca <- unlist(predict.sgcca(invariant.linn.pls.sgcca, newdata), recursive = FALSE)

names(pred.invariant.linn.pls)

max(abs(pred.invariant.linn.pls$predict) - abs(pred.invariant.linn.pls.sgcca$predict))
max(abs(pred.invariant.linn.pls$variates) - abs(pred.invariant.linn.pls.sgcca$variates))
max(abs(pred.invariant.linn.pls$B.hat) - abs(pred.invariant.linn.pls.sgcca$B.hat))

### Canonical
canonical.linn.pls <- pls(X, Y, mode = "canonical", tol = 1e-25)
canonical.linn.pls.sgcca <- pls.sgcca(X, Y, mode = "canonical", tol = 1e-25)

canonical.linn.pls.abs = abslist(canonical.linn.pls, c("variates", "loadings", "mat.c", "mat.e"))
canonical.linn.pls.sgcca.abs = abslist(canonical.linn.pls.sgcca, c("variates", "loadings", "mat.c", "mat.e"))

lapply(checklist, function(x){all.equal(canonical.linn.pls.abs[names(canonical.linn.pls.abs) == x], canonical.linn.pls.sgcca.abs[names(canonical.linn.pls.sgcca.abs) == x])})

head(canonical.linn.pls$loadings$X); head(canonical.linn.pls.sgcca$loadings$X);

pred.canonical.linn.pls <- predict(canonical.linn.pls, newdata)
pred.canonical.linn.pls.sgcca <- unlist(predict.sgcca(canonical.linn.pls.sgcca, newdata), recursive = FALSE)

names(pred.canonical.linn.pls)

max(abs(pred.canonical.linn.pls$predict) - abs(pred.canonical.linn.pls.sgcca$predict))
max(abs(pred.canonical.linn.pls$variates) - abs(pred.canonical.linn.pls.sgcca$variates))
max(abs(pred.canonical.linn.pls$B.hat) - abs(pred.canonical.linn.pls.sgcca$B.hat))

### Regression
regression.linn.pls <- pls(X, Y, mode = "regression", tol = 1e-25, ncomp = 3)
regression.linn.pls.sgcca <- pls.sgcca(X, Y, mode = "regression", tol = 1e-25, ncomp = 3)

regression.linn.pls.abs = abslist(regression.linn.pls, c("variates", "loadings", "mat.c", "mat.d"))
regression.linn.pls.sgcca.abs = abslist(regression.linn.pls.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(regression.linn.pls.abs[names(regression.linn.pls.abs) == x],
                                        regression.linn.pls.sgcca.abs[names(regression.linn.pls.sgcca.abs) == x])})

pred.regression.linn.pls <- predict(regression.linn.pls, newdata)
pred.regression.linn.pls.sgcca <- unlist(predict.sgcca(regression.linn.pls.sgcca, newdata), recursive = FALSE)

names(pred.regression.linn.pls)

max(abs(pred.regression.linn.pls$predict) - abs(pred.regression.linn.pls.sgcca$predict))
max(abs(pred.regression.linn.pls$variates) - abs(pred.regression.linn.pls.sgcca$variates))
max(abs(pred.regression.linn.pls$B.hat) - abs(pred.regression.linn.pls.sgcca$B.hat))

system.time(perf(regression.linn.pls, method.predict = "all", validation = "loo", progressBar = FALSE))
system.time(perf.sgcca(regression.linn.pls.sgcca, validation = "loo"))

### Example 2 ###
rm(list=ls())
source("pls.sgcca.R")

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
row.names(Y) <- paste("X", row.names(Y), sep = "")

### Classic
classic.tox.pls <- pls(X, Y, ncomp = 3, mode = "classic", tol = 1e-25)
classic.tox.pls.sgcca <- pls.sgcca(X, Y, ncomp = 3, mode = "classic", tol = 1e-25)

checklist = names(classic.tox.pls)
checklist = checklist[!checklist %in% c("iter")]

classic.tox.pls = abslist(classic.tox.pls, c("variates", "loadings", "mat.c"))
classic.tox.pls.sgcca = abslist(classic.tox.pls.sgcca, c("variates", "loadings", "mat.c"))

lapply(checklist, function(x){all.equal(classic.tox.pls[names(classic.tox.pls) == x], classic.tox.pls.sgcca[names(classic.tox.pls.sgcca) == x])})

### Invariant
invariant.tox.pls <- pls(X, Y, ncomp = 3, mode = "invariant", tol = 1e-25)
invariant.tox.pls.sgcca <- pls.sgcca(X, Y, ncomp = 3, mode = "invariant", tol = 1e-25)

invariant.tox.pls = abslist(invariant.tox.pls, c("variates", "loadings", "mat.c"))
invariant.tox.pls.sgcca = abslist(invariant.tox.pls.sgcca, c("variates", "loadings", "mat.c"))

lapply(checklist, function(x){all.equal(invariant.tox.pls[names(invariant.tox.pls) == x], invariant.tox.pls.sgcca[names(invariant.tox.pls.sgcca) == x])})

### Canonical
canonical.tox.pls <- pls(X, Y, ncomp = 3, mode = "canonical", tol = 1e-25)
canonical.tox.pls.sgcca <- pls.sgcca(X, Y, ncomp = 3, mode = "canonical", tol = 1e-25)

canonical.tox.pls = abslist(canonical.tox.pls, c("variates", "loadings", "mat.c", "mat.e"))
canonical.tox.pls.sgcca = abslist(canonical.tox.pls.sgcca, c("variates", "loadings", "mat.c", "mat.e"))

lapply(checklist, function(x){all.equal(canonical.tox.pls[names(canonical.tox.pls) == x], canonical.tox.pls.sgcca[names(canonical.tox.pls.sgcca) == x])})

### Regression
regression.tox.pls <- pls(X, Y, ncomp = 3, mode = "regression", tol = 1e-25)
regression.tox.pls.sgcca <- pls.sgcca(X, Y, ncomp = 3, mode = "regression", tol = 1e-25)

regression.tox.pls = abslist(regression.tox.pls, c("variates", "loadings", "mat.c", "mat.d"))
regression.tox.pls.sgcca = abslist(regression.tox.pls.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(regression.tox.pls[names(regression.tox.pls) == x], regression.tox.pls.sgcca[names(regression.tox.pls.sgcca) == x])})


### Example 3 ### Book M. Tenenhaus "La regression PLS - Theorie et pratique" (page 78)
rm(list=ls())
source("pls.sgcca.R"); source("predict.sgcca.R"); source("predict.sgcca.R");
source("perf.sgcca.R"); source("mix.options.R")

X <- matrix(c(0.00, 0.23, 0.00, 0.00, 0.00, 0.74, 0.03,
              0.00, 0.10, 0.00, 0.00, 0.12, 0.74, 0.04,
              0.00, 0.00, 0.00, 0.10, 0.12, 0.74, 0.04,
              0.00, 0.49, 0.00, 0.00, 0.12, 0.37, 0.02,
              0.00, 0.00, 0.00, 0.62, 0.12, 0.18, 0.08,
              0.00, 0.62, 0.00, 0.00, 0.00, 0.37, 0.01,
              0.17, 0.27, 0.10, 0.38, 0.00, 0.00, 0.08,
              0.17, 0.19, 0.10, 0.38, 0.02, 0.06, 0.08,
              0.17, 0.21, 0.10, 0.38, 0.00, 0.06, 0.08,
              0.17, 0.15, 0.10, 0.38, 0.02, 0.10, 0.08, 
              0.21, 0.36, 0.12, 0.25, 0.00, 0.00, 0.06,
              0.00, 0.00, 0.00, 0.55, 0.00, 0.37, 0.08), nrow = 12, ncol = 7, byrow = TRUE,
            dimnames = list(paste("Ind", 1:12), paste("Col", 1 : 7)))

Y <- c(98.7, 97.8, 96.6, 92.0, 86.6, 91.2, 81.9, 83.1, 82.4, 83.2, 81.4, 88.1)

example3.pls.reg <- spls(X, Y, mode = "regression", ncomp = 4)
example3.pls.reg.perf = perf(example3.pls.reg, validation = "loo")
example3.pls.reg.perf$Q2
example3.pls.reg.perf$Q2.total

example3.pls.sgcca.reg <- pls.sgcca(X, Y, mode = "regression", tol = 1e-25, ncomp = 4)
example3.pred.sgcca.reg <- predict.sgcca(example3.pls.sgcca.reg, X)
example3.pls.sgcca.reg.perf = perf.sgcca(example3.pls.sgcca.reg, validation = "loo")
example3.pls.sgcca.reg.perf$Q2
example3.pls.sgcca.reg.perf$Q2.total

example3.pls.sgcca.reg$loadings
example3.pls.sgcca.reg$variates

################################
############# SPLS #############
################################
rm(list=ls())
source("spls.sgcca.R"); source("helpers.R")
library(mixOmics);

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
row.names(Y) <- paste("X", row.names(Y), sep = "")

### Regression
regression.tox.spls <- spls(X, Y, ncomp = 2, mode = "regression", keepX = c(50, 50), keepY = c(10, 10))
regression.tox.spls.sgcca <- spls.sgcca(X, Y, ncomp = 2, mode = "regression", keepX = c(50, 50), keepY = c(10, 10))

checklist = names(regression.tox.spls)
checklist = checklist[!checklist %in% c("iter")]

regression.tox.spls = abslist(regression.tox.spls, c("variates", "loadings", "mat.c", "mat.d"))
regression.tox.spls.sgcca = abslist(regression.tox.spls.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(regression.tox.spls[names(regression.tox.spls) == x], regression.tox.spls.sgcca[names(regression.tox.spls.sgcca) == x])})

### Canonical
canonical.tox.spls <- spls(X, Y, ncomp = 2, mode = "canonical", keepX = c(50, 50), keepY = c(10, 10))
canonical.tox.spls.sgcca <- spls.sgcca(X, Y, ncomp = 2, mode = "canonical", keepX = c(50, 50), keepY = c(10, 10))

canonical.tox.spls = abslist(canonical.tox.spls, c("variates", "loadings", "mat.c", "mat.e"))
canonical.tox.spls.sgcca = abslist(canonical.tox.spls.sgcca, c("variates", "loadings", "mat.c", "mat.e"))

lapply(checklist, function(x){all.equal(canonical.tox.spls[names(canonical.tox.spls) == x], canonical.tox.spls.sgcca[names(canonical.tox.spls.sgcca) == x])})

##################################
############# PLS-DA #############
##################################
rm(list=ls())
source("splsda.sgcca.R"); source("sgccda.sgcca.R"); source("predict.sgcca.R");
source("perf.sgcca.R"); source("mix.options.R")
library(mixOmics);

### Example 1 ###
data(breast.tumors)
X <- breast.tumors$gene.exp
# X <- t(na.omit(t(X)))  ### discrepancy with missing data
Y <- breast.tumors$sample$treatment

breast.plsda <- plsda(X, Y, ncomp = 1)
breast.plsda.sgcca <- plsda.sgcca(X, Y, ncomp = 1)

checklist = names(breast.plsda)
checklist = checklist[!checklist %in% c("iter")]

breast.plsda.abs = abslist(breast.plsda, c("variates", "loadings", "mat.c", "mat.d"))
breast.plsda.sgcca.abs = abslist(breast.plsda.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(breast.plsda.abs[names(breast.plsda.abs) == x], 
                                        breast.plsda.sgcca.abs[names(breast.plsda.sgcca.abs) == x])})

breast.plsda.predict = predict(breast.plsda, X)
breast.plsda.predict.sgcca <- predict.sgcca(breast.plsda.sgcca, newdata = X, method = "all")

names(breast.plsda.predict); names(breast.plsda.predict.sgcca)
max(abs(breast.plsda.predict$predict) - abs(breast.plsda.predict.sgcca$predict[[1]]))
max(abs(breast.plsda.predict$variates) - abs(breast.plsda.predict.sgcca$variates[[1]]))
max(abs(breast.plsda.predict$B.hat) - abs(breast.plsda.predict.sgcca$B.hat[[1]]))
max(abs(breast.plsda.predict$centroids) - abs(breast.plsda.predict.sgcca$centroids[[1]]))
breast.plsda.predict$method == breast.plsda.predict.sgcca$method
max(abs(breast.plsda.predict$class[[1]]) - abs(unlist(breast.plsda.predict.sgcca$class, recursive = FALSE)[[1]]))

breast.plsda.perf = perf(breast.plsda, validation = "loo")
breast.plsda.perf.sgcca = perf.sgcca(breast.plsda.sgcca, validation = "loo")
max(breast.plsda.perf$error.rate - breast.plsda.perf.sgcca$error.rate.X$X)

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "breast.plsda.predict.sgcca", "map", "plsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "crossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

for (i in 1 : 10){
  folds = breast.plsda.perf = breast.plsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); breast.plsda.perf = perf(breast.plsda, validation = "Mfold", folds = folds)
  set.seed(i); breast.plsda.perf.sgcca = perf.sgcca(breast.plsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(breast.plsda.perf$error.rate - breast.plsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(breast.plsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(breast.plsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

### Example 2 ###
rm(list = ls())
source("plsda.sgcca.R"); source("sgccda.sgcca.R"); source("predict.sgcca.R")
source("perf.sgcca.R"); source("mix.options.R")
require(mixOmics)

data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.factor(liver.toxicity$treatment[, 4])

#### Difference because not the same starting point in the algo between function splsda and plsda in mixOmics package (indeed between pls and spls, see above)
#### To play down this discrepancy, set up tol = 1e-25 instead of 1e-06 
liver.plsda <- plsda(X, Y, ncomp = 2, tol = 1e-25)               
liver.plsda.sgcca <- plsda.sgcca(X, Y, ncomp = 2, tol = 1e-25)

checklist = names(liver.plsda)
checklist = checklist[!checklist %in% c("iter")]

liver.plsda.abs = abslist(liver.plsda, c("variates", "loadings", "mat.d", "mat.c"))
liver.plsda.sgcca.abs = abslist(liver.plsda.sgcca, c("variates", "loadings", "mat.d", "mat.c"))

lapply(checklist, function(x){all.equal(liver.plsda.abs[names(liver.plsda.abs) == x], 
                                        liver.plsda.sgcca.abs[names(liver.plsda.sgcca.abs) == x])})

# Comparison prediction
samp <- sample(1:5, nrow(X), replace = TRUE)  
test <- which(samp == 1)   # testing on the first fold
train <- setdiff(1:nrow(X), test)

plsda.train <- plsda(X[train, ], Y[train], ncomp = 2, tol = 1e-25)
test.predict <- predict(plsda.train, X[test, ], method = "max.dist")

plsda.train.sgcca <- plsda.sgcca(X[train, ], Y[train], ncomp = 2, tol = 1e-25)
test.predict.sgcca <- predict.sgcca(plsda.train.sgcca, X[test, ], method = "max.dist")

names(test.predict); names(test.predict.sgcca)
max(abs(test.predict$predict) - abs(test.predict.sgcca$predict[[1]]))
max(abs(test.predict$variates) - abs(test.predict.sgcca$variates[[1]]))
max(abs(test.predict$B.hat) - abs(test.predict.sgcca$B.hat[[1]]))
max(abs(test.predict$centroids) - abs(test.predict.sgcca$centroids[[1]]))
test.predict$method == test.predict.sgcca$method
max(abs(test.predict$class[[1]]) - abs(unlist(test.predict.sgcca$class, recursive = FALSE)[[1]]))

Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
cbind(Y = as.character(Y[test]), Prediction)

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "liver.plsda.sgcca", "map", "plsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "crossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

liver.plsda.perf = perf(liver.plsda, validation = "loo")
liver.plsda.perf.sgcca = perf.sgcca(liver.plsda.sgcca, validation = "loo")
max(liver.plsda.perf$error.rate - liver.plsda.perf.sgcca$error.rate.X$X)

for (i in 1 : 5){
  folds = liver.plsda.perf = liver.plsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); liver.plsda.perf = perf(liver.plsda, validation = "Mfold", folds = folds)
  set.seed(i); liver.plsda.perf.sgcca = perf.sgcca(liver.plsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(liver.plsda.perf$error.rate - liver.plsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(liver.plsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(liver.plsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

###################################
############# SPLS-DA #############
###################################
rm(list=ls())
source("splsda.sgcca.R"); source("sgccda.sgcca.R"); source("predict.sgcca.R")
source("perf.sgcca.R"); source("mix.options.R")
library(mixOmics);

### Example 1 ###
data(breast.tumors)
X <- breast.tumors$gene.exp
#X <- t(na.omit(t(X)))  ### issue with missing data
Y <- breast.tumors$sample$treatment

breast.splsda <- splsda(X, Y, ncomp = 3)#, keepX = c(25, 25))
breast.splsda.sgcca <- splsda.sgcca(X, Y, ncomp =2,  keepX = c(20, 20))
breast.plsda.sgcca <- plsda.sgcca(X, Y, ncomp =2)

checklist = names(breast.splsda)
checklist = checklist[!checklist %in% c("iter")]

breast.splsda.abs = abslist(breast.splsda, c("variates", "loadings", "mat.c", "mat.d"))
breast.splsda.sgcca.abs = abslist(breast.splsda.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(breast.splsda.abs[names(breast.splsda.abs) == x], 
                                        breast.splsda.sgcca.abs[names(breast.splsda.sgcca.abs) == x])})

# Comparison prediction
breast.splsda.predict = predict(breast.splsda, X)
breast.splsda.predict.sgcca <- predict.sgcca(breast.splsda.sgcca, newdata = X, method = "all")

names(breast.splsda.predict); names(breast.splsda.predict.sgcca)
max(abs(breast.splsda.predict$predict) - abs(breast.splsda.predict.sgcca$predict[[1]]))
max(abs(breast.splsda.predict$variates) - abs(breast.splsda.predict.sgcca$variates[[1]]))
max(abs(breast.splsda.predict$B.hat) - abs(breast.splsda.predict.sgcca$B.hat[[1]]))
max(abs(breast.splsda.predict$centroids) - abs(breast.splsda.predict.sgcca$centroids[[1]]))
breast.splsda.predict$method == breast.splsda.predict.sgcca$method
max(abs(breast.splsda.predict$class[[1]]) - abs(unlist(breast.splsda.predict.sgcca$class, recursive = FALSE)[[1]]))

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "breast.splsda.predict.sgcca", "map", "splsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "crossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

breast.splsda.perf = perf(breast.splsda, validation = "loo")
breast.splsda.perf.sgcca = perf.sgcca(breast.splsda.sgcca, validation = "loo")
max(breast.splsda.perf$error.rate - breast.splsda.perf.sgcca$error.rate.X$X)

for (i in 1 : 5){
  folds = breast.splsda.perf = breast.splsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); breast.splsda.perf = perf(breast.splsda, validation = "Mfold", folds = folds)
  set.seed(i); breast.splsda.perf.sgcca = perf.sgcca(breast.splsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(breast.splsda.perf$error.rate - breast.splsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(breast.splsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(breast.splsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

### Example 2 ###
rm(list=ls())
source("splsda.sgcca.R"); source("sgccda.sgcca.R"); source("predict.sgcca.R")
source("perf.sgcca.R"); source("mix.options.R")
library(mixOmics);

data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.factor(liver.toxicity$treatment[, 4])

liver.splsda = splsda(X, Y, ncomp = 2, keepX = c(5, 5))
liver.splsda.sgcca = splsda.sgcca(X, Y, ncomp = 2, keepX = c(5, 5))

checklist = names(liver.splsda)
checklist = checklist[!checklist %in% c("iter")]

liver.splsda.abs = abslist(liver.splsda, c("variates", "loadings", "mat.d", "mat.c"))
liver.splsda.sgcca.abs = abslist(liver.splsda.sgcca, c("variates", "loadings", "mat.d", "mat.c"))

lapply(checklist, function(x){all.equal(liver.splsda.abs[names(liver.splsda.abs) == x],
                                        liver.splsda.sgcca.abs[names(liver.splsda.sgcca.abs) == x])})

# Comparison prediction
samp <- sample(1:5, nrow(X), replace = TRUE)  
test <- which(samp == 1)   # testing on the first fold
train <- setdiff(1:nrow(X), test)

splsda.train <- splsda(X[train, ], Y[train], ncomp = 2, keepX = c(30, 30))
test.predict <- predict(splsda.train, X[test, ], method = "max.dist")

splsda.train.sgcca <- splsda.sgcca(X[train, ], Y[train], ncomp = 2, keepX = c(30, 30))
test.predict.sgcca <- predict.sgcca(splsda.train.sgcca, X[test, ], method = "max.dist")

names(test.predict); names(test.predict.sgcca)
max(abs(test.predict$predict) - abs(test.predict.sgcca$predict[[1]]))
max(abs(test.predict$variates) - abs(test.predict.sgcca$variates[[1]]))
max(abs(test.predict$B.hat) - abs(test.predict.sgcca$B.hat[[1]]))
max(abs(test.predict$centroid) - abs(test.predict.sgcca$centroid[[1]]))
test.predict$method == test.predict.sgcca$method
max(abs(test.predict$class[[1]]) - abs(unlist(test.predict.sgcca$class, recursive = FALSE)[[1]]))

Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
cbind(Y = as.character(Y[test]), Prediction)

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "liver.splsda.sgcca", "map", "splsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "crossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

liver.splsda.perf = perf(liver.splsda, validation = "loo")
liver.splsda.perf.sgcca = perf.sgcca(liver.splsda.sgcca, validation = "loo")
max(liver.splsda.perf$error.rate - liver.splsda.perf.sgcca$error.rate.X$X)

for (i in 1 : 5){
  folds = liver.splsda.perf = liver.splsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); liver.splsda.perf = perf(liver.splsda, validation = "Mfold", folds = folds)
  set.seed(i); liver.splsda.perf.sgcca = perf.sgcca(liver.splsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(liver.splsda.perf$error.rate - liver.splsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(liver.splsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(liver.splsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

### Test perf.splsda ###
rm(list=ls())
source("splsda.sgcca.R"); source("perf.sgcca.R"); source("predict.sgcca.R")
source("perf.sgcca.R"); source("mix.options.R")
library(mixOmics);

data(srbct)
X <- srbct$gene
Y <- srbct$class  

ncomp = 5

srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))  
srbct.splsda.sgcca <- splsda.sgcca(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))  

# with Mfold
# ---------
set.seed(2)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8, method.predict = "all")
error$error.rate

set.seed(2)
error.sgcca <- perf.sgcca(srbct.splsda.sgcca, validation = "Mfold", folds = 8, method.predict = "all")
error.sgcca$error.rate.X

max(abs(error$error.rate - error.sgcca$error.rate.X$X))

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "srbct.splsda.sgcca", "map", "splsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "crossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

srbct.splsda.perf = perf(srbct.splsda, validation = "loo")
srbct.splsda.perf.sgcca = perf.sgcca(srbct.splsda.sgcca, validation = "loo")
max(srbct.splsda.perf$error.rate - srbct.splsda.perf.sgcca$error.rate.X$X)

all.equal(srbct.splsda.perf$features$stable, srbct.splsda.perf.sgcca$features$stable$X)
all.equal(srbct.splsda.perf$features$final, srbct.splsda.perf.sgcca$features$final$X)

for (i in 1 : 5){
  folds = srbct.splsda.perf = srbct.splsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); srbct.splsda.perf = perf(srbct.splsda, validation = "Mfold", folds = folds)
  set.seed(i); srbct.splsda.perf.sgcca = perf.sgcca(srbct.splsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(srbct.splsda.perf$error.rate - srbct.splsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(srbct.splsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(srbct.splsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

#########################################
############# wrapper.rgcca #############
#########################################
rm(list=ls())
source("/Users/bgautier/Desktop/rGCCA-DA/rGCCA codes/sgccdaFunctions_BG3.R")
library(mixOmics); #library(RGCCA)

data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid, Y)
# with this design, gene expression and lipids are connected to the diet factor
design = matrix(c(0,0,1,
                  0,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

### centroid ###
centroid.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0.5), ncomp = c(2, 2, 1), scheme = "centroid", verbose = FALSE, tol = 1e-25)
centroid.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0.5), verbose = FALSE, scheme = "centroid", tol = 1e-25) 

checklist = names(centroid.wrap.result.rgcca)
checklist = checklist[!checklist %in% c("class", "data", "names")]

lapply(checklist, function(x){all.equal(centroid.wrap.result.rgcca[names(centroid.wrap.result.rgcca) == x], centroid.sgcca.nutrimouse.rgccak[names(centroid.sgcca.nutrimouse.rgccak) == x])})

### horst ###
horst.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0.5), ncomp = c(2, 2, 1), scheme = "horst", verbose = FALSE, tol = 1e-25)
horst.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0.5), verbose = FALSE, scheme = "horst", tol = 1e-25) 

lapply(checklist, function(x){all.equal(horst.wrap.result.rgcca[names(horst.wrap.result.rgcca) == x], horst.sgcca.nutrimouse.rgccak[names(horst.sgcca.nutrimouse.rgccak) == x])})

### factorial ###
factorial.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0.5), ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE, tol = 1e-25)
factorial.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0.5), verbose = FALSE, scheme = "factorial", tol = 1e-25) 

lapply(checklist, function(x){all.equal(factorial.wrap.result.rgcca[names(factorial.wrap.result.rgcca) == x], factorial.sgcca.nutrimouse.rgccak[names(factorial.sgcca.nutrimouse.rgccak) == x])})

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#note: the tau parameter is the regularization parameter

### centroid ###
centroid.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0), ncomp = c(2, 2, 1), scheme = "centroid", verbose = FALSE, tol = 1e-25)
centroid.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0), verbose = FALSE, scheme = "centroid", tol = 1e-25) 

lapply(checklist, function(x){all.equal(centroid.wrap.result.rgcca[names(centroid.wrap.result.rgcca) == x], centroid.sgcca.nutrimouse.rgccak[names(centroid.sgcca.nutrimouse.rgccak) == x])})

### horst ###
horst.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0), ncomp = c(2, 2, 1), scheme = "horst", verbose = FALSE, tol = 1e-25)
horst.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0), verbose = FALSE, scheme = "horst", tol = 1e-25) 

lapply(checklist, function(x){all.equal(horst.wrap.result.rgcca[names(horst.wrap.result.rgcca) == x], horst.sgcca.nutrimouse.rgccak[names(horst.sgcca.nutrimouse.rgccak) == x])})

### factorial ###
factorial.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0), ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE, tol = 1e-25)
factorial.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0), verbose = FALSE, scheme = "factorial", tol = 1e-25) 

lapply(checklist, function(x){all.equal(factorial.wrap.result.rgcca[names(factorial.wrap.result.rgcca) == x], factorial.sgcca.nutrimouse.rgccak[names(factorial.sgcca.nutrimouse.rgccak) == x])})

#########################################
############# wrapper.sgcca #############
#########################################
rm(list=ls())
source("/Users/bgautier/Desktop/rGCCA-DA/rGCCA codes/sgccdaFunctions_BG.R")
library(mixOmics); #library(RGCCA)

data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
# with this design, gene expression and lipids are connected to the diet factor
design = matrix(c(0,0,1,
                  0,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

### centroid ###
centroid.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "centroid", verbose = FALSE, tol = 1e-25)
centroid.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(.3, .5, 1), verbose = FALSE, scheme = "centroid", tol = 1e-25, sparse = TRUE) 

checklist = names(centroid.wrap.result.sgcca)
checklist = checklist[!checklist %in% c("class", "data", "names")]

lapply(checklist, function(x){all.equal(centroid.wrap.result.sgcca[names(centroid.wrap.result.sgcca) == x], centroid.sgcca.nutrimouse.sgccak[names(centroid.sgcca.nutrimouse.sgccak) == x])})

### horst ###
horst.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "horst", verbose = FALSE, tol = 1e-25)
horst.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(.3, .5, 1), verbose = FALSE, scheme = "horst", tol = 1e-25, sparse = TRUE) 

lapply(checklist, function(x){all.equal(horst.wrap.result.sgcca[names(horst.wrap.result.sgcca) == x], horst.sgcca.nutrimouse.sgccak[names(horst.sgcca.nutrimouse.sgccak) == x])})

### factorial ###
factorial.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE, tol = 1e-25)
factorial.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(.3, .5, 1), verbose = FALSE, scheme = "factorial", tol = 1e-25, sparse = TRUE) 

lapply(checklist, function(x){all.equal(factorial.wrap.result.sgcca[names(factorial.wrap.result.sgcca) == x], factorial.sgcca.nutrimouse.sgccak[names(factorial.sgcca.nutrimouse.sgccak) == x])})

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#note: the penalty parameters will need to be tuned

### centroid ###
centroid.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "centroid", verbose = FALSE, tol = 1e-25)
centroid.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(0.3, 0.5, 1), verbose = FALSE, scheme = "centroid", tol = 1e-25, sparse = TRUE) 

lapply(checklist, function(x){all.equal(centroid.wrap.result.sgcca[names(centroid.wrap.result.sgcca) == x], centroid.sgcca.nutrimouse.sgccak[names(centroid.sgcca.nutrimouse.sgccak) == x])})

### horst ###
horst.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "horst", verbose = FALSE, tol = 1e-25)
horst.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(0.3, 0.5, 1), verbose = FALSE, scheme = "horst", tol = 1e-25, sparse = TRUE) 

lapply(checklist, function(x){all.equal(horst.wrap.result.sgcca[names(horst.wrap.result.sgcca) == x], horst.sgcca.nutrimouse.sgccak[names(horst.sgcca.nutrimouse.sgccak) == x])})

### factorial ###
factorial.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE, tol = 1e-25)
factorial.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(0.3, 0.5, 1), verbose = FALSE, scheme = "factorial", tol = 1e-25, sparse = TRUE, max.iter = 1000) 

lapply(checklist, function(x){all.equal(factorial.wrap.result.sgcca[names(factorial.wrap.result.sgcca) == x], factorial.sgcca.nutrimouse.sgccak[names(factorial.sgcca.nutrimouse.sgccak) == x])})