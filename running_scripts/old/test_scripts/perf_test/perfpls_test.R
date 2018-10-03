library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

data(nutrimouse)
X <- nutrimouse$gene
Y <- nutrimouse$lipid#[,1]
sp=spls(X,Y,keepX=c(10,10,10,10),ncomp=4)

ind=which(sp$loadings$X[,1]!=0)
ind=c(ind,which(sp$loadings$X[,2]!=0))
ind=c(ind,which(sp$loadings$X[,3]!=0))



toxicity.pls <- pls(X[,ind], Y, ncomp = 10)
#toxicity.pls <- pls(X, Y, ncomp = 10)

validation="loo"
folds=10
progressBar = TRUE
near.zero.var=FALSE
object=toxicity.pls

X = object$X
Y = object$Y

tol = object$tol
max.iter = object$max.iter
mode = object$mode
ncomp = object$ncomp
n = nrow(X)
p = ncol(X)
q = ncol(Y)
res = list()

nzv = nearZeroVar(X)
if (length(nzv$Position > 0))
{
    warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
    X = X[, -nzv$Position, drop=FALSE]
    
    if(ncol(X)==0) {stop("No more predictors after Near Zero Var has been applied!")}
    
}


#-- define the folds --#
if (validation == "Mfold") {
    if (is.list(folds)) {
        
        if (length(folds) < 2 || length(folds) > n)
        stop("Invalid number of folds.", call. = FALSE)
        
        if (length(unlist(folds)) != n)
        stop("Invalid folds. The total number of samples in folds must be equal to ", n, ".", call. = FALSE)
        
        if (length(unique(unlist(folds))) != n)
        stop("Invalid folds. Repeated samples in folds.", call. = FALSE)
        
        M = length(folds)
    } else {
        if (is.null(folds) || !is.finite(folds) || folds < 2 || folds > n) {
            stop("Invalid number of folds.", call. = FALSE)
        } else {
            M = round(folds)
            folds = split(sample(1:n), rep(1:M, length = n))
        }
    }
} else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
}

## add Q2
RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
RSS.indiv = array(0, c(n, q, ncomp+1))
PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)

RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = q)})
RSS.indiv[[1]] = X

press.mat = Ypred = MSEP.mat = array(0, c(n, q, ncomp))
MSEP = R2 = matrix(0, nrow = q, ncol = ncomp)

# set up dimnames
rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
dimnames(press.mat)[[2]] = colnames(Y)

# in case the test set only includes one sample, it is better to advise the user to perform loocv
stop.user = FALSE
if (progressBar == TRUE) pb <- txtProgressBar(style = 3)

for (i in 1:M) {
    if (progressBar == TRUE) setTxtProgressBar(pb, i/M)
    
    omit = folds[[i]]
    # see below, we stop the user if there is only one sample drawn on the test set using MFold
    if(length(omit) == 1) stop.user = TRUE
    
    # the training set is NOT scaled
    X.train = X[-omit, ]
    Y.train = Y[-omit, ]
    # the test set is scaled either in the predict function directly (for X.test)
    # or below for Y.test
    X.test = matrix(X[omit, ], nrow = length(omit))
    Y.test = matrix(Y[omit, ], nrow = length(omit))
    
    
    #-- pls --#
    pls.res = pls(X = X.train, Y = Y.train, ncomp = ncomp,
    mode = mode, max.iter = max.iter, tol = tol, near.zero.var = near.zero.var)
    
    if (!is.null( pls.res$nzv$Position)) X.test = X.test[, - pls.res$nzv$Position,drop=FALSE]
    
    # in the predict function, X.test is already normalised w.r.t to training set X.train, so no need to do it here
    Y.hat = predict( pls.res, X.test)$predict #used for Ypred and MSEP.mat
    
    #update 5.0-4, another way to calculate Y.hat is used for press.mat
    
    for (h in 1:ncomp) {
        Ypred[omit, , h] = Y.hat[, , h]
        
        # compute the press and the RSS
        # by definition (tenenhaus), RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
        #RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
        
        #from IG, RSS calculated on the whole data, nothing to do with the folds
        if(TRUE)
        {
            d = object$mat.d[, h]
            tt = object$variates$X[, h]
            RSS.indiv[[h + 1]] = Y - tt %*% t(d)
            RSS[h + 1, ] = colSums((Y - tt %*% t(d))^2)
        }
        
        
        MSEP.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
        
        a.cv=pls.res$loadings$X[,h]
        d.cv = t(Y.train) %*% (X.train %*% a.cv) / norm((X.train %*% a.cv), type = "2")^2
        Y.hat.cv = (X.test %*% a.cv) %*% t(d.cv)
        press.mat[omit, , h] = (Y.test - Y.hat.cv)^2
        #        press.mat[omit, , h] = (Y.test - Y.hat[,,h])^2
        
    } # end h
} #end i (cross validation)


# these criteria are computed across all folds.
for (h in 1:ncomp) {
    MSEP[, h] = apply(as.matrix(MSEP.mat[, , h]), 2, mean, na.rm = TRUE)
    R2[, h] = (diag(cor(Y, Ypred[, , h], use = "pairwise")))^2
    
    # on en profite pour calculer le PRESS par composante et le Q2
    if(q>1){
        #RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
        PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
    }else{
        #RSS[h+1,q] = sum(RSS.indiv[,q,h+1])
        PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
    }
    
    Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
    
}

Q2.total = matrix(1 - rowSums(PRESS.inside) / rowSums(RSS[-(ncomp+1), , drop = FALSE]), nrow = 1, ncol = ncomp,
dimnames = list("Q2.total", paste0(1:ncomp, " comp")))

if(FALSE)
{
colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")

if (q == 1){
    rownames(MSEP) = rownames(R2) = ""
}

if(ncomp>1){
    # compute Q2 total
    if(q>1){
        Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
    }else{ # for q == 1
        Q2.total = t(1 - PRESS.inside/RSS[-(ncomp+1), ])
    }
}else{
    Q2.total = NA
}

Q2.total
}

#a=perf.pls(toxicity.pls,validation="loo")
