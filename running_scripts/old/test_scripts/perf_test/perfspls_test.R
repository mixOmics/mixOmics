library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

 ncomp = 7
toxicity.pls <- spls(X, Y, ncomp = ncomp, mode = 'regression',
keepX = c(rep(10, ncomp)), keepY = c(rep(3,ncomp)))

# first, learn the model on the whole data set
#spls(X, Y,keepX=rep(10,10), ncomp = 10)


validation="loo"
folds=10
progressBar = TRUE
near.zero.var=FALSE
object=toxicity.pls



featuresX  = featuresY =  list()
for(k in 1:ncomp){
    featuresX[[k]] = featuresY[[k]] = NA
}
#-- cross-validation approach  ---------------------------------------------#
#---------------------------------------------------------------------------#

#-- initialising arguments --#
# these are the centered and scaled matrices output from pls
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

if (any(is.na(X)) || any(is.na(Y)))
stop("missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.", call. = FALSE)


#-- tells which variables are selected in X and in Y --#
keepX = object$keepX
keepY = object$keepY

# -------------------------------------
# added: first check for near zero var on the whole data set
nzv = nearZeroVar(X)
if (length(nzv$Position > 0))
{
    warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
    X = X[, -nzv$Position, drop=TRUE]
    
    if(ncol(X)==0) {stop("No more predictors after Near Zero Var has been applied!")}
    
}
# and then we start from the X data set with the nzv removed



#-- M fold or loo cross validation --#
#- define the folds
if (validation == "Mfold") {
    if (is.list(folds)) {
        
        if (length(folds) < 2 || length(folds) > n)
        stop("Invalid number of folds.", call. = FALSE)
        
        if (length(unlist(folds)) != n)
        stop("Invalid folds. The total number of samples in folds must be equal to ",
        n, ".", call. = FALSE)
        
        if (length(unique(unlist(folds))) != n)
        stop("Invalid folds. Repeated samples in folds.", call. = FALSE)
        
        M = length(folds)
    } else {
        if (is.null(folds) || !is.finite(folds) || folds < 2 || folds > n)
        stop("Invalid number of folds.", call. = FALSE)
        else {
            M = round(folds)
            folds = split(sample(1:n), rep(1:M, length = n))
        }
    }
} else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
}

# =================================================================
#                           from DEVEL
# =================================================================
## add Q2
RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = q)})
RSS.indiv[[1]] = X
PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)


#KA: all criteria included in the computation

press.mat = Ypred = MSEP.mat = array(0, c(n, q, ncomp))
MSEP = R2 = matrix(0, nrow = q, ncol = ncomp)

# set up dimnames
rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
dimnames(press.mat)[[2]] = colnames(Y)

# in case the test set only includes one sample, it is better to advise the user to perform loocv
stop.user = FALSE
if (progressBar == TRUE) pb <- txtProgressBar(style = 3)



for (i in 1:M)
{
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
    
    #-- spls --#
    # added the near.zero.var option
    spls.res = spls(X.train, Y.train, ncomp=ncomp, mode=mode,
    max.iter=max.iter, tol=tol, keepX=keepX, keepY=keepY, near.zero.var = near.zero.var)
    
    # added: record selected features in each set
    for(k in 1:ncomp)
    {
        featuresX[[k]] = c(unlist(featuresX[[k]]), select.var(spls.res, comp = k)$name.X)
        featuresY[[k]] = c(unlist(featuresY[[k]]), select.var(spls.res, comp = k)$name.Y)
    }
    
    if (!is.null( spls.res$nzv$Position)) X.test = X.test[, - spls.res$nzv$Position,drop=FALSE]
    
    # in the predict function, X.test is already normalised w.r.t to training set X.train, so no need to do it here
    Y.hat = predict( spls.res, X.test)$predict #used for Ypred and MSEP.mat
    
    
    for (h in 1:ncomp)
    {
        Ypred[omit, , h] = Y.hat[, , h]
        # compute the MSEP
        MSEP.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
        
        # update 5.0-4, another way to calculate Y.hat is used for press.mat
        a.cv=spls.res$loadings$X[,h]
        d.cv = t(Y.train) %*% (X.train %*% a.cv) / norm((X.train %*% a.cv), type = "2")^2
        Y.hat.cv = (X.test %*% a.cv) %*% t(d.cv)
        press.mat[omit, , h] = (Y.test - Y.hat.cv)^2
        
    } # end h
} #end i (cross validation)


# these criteria are computed across all folds.
for (h in 1:ncomp)
{
    MSEP[, h] = apply(as.matrix(MSEP.mat[, , h]), 2, mean, na.rm = TRUE)
    R2[, h] = (diag(cor(Y, Ypred[, , h], use = "pairwise")))^2
    
    # RSS calculated on the whole data, nothing to do with the folds
    d = object$mat.d[, h]
    tt = object$variates$X[, h]
    RSS.indiv[[h + 1]] = Y - tt %*% t(d)
    RSS[h + 1, ] = colSums((Y - tt %*% t(d))^2)
    
    # on en profite pour calculer le PRESS par composante et le Q2
    if(q>1){
        PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
    }else{
        PRESS.inside[h, q] = sum(press.mat[,q, h], na.rm = TRUE)
    }
    
    Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
    
}

Q2.total = matrix(1 - rowSums(PRESS.inside) / rowSums(RSS[-(ncomp+1), , drop = FALSE]), nrow = 1, ncol = ncomp,
dimnames = list("Q2.total", paste0(1:ncomp, " comp")))

# set up dimnames
rownames(MSEP) = rownames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
colnames(MSEP) = colnames(R2) = colnames(Q2.inside) = object$names$Y

if (progressBar == TRUE) cat('\n')


#---- extract stability of features -----#
list.features.X = features.final.X = list()
list.features.Y = features.final.Y = list()

for(k in 1:ncomp){
    #remove the NA value that was added for initialisation
    remove.naX = which(is.na(featuresX[[k]]))
    remove.naY = which(is.na(featuresY[[k]]))
    # then summarise as a factor and output the percentage of appearance
    list.features.X[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
    list.features.Y[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
    
    #-- extract features selected from the full model
    features.final.X[[k]] = row.names(object$loadings$X)[object$loadings$X[, 1, drop = FALSE] != 0]
    features.final.Y[[k]] = row.names(object$loadings$Y)[object$loadings$Y[, 1, drop = FALSE] != 0]
}

names(features.final.X) = names(list.features.X) = names(features.final.Y) = names(list.features.Y) = paste('comp', 1:ncomp)



res = list()
res$MSEP = t(MSEP)
res$R2 = t(R2)
res$Q2 = t(Q2)
res$Q2.total =  t(Q2.total)
res$RSS = RSS
res$PRESS.inside = PRESS.inside
res$press.mat = press.mat
res$RSS.indiv = RSS.indiv

# features
res$features$stable.X = list.features.X
res$features$stable.Y = list.features.Y
res$features$final.X = features.final.X
res$features$final.Y = features.final.Y
res$nzvX = nzv$Position

res$Q2.total