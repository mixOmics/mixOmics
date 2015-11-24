library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}

sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE)


 ncomp = 2
toxicity.pls <- spls(X, Y, ncomp = ncomp, mode = 'regression',
keepX = c(rep(10, ncomp)), keepY = c(rep(3,ncomp)))

out=perf(toxicity.pls,validation="Mfold",folds=10)

# first, learn the model on the whole data set
#spls(X, Y,keepX=rep(10,10), ncomp = 10)


validation="Mfold"
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


#-- define the folds --#
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

#-- initialize new objects --#
RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
PRESS.inside = Q2 = MSEP = R2 = matrix(nrow = ncomp, ncol = q)
MSEP.mat = Ypred = array(0, c(n, q, ncomp))

press.mat = lapply(1 : ncomp, function(x){matrix(NA, nrow = n, ncol = q)})
RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = q)})
RSS.indiv[[1]] = X

#-- record feature stability --#
# initialize new objects:= to record feature stability
featuresX  = featuresY =  list()
for(k in 1:ncomp){
    featuresX[[k]] = featuresY[[k]] = NA
}

pb = txtProgressBar(style = 3)
nBar = 1


#-- loop on h = ncomp --#
for (h in 1:ncomp) {
    
    #-- initialising arguments --#
    tt = object$variates$X[, h]
    u = object$variates$Y[, h]
    b = object$loadings$Y[, h]
    c = object$mat.c[, h]
    d = object$mat.d[, h]
    nx = p - keepX[h]
    ny = q - keepY[h]
    
    RSS.indiv[[h + 1]] = Y - tt %*% t(d)
    RSS[h + 1, ] = colSums((Y - tt %*% t(d))^2)
    
    #-- loop on i (cross validation) --#
    for (i in 1:M)
    {
        if (progressBar == TRUE)
        {
            setTxtProgressBar(pb, nBar/(ncomp * M))
            nBar = nBar + 1
        }
        
        omit = folds[[i]]
        X.train = X[-omit, , drop = FALSE]
        Y.train = Y[-omit, , drop = FALSE]
        X.test = X[omit, , drop = FALSE]
        Y.test = Y[omit, , drop = FALSE]
        u.cv = u[-omit]
        
        #-- Q2 criterion
        a.old.cv = 0
        iter.cv = 1
        
        repeat{
            a.cv = crossprod(X.train, u.cv)
            if (nx != 0) {
                a.cv = ifelse(abs(a.cv) > abs(a.cv[order(abs(a.cv))][nx]),
                (abs(a.cv) - abs(a.cv[order(abs(a.cv))][nx])) * sign(a.cv), 0)
            }
            a.cv = a.cv / drop(sqrt(crossprod(a.cv)))
            t.cv = X.train %*% a.cv
            
            b.cv = crossprod(Y.train, t.cv)
            if (ny != 0) {
                b.cv = ifelse(abs(b.cv) > abs(b.cv[order(abs(b.cv))][ny]),
                (abs(b.cv) - abs(b.cv[order(abs(b.cv))][ny])) * sign(b.cv), 0)
            }
            b.cv = b.cv / drop(sqrt(crossprod(b.cv)))
            u.cv = Y.train %*% b.cv
            
            if ((crossprod(a.cv - a.old.cv) < tol) || (iter.cv == max.iter))
            break
            
            a.old.cv = a.cv
            iter.cv = iter.cv + 1
        }
        
        d.cv = t(Y.train) %*% (X.train %*% a.cv) / norm((X.train %*% a.cv), type = "2")^2
        Y.hat.cv = (X.test %*% a.cv) %*% t(d.cv)
        press.mat[[h]][omit, ] = Y.test - Y.hat.cv
        
        #-- for MSEP and R2 criteria
        if (h == 1)
        {
            nzv = (apply(X.train, 2, var) > .Machine$double.eps)
            spls.res = spls(X.train[, nzv],Y.train, ncomp = ncomp, mode = mode, max.iter = max.iter, tol = tol, keepX = keepX, keepY = keepY, near.zero.var = FALSE)
            
            X.test = X.test[, nzv, drop = FALSE]
            Y.hat = predict(spls.res, X.test)$predict
            
            for (k in 1:ncomp)
            {
                Ypred[omit, , k] = Y.hat[, , k]
                MSEP.mat[omit, , k] = (Y.test - Y.hat[, , k])^2
                
                # added: record selected features in each set
                featuresX[[k]] = c(unlist(featuresX[[k]]), selectVar(spls.res, comp = k)$name.X)
                featuresY[[k]] = c(unlist(featuresY[[k]]), selectVar(spls.res, comp = k)$name.Y)
                
            } # end loop on k
        }
    } # end i (cross validation)
    
    #-- compute the Q2 creterion --#
    PRESS.inside[h, ] = apply(press.mat[[h]], 2, function(x){norm(x, type = "2")^2})
    Q2[h, ] = 1 - PRESS.inside[h, ] / RSS[h, ]
    
    #-- deflation des matrices (for Q2 criterion)
    X = X - tt %*% t(c)
    Y = Y - tt %*% t(d)
    
    #-- compute the MSEP creterion --#
    MSEP[h, ] = apply(as.matrix(MSEP.mat[, , h]), 2, mean)
    
    #-- compute the R2 creterion --#
    R2[h, ] = (diag(cor(object$Y, Ypred[, , h])))^2
    
} #-- end loop on h --#

if (progressBar == TRUE) cat('\n')

#---- extract stability of features -----#
list.features.X = list()
list.features.Y = list()

for(k in 1:ncomp){
    #remove the NA value that was added for initialisation
    remove.naX = which(is.na(featuresX[[k]]))
    remove.naY = which(is.na(featuresY[[k]]))
    # then summarise as a factor and output the percentage of appearance
    list.features.X[[k]] = sort(table(as.factor(featuresX[[k]][-remove.naX]))/M, decreasing = TRUE)
    list.features.Y[[k]] = sort(table(as.factor(featuresY[[k]][-remove.naY]))/M, decreasing = TRUE)
    
}
names(list.features.X)  = names(list.features.Y) = paste('comp', 1:ncomp)

#-- output -----------------------------------------------------------------#
#---------------------------------------------------------------------------#
Q2.total = matrix(1 - rowSums(PRESS.inside) / rowSums(RSS[-(ncomp+1), , drop = FALSE]), nrow = 1, ncol = ncomp,
dimnames = list("Q2.total", paste0(1:ncomp, " comp")))

# set up dimnames
rownames(MSEP) = rownames(R2) = rownames(Q2) = paste0(1:ncomp, " comp")
colnames(MSEP) = colnames(R2) = colnames(Q2) = object$names$Y

res = list()
res$MSEP = t(MSEP)
res$R2 = t(R2)
res$Q2 = t(Q2)
res$Q2.total =  t(Q2.total)
res$RSS = RSS
res$PRESS = PRESS.inside
res$press.mat = press.mat
res$RSS.indiv = RSS.indiv

# features
res$features$stable.X = list.features.X
res$features$stable.Y = list.features.Y
