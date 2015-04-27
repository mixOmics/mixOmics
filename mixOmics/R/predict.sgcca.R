# ----------------------------------------------------------------------------------------------------------
# predict.sgcca - Classify subejcts of the test data set using training model
#   inputs: object - object obtain from running wrapper.sGCCA()
#           response.train - response class of the training data
#           newdata - test dataset
#           g - index of dataset (gth dataset used to make prediction)
# ----------------------------------------------------------------------------------------------------------
predict.sgcca=function(object = NULL, newdata = NULL, method = NULL) {
    
    if (is.null(object) || is.null(newdata))
    stop("New data or object missing.")
    
    X = object$X; Y = object$Y[[1]]
    q = ncol(Y); p = lapply(X, ncol)
    ncomp = object$ncomp
    means.Y = attr(Y, "scaled:center"); sigma.Y = attr(Y, "scaled:scale")
    J = length(X)
    
    variatesX = object$variates[-(J + 1)]; loadingsX = object$loadings[-(J + 1)]
    
    ### Start: Check input newdata
    if (!is.list(newdata))
    newdata = list(newdata = newdata)
    
    newdata = lapply(newdata, function(x) {as.matrix(x, rownames.force = ifelse(is.null(row.names(x)), FALSE, TRUE))})
    
    if (any(lapply(newdata, function(x){length(dim(x))}) != 2)) {
        if (any(unlist(lapply(newdata, ncol)) != unlist(p)))
        stop("'newdata' must be a list with ", length(p), " numeric matrix and ncol respectively equal to ", paste(p, collapse = ", "), ".")
    }
    
    if (any(lapply(newdata, function(x){length(dim(x))}) == 0)) {
        if (any(unlist(lapply(newdata, ncol)) != unlist(p)))
        stop("'newdata' must be a list with ", length(p), " numeric matrix and ncol respectively equal to ", paste(p, collapse = ", "), ".")
        dim(newdata) = c(1, p)
    }
    
    if (any(sapply(newdata, function(x){any(is.na(x))})))
    stop("Some missing data are present in the matrix")
    ### End: Check input newdata
    
    ### Start: Check input method
    if (is.null(attr(Y,"levels")) & !is.null(method)) {
        warning("'method' available only at discriminant analysis level")
    }
    
    if (!is.null(attr(Y,"levels")))
    if ((method != "all") & (method != "max.dist") & (method != "centroids.dist") & (method != "mahalanobis.dist") || is.null(method))
    stop("ERROR : choose one of the four following modes: all, max.dist, centroids.dist or mahalanobis.dist")
    ### End: Check input method
    
    ### Start scaling and/or centering newdata
    if (object$scale == TRUE){
        if (!is.null(attr(X[[1]], "scaled:center")))
        newdata = lapply(1:J, function(x){sweep(newdata[[x]], 2, STATS = attr(X[[x]], "scaled:center"))})
        if (!is.null(attr(X[[1]], "scaled:scale")))
        newdata = lapply(1:J, function(x){sweep(newdata[[x]], 2, FUN = "/", STATS = attr(X[[x]], "scaled:scale"))})
    } else {
        stop("The object data are not scaled. Retry with the option: scale = TRUE.")
    }
    ### End scaling and/or centering newdata
    
    B.hat = t.pred = Y.hat = list() #= betay
    for (i in 1 : J){
        Pmat = Cmat = Wmat = NULL
        
        ### Start estimation using formula Y = XW(P'W)C (+ Yr, residuals on Y) See page 136 La regression PLS Theorie et pratique Tenenhaus
        # Estimation matrix W, P and C
        Pmat = crossprod(X[[i]], variatesX[[i]])
        Cmat = crossprod(Y, variatesX[[i]])
        Wmat = loadingsX[[i]]
        
        # Prediction Y.hat, B.hat and t.pred
        Ypred = lapply(1 : ncomp[i], function(x){newdata[[i]] %*% Wmat[, 1:x] %*% solve(t(Pmat[, 1:x]) %*% Wmat[, 1:x]) %*% t(Cmat)[1:x, ]})
        Ypred = sapply(Ypred, function(x){scale(scale(x, center = FALSE, scale = 1/sigma.Y), center = - means.Y, scale = FALSE)}, simplify = "array")
        
        Y.hat[[i]] = Ypred
        
        t.pred[[i]] = newdata[[i]] %*% Wmat %*% solve(t(Pmat) %*% Wmat)
        t.pred[[i]] = matrix(data = sapply(1:ncol(t.pred[[i]]),
        function(x) {t.pred[[i]][, x] * apply(variatesX[[i]], 2,
            function(y){(norm(y, type = "2"))^2})[x]}), nrow = nrow(newdata[[i]]), ncol = ncol(t.pred[[i]]))
        
        B.hat[[i]] = sapply(1 : ncomp[i], function(x){Wmat[, 1:x] %*% solve(t(Pmat[, 1:x]) %*% Wmat[, 1:x]) %*% t(Cmat)[1:x, ]}, simplify = "array")
        
        # betay[[i]] = Cmat
        ### End estimation using formula Y = XW(P'W)C (+ Yr, residuals on Y) See page 136 La regression PLS Theorie et pratique Tenenhaus
    }
    
    ### Start if discriminant analysis is performed
    if (!is.null(attr(Y,"levels"))) {
        Y.prim = object$ind.mat
        G = cls = list()
        for (i in 1 : J){
            G[[i]] = sapply(1:q, function(x) {apply(as.matrix(variatesX[[i]][Y.prim[, x] == 1,]), 2, mean)})
            if (ncomp[i] == 1)
            G[[i]] = t(t(G[[i]]))
            else
            G[[i]] = t(G[[i]])
        }
        
        ### Start: Maximum distance
        if (any(method == "all") || any(method == "max.dist")) {
            cls$max.dist = lapply(1:J, function(x){matrix(sapply(1:ncomp[x], ### List level
                function(y){apply(Y.hat[[x]][, , y, drop = FALSE], 1,  ### component level
                    function(z){which(z == max(z))}) ### matrix level
                }), nrow = nrow(newdata[[x]]), ncol = ncomp[x])
            })
            cls$max.dist = lapply(1:J, function(x){colnames(cls$max.dist[[x]]) = paste(rep("comp", ncomp[x]), 1 : ncomp[[x]], sep = " "); return(cls$max.dist[[x]])})
        }
        ### End: Maximum distance
        
        ### Start: Centroids distance
        if (any(method == "all") || any(method == "centroids.dist")) {
            cl = list()
            centroids.fun = function(x, G, h, i) {
                q = nrow(G[[i]])
                x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
                if (h > 1) {
                    d = apply((x - G[[i]][, 1:h])^2, 1, sum)
                }
                else {
                    d = (x - G[[i]][, 1])^2
                }
                cl.id = which.min(d)
            }
            
            for (i in 1 : J) {
                cl[[i]] = matrix(nrow = nrow(newdata[[1]]), ncol = ncomp[i])
                
                for (h in 1 : ncomp[[i]]) {
                    cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 1, function(x) {centroids.fun(x = x, G = G, h = h, i = i)})
                    cl[[i]][, h] = cl.id
                }
            }
            
            cls$centroids.dist = lapply(1:J, function(x){colnames(cl[[x]]) = paste(rep("comp", ncomp[x]), 1 : ncomp[[x]], sep = " "); return(cl[[x]])})
        }
        ### End: Centroids distance
        
        ### Start: Mahalanobis distance
        if (any(method == "all") || any(method == "mahalanobis.dist")) {
            cl = list()
            Sr.fun = function(x, G, Yprim, h, i) {
                q = nrow(G[[i]])
                Xe = Yprim %*% G[[i]][, 1:h]
                #Xr = object$variates$X[, 1:h] - Xe
                Xr = variatesX[[i]][, 1:h] - Xe
                Sr = t(Xr) %*% Xr/nrow(Y)
                Sr.inv = solve(Sr)
                x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
                if (h > 1) {
                    mat = (x - G[[i]][, 1:h]) %*% Sr.inv %*% t(x - G[[i]][, 1:h])
                    d = apply(mat^2, 1, sum)
                } else {
                    d = drop(Sr.inv) * (x - G[[i]][, 1])^2
                }
                cl.id = which.min(d)
            }
            
            for (i in 1 : J){
                cl[[i]] = matrix(nrow = nrow(newdata[[1]]), ncol = ncomp[i])
                
                for (h in 1:ncomp[[i]]) {
                    cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 1, Sr.fun, G = G, Yprim = Y.prim, h = h, i = i)
                    cl[[i]][, h] = cl.id
                }
            }
            
            cls$mahalanobis.dist = lapply(1:J, function(x){colnames(cl[[x]]) = paste(rep("comp", ncomp[x]), 1 : ncomp[[x]], sep = " "); return(cl[[x]])})
        }
        ### End: Mahalanobis distance  
    }
    ### End if discriminant analysis is performed  
    
    if (!is.null(X))
    names(Y.hat) = names(t.pred) = names(X) = names(B.hat) #names(betay) = 
    
    if (is.null(attr(Y,"levels")))
    return(invisible(list(predict = Y.hat, variates = t.pred, B.hat = B.hat)))
    else {
        return(invisible(list(predict = Y.hat, variates = t.pred, B.hat = B.hat, centroids = G, method = method, class = cls)))
    }
}