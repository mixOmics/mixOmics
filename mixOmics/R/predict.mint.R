# Copyright (C) 2015
# Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# created 27-05-2015
# last modified 18-02-2016

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


# -------------------------- for plsda and splsda ---------------------------------------
predict.mint.pls <- predict.mint.spls <- predict.mint.plsda <-predict.mint.splsda<-
function(object, newdata,study.test,
method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), ...)
#method = "all", ...)

{
    #-- validation des arguments --#
    if (missing(newdata))
    stop("No new data available.")
    
    X = object$X
    Y = object$Y
    study.learn=object$study
    concat.Y.init=object$concat.Y.init

    if(any(class(object)%in%c("mint.plsda","mint.splsda")))
    {
        q = ncol(object$ind.mat)
    }else{
        q=ncol(object$Y)
    }
    p = ncol(X)
    
    if (length(dim(newdata)) == 2) {
        if (ncol(newdata) != p)
        stop("'newdata' must be a numeric matrix with ncol = ", p,
        " or a vector of length = ", p, ".")
    }
    
    if (length(dim(newdata)) == 0) {
        if (length(newdata) != p)
        stop("'newdata' must be a numeric matrix with ncol = ", p,
        " or a vector of length = ", p, ".")
        dim(newdata) = c(1, p)
    }
    
    X.test=newdata
    n=nrow(X.test)
    if(length(rownames(X.test))==0) rownames(X.test)=1:nrow(X.test)
    if(max(table(rownames(X.test)))>1) stop('samples should have a unique identifier/rowname')

	X.test.list.study = study_split(X.test,study.test)$data.list.study
  	study.name = split(rownames(X.test),study.test)
  	M=length(levels(study.test))
    P=ncol(X.test)
  	
   
    
    #  for(m in 1:M){
    #
    # 	X.test.list.study[[m]]=matrix(X.test.list.study[[m]], ncol=P)
    # 	colnames(X.test.list.study[[m]]) = colnames(X.test)
    # 	rownames(X.test.list.study[[m]]) = study.name[[m]]
   	#   }
    
    scale =object$normaliz$scale

    means.Y=object$normaliz$means.Y
    means.Y.study=apply(means.Y,1,mean)
    intercept.Y=matrix(0,nrow=n,ncol=q)

    # each study is normalized, dependin on how the normalization was performed on the learning set (center/scale)
    X.test.list.study.scale = list()
    concat.X.test  = NULL
    for (m in 1:M)
    {
    	# case 1: sample test is part of the learning data set
        if(any(levels(study.test)[m]==levels(study.learn))) #if some study of study.test were already in the learning set
        {
            
               ind=which(levels(study.test)[m]==levels(study.learn))
            
            if(scale==TRUE)
    		X.test.list.study.scale[[m]]=scale(X.test.list.study[[m]], center = object$normaliz$means.X[ind,], scale = object$normaliz$sigma.X[m,])

            if(scale==FALSE)
            X.test.list.study.scale[[m]]=scale(X.test.list.study[[m]], center = object$normaliz$means.X[ind,],scale=FALSE)
            
            intercept.Y[which(study.test%in%levels(study.learn)[m]),]=means.Y.study[m]

            
        }else{ 
        # case 2: sample test is a new study # a new study which was not in the learning set
            X.test.list.study.scale[[m]]=scale(X.test.list.study[[m]],center=TRUE,scale=scale)
        }
        concat.X.test = rbind(concat.X.test, unlist(X.test.list.study.scale[[m]]))
    }
    #print(intercept.Y)
    #sort the samples of concat.X.test to match the ones of X.test
    indice.match=match(rownames(X.test),rownames(concat.X.test))
    concat.X.test[which(is.na(concat.X.test))]=0 # taking care of the NA due to normalisation: put to 0 so they don't influence the product below (in Y.hat), reviens au meme que de supprimer les colonnes sans variances au depart.
    

    #-- initialisation des matrices --#
    ncomp = object$ncomp
    c = object$mat.c
    a = object$loadings.global$X
    concat.Y.init= object$concat.Y.init
     variates.global.X = object$variates.global$X
    
    betay = list()
   
    newdata = as.matrix(concat.X.test)
    ones = matrix(rep(1, nrow(newdata)), ncol = 1)
    ##- coeff de regression
    B.hat = array(0, dim = c(p, q, ncomp))
    ##- prediction
    Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
    Y.hat.2 = array(0, dim = c(nrow(newdata), q, ncomp))
    ##- variates
    t.pred = array(0, dim = c(nrow(newdata), ncomp))
    
    #-- calcul de la prediction --#
    for(h in 1:ncomp)
    {
        # compute coefficients matrix B
        dd= coefficients(lm(concat.Y.init~variates.global.X[,1:h,drop=FALSE])) #regression of concat.Y.init (Y centered by group) on variates.global.X => =loadings.global.Y at a scale factor
         if(q==1){betay[[h]]=(dd[-1])}
         if(q>=2){betay[[h]]=(dd[-1,])}
        
        W = a[, 1:h,drop=FALSE] %*% solve(t(c[, 1:h,drop=FALSE]) %*% a[, 1:h,drop=FALSE])
        B = W %*% drop(betay[[h]])
        
        #intercept=
        

Y.hat[, , h] = concat.X.test %*% as.matrix(B) +intercept.Y#+ ones %*% intercept     #so far: gives a prediction of Ycentered, need to add the Ybar(=means.Y) (cf these aida)
Y.hat.2[, , h] = concat.X.test %*% as.matrix(B) #+ ones %*% intercept     #so far: gives a prediction of Ycentered, need to add the Ybar(=means.Y) (cf these aida)
        t.pred[, h] = concat.X.test %*% W[, h]
        
        B.hat[, , h] = B
    
    #try to see if lm gives the same as loadgins.global.Y => seems NOOO
    #  betay2=matrix(object$loadings.global$Y[,h],nrow=1)
    #   B = W %*% drop(betay2)
    #   Y.hat.2[, , h] = concat.X.test %*% as.matrix(B) #+ ones %*% intercept     #so far: gives a prediction of Ycentered, need to add the Ybar(=means.Y) (cf these aida)


    } # end h

    
    #match indice
    concat.X.test=concat.X.test[indice.match,]
    Y.hat=Y.hat[indice.match,,,drop=FALSE]
    Y.hat.2=Y.hat.2[indice.match,,,drop=FALSE]
    newdata=newdata[indice.match,]
       
    #     for(h in 1:ncomp){
    #       W = a[, 1:h] %*% solve(t(c[, 1:h]) %*% a[, 1:h])
    #       B = W %*% drop(t(b[, 1:h]))
    #       B = scale(B, center = FALSE, scale = 1 / sigma.Y)
    #       B = as.matrix(scale(t(B), center = FALSE, scale = sigma.X))
    #       intercept = -scale(B, center = FALSE, scale = 1 / means.X)
    #       intercept = matrix(apply(intercept, 1, sum) + means.Y, nrow = 1)
    #       Y.hat[, , h] = newdata %*% t(B) + ones %*% intercept
    #       t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[, h]
    #       B.hat[, , h] = B
    #     }  #end h
    #

    # ----  if Y is a factor  -----------------
    
    if(any(class(object)%in%c("mint.plsda","mint.splsda")))
    {
        # ----    max distance -----------------
        cls = list()
        if (any(method == "all") || any(method == "max.dist")) {
            
            function.pred = function(x){
                nr = nrow(x)
                tmp = vector("numeric", nr)
                for(j in 1:nr){
                    tmp[j] = (which(x[j, ] == max(x[j, ]))[1])
                }
                return(tmp)
            }
            cls$max.dist = matrix(apply(Y.hat, 3, function.pred), ncol = ncomp)
            colnames(cls$max.dist) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
        }
    
    
        G = matrix(0, nrow = q, ncol = ncomp)
        
        for (i in 1:q) {
            if(ncomp > 1) {
                G[i, ] = apply(object$variates.global$X[object$ind.mat[, i] == 1, ], 2, mean)
            }
            else {
                G[i, ] = mean(object$variates.global$X[object$ind.mat[, i] == 1, ])
            }
        }
        
        # ----    centroids distance -----------------
        
        if (any(method == "all") || any(method == "centroids.dist")) {
            
            cl = matrix(nrow = nrow(newdata), ncol = ncomp)
            
            centroids.fun = function(x, G, h) {
                q = nrow(G)
                x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
                if (h > 1) {
                    d = apply((x - G[, 1:h])^2, 1, sum)
                }
                else {
                    d = (x - G[, 1])^2
                }
                cl.id = which.min(d)
            }
            
            for (h in 1:ncomp) {
                cl.id = apply(matrix(t.pred[, 1:h], ncol = h), 1, centroids.fun, G = G, h = h)
                cl[, h] = cl.id
            }
            colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
            cls$centroids.dist = cl
        }
        
        # ----    mahalanobis distance -----------------
        
        if (any(method == "all") || any(method == "mahalanobis.dist")) {
            
            cl = matrix(nrow = nrow(newdata), ncol = ncomp)
            
            Sr.fun = function(x, G, concat.Y.init, h) {
                q = nrow(G)
                Xe = concat.Y.init %*% G[, 1:h]
                Xr = object$variates.global$X[, 1:h] - Xe
                Sr = t(Xr) %*% Xr / n
                Sr.inv = solve(Sr)
                x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
                if (h > 1) {
                    mat = (x - G[, 1:h]) %*% Sr.inv %*% t(x - G[, 1:h])
                    d = apply(mat^2, 1, sum)
                }
                else {
                    d = drop(Sr.inv) * (x - G[, 1])^2
                }
                cl.id = which.min(d)
            }
            
            for (h in 1:ncomp) {
                cl.id = apply(matrix(t.pred[, 1:h], ncol = h), 1, Sr.fun, G = G, concat.Y.init = object$ind.mat, h = h)
                cl[, h] = cl.id
            }
            colnames(cl) = paste(rep("comp", ncomp), 1:ncomp, sep = " ")
            cls$mahalanobis.dist = cl
        }
        
        
        #-- valeurs sortantes --#
        if (any(method == "all")) method = "all"
        rownames(t.pred) = rownames(newdata)
        colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
        rownames(Y.hat) = rownames(newdata)
        colnames(Y.hat) = colnames(Y)
        colnames(G) = paste("dim", c(1:ncomp), sep = " ")
        
        return(invisible(list(predict = Y.hat,
        variates = t.pred,
        B.hat = B.hat,
        centroids = G,
        method = method,
                              betay=betay,
        class = cls,
        X.test.list.study.scale=X.test.list.study.scale,
        concat.X.test=concat.X.test)))
        
    }else{    #end condition on Y is a factor
        
        
        #-- valeurs sortantes --#
        rownames(t.pred) = rownames(newdata)
        colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
        rownames(Y.hat) = rownames(newdata)
        colnames(Y.hat) = colnames(Y)
        
        return(invisible(list(predict = Y.hat,predict.2 = Y.hat.2,
        variates = t.pred,
        B.hat = B.hat,
        method = method,
        t.pred=t.pred,
                              betay=betay,
                              X.test.list.study.scale=X.test.list.study.scale,
                              concat.X.test=concat.X.test)))
        
        
    }

}
