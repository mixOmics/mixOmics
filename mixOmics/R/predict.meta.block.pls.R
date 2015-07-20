# Copyright (C) 2014
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Amrit Singh, University of British Columbia, Vancouver.
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
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

# created: 27/05/15
# modified: 29/05/15
# author: Florian Rohart

# This function makes a prediction of a 'newdata' by using the results of 'object'.
# Depending on the class of the object (meta).(block).(s)pls(da) (16 classes so far), the input data is different
# and the preparation of the data is different - different scaling for instance if object="meta...."
# However, the prediction formula is the same for all classes, thus only one code

predict.mixOmics <-
predict.pls <-  predict.spls<- predict.plsda <- predict.splsda <-
predict.meta.pls <- predict.meta.spls <- predict.meta.plsda <- predict.meta.splsda <-
predict.block.pls <- predict.block.spls <- predict.block.plsda <- predict.block.splsda <-
predict.meta.block.pls <- predict.meta.block.spls <- predict.meta.block.plsda <- predict.meta.block.splsda <-

function(object, newdata,study.test,method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),  ...)
{
    
    if(any(class(object)%in%c("rgcca","sparse.rgcca")))
    stop("no prediction for RGCCA methods")
    
    #check on method
    if (!any(method %in% c("all", "max.dist", "centroids.dist", "mahalanobis.dist")) & length(grep("plsda",class(object)))>0)
    stop("ERROR : choose one of the four following modes: 'all', 'max.dist', 'centroids.dist' or 'mahalanobis.dist'")
    #end general checks
    
    ncomp = object$ncomp
    newdata.input=newdata
    
    
    if(length(grep("plsda",class(object)))>0) # a DA analysis (meta).(block).(s)plsda
    {
        #if DA analysis, the unmap Y is in ind.mat
        Y.factor=object$Y
        Y=object$ind.mat
    }else{
        #if not DA, Y is in object$Y
        Y=object$Y
    }
    q=ncol(Y)
    

    ### if the object is a block, the input newdata is different, we check newdata, make sure it's a list and check newdata/X
    if(length(grep("block",class(object)))==0) # not a block (pls/spls/plsda/splsda/meta...)
    {
        p=ncol(object$X)
        if(is.list(object$X))
        stop("Something is wrong, object$X should be a matrix and it appears to be a list") #this should never happen/intern check
        
        if(is.list(newdata) & !is.data.frame(newdata))
        stop("'newdata' must be a numeric matrix")
        
        # deal with near.zero.var in object, to remove the same variable in newdata as in object$X (already removed in object$X)
        if(length(object$nzv$Position) > 0)
        newdata = newdata[, -object$nzv$Position,drop=FALSE]
        
        if(all.equal(colnames(newdata),colnames(object$X))!=TRUE)
        stop("'newdata' must include all the variables of 'object$X'")
        
        #not a block, the input newdata should be a matrix
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
        
        #check col/rownames of newdata
        check=Check.entry.single(newdata, ncomp,q=1)
        newdata=check$X
        
        if(length(rownames(newdata))==0) rownames(newdata)=1:nrow(newdata)
        if(max(table(rownames(newdata)))>1) stop('samples should have a unique identifier/rowname')
        
        # we transform everything in lists
        X=list(X=object$X)
        newdata=list(newdata=newdata)
        
    }else{

        # a block, newdata should be a list, each blocks should have the same number of samples
        if(!is.list(newdata))
        stop("'newdata' should be a list")
        
        X = object$X
        p = lapply(X, ncol)

        # deal with near.zero.var in object, to remove the same variable in newdata as in object$X (already removed in object$X)
        if(!is.null(object$nzv))
        {
            newdata = lapply(1:(length(object$nzv)-1),function(x){if(length(object$nzv[[x]]$Position>0)) {newdata[[x]][, -object$nzv[[x]]$Position,drop=FALSE]}else{newdata[[x]]}})
        }
        names(newdata)=names(X)

        if(length(newdata)!=length(object$X)) stop("'newdata' must have as many blocks as 'object$X'")


        
        if (any(lapply(newdata, function(x){length(dim(x))}) != 2)) {
            if (any(unlist(lapply(newdata, ncol)) != unlist(p)))
            stop("'newdata' must be a list with ", length(p), " numeric matrix and ncol respectively equal to ", paste(p, collapse = ", "), ".")
        }
        
        if (any(lapply(newdata, function(x){length(dim(x))}) == 0)) {
            if (any(unlist(lapply(newdata, ncol)) != unlist(p)))
            stop("'newdata' must be a list with ", length(p), " numeric matrix and ncol respectively equal to ", paste(p, collapse = ", "), ".")
            dim(newdata) = c(1, p) #don't understand that, Benoit?
        }
        
        if (any(sapply(newdata, function(x){any(is.na(x))})))
        stop("Some missing data are present in the matrix")
        
        
        #check dimnames and ncomp per block of A
        for(q in 1:length(newdata))
        {
            check=Check.entry.single(newdata[[q]], ncomp[q],q=q)
            newdata[[q]]=check$X
        }
        names(newdata)=names(X)
        #check that newdata and X have the same variables
        if(all.equal(lapply(newdata,colnames),lapply(X,colnames))!=TRUE)
        stop("Each 'newdata[[i]]' must include all the variables of 'object$X[[i]]'")
        
    }
    
    
    p = lapply(X, ncol)
    q = ncol(Y)
    J=length(X) #at this stage we have a list of blocks
    variatesX = object$variates[-(J + 1)]; loadingsX = object$loadings[-(J + 1)]
    
    scale=object$scale # X and Y are both mean centered by groups and if scale=TRUE they are scaled by groups

    ### if the object is not a meta analysis, the input study.test is missing and we can go faster to scale the data
    if(length(grep("meta",class(object)))==0 )#| nlevels(factor(object$study))<=1) #not a meta object or just one level in the study
    {   # not a meta (pls/spls/plsda/splsda/block...)

        # scale newdata if just one study
        if (!is.null(attr(X[[1]], "scaled:center")))
        newdata = lapply(1:J, function(x){sweep(newdata[[x]], 2, STATS = attr(X[[x]], "scaled:center"))})
        if (scale)
        newdata = lapply(1:J, function(x){sweep(newdata[[x]], 2, FUN = "/", STATS = attr(X[[x]], "scaled:scale"))})

        means.Y = matrix(attr(Y, "scaled:center"),nrow=nrow(newdata.input),ncol=q,byrow=TRUE);
        if (scale)
        {sigma.Y = matrix(attr(Y, "scaled:scale"),nrow=nrow(newdata.input),ncol=q,byrow=TRUE)}else{sigma.Y=matrix(1,nrow=nrow(newdata.input),ncol=q)}
        concat.newdata=newdata
        
    }else{
        # a meta analysis
        
        #check study.test
        if(missing(study.test))
        {
            if(nlevels(object$study)==1)
            {study.test=factor(rep(1,nrow(newdata[[1]])))}else{
                stop("'study.test' is missing")}
        }else{
            study.test=as.factor(study.test)
        }
        if (any(unlist(lapply(newdata, nrow)) != length(study.test)))
        {
            stop(paste0("'study' must be a factor of length ",nrow(X),"."))
        }
        
        #scale newdata if more than one study. If some levels of study.test are included in study.learn, the means/sigmas of study.learn are used to scale
        M=nlevels(study.test)
        study.learn=factor(object$study)
        names.study.learn=levels(study.learn);names.study.test=levels(study.test)
        match.study=match(names.study.test,names.study.learn)
        match.study.indice=which(!is.na(match.study))
        
        newdata.list.study = lapply(newdata,study_split,study.test) #a list of lists. [[j]][[m]]: j is for the blocks, m is for the levels(study)
        
        # each study is normalized, depending on how the normalization was performed on the learning set (center/scale)
        newdata.list.study.scale.temp =NULL#vector("list",length=J) #list(list())
        concat.newdata  = vector("list",length=J)
        for(j in 1:J) # loop on the blocks
        {
            for (m in 1:M) #loop on the groups of study.test
            {
                # case 1: sample test is part of the learning data set
                if(m%in%match.study.indice) #if some study of study.test were already in the learning set
                {
                    
                    if(scale==TRUE)
                    {
                        if(nlevels(object$study)>1)
                        {
                            newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]], center = attr(X[[j]],paste0("means:", levels(study.test)[m])), scale = attr(X[[j]],paste0("sigma:", levels(study.test)[m])))
                        }else{#if just one level in object$study, the normalisation attributes are not named the same
                            newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]], center = attr(X[[j]],"scaled:center"), scale = attr(X[[j]],"scaled:scale"))
                        }
                    }
                    
                    if(scale==FALSE)
                    {
                        if(nlevels(object$study)>1)
                        {
                            newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]], center = attr(X[[j]],paste0("means:", levels(study.test)[m])),scale=FALSE)
                        }else{#if just one level in object$study, the normalisation attributes are not named the same
                            newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]], center = attr(X[[j]],"scaled:center"),scale=FALSE)
                        }
                        
                        
                    }
                    
                }else{
                # case 2: sample test is a new study # a new study which was not in the learning set
                    newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]],center=TRUE,scale=scale)
                }
                #concatenation of the scaled data
                concat.newdata[[j]] = rbind(concat.newdata[[j]], unlist(newdata.list.study.scale.temp))#[[j]][[m]]))
            }
        }
        
        # we now need to reorganise concat.newdata as newdata. Indeed, the concatenation was done without taking care of the order of the samples in newdata
        for(j in 1:J) # loop on the blocks
        {
            indice.match=match(rownames(newdata[[j]]),rownames(concat.newdata[[j]]))
            #match indice
            concat.newdata[[j]]=concat.newdata[[j]][indice.match,]
            concat.newdata[[j]][which(is.na(concat.newdata[[j]]))]=0 # taking care of the NA due to normalisation: put to 0 so they don't influence the product below (in Y.hat), reviens au meme que de supprimer les colonnes sans variances au depart.
            
        }
        names(concat.newdata)=names(X)

        means.Y=matrix(0,nrow=nrow(concat.newdata[[1]]),ncol=q)
        sigma.Y=matrix(1,nrow=nrow(concat.newdata[[1]]),ncol=q)
        
        #loop on the blocks to define means.Y and sigma.Y for meta analysis
        for(m in 1:M)
        {
            if(m%in%match.study.indice) #if some study of study.test were already in the learning set
            {
                if(nlevels(object$study)>1)
                {
                    means.Y[which(study.test%in%levels(study.learn)[match.study[m]]),]=matrix(attr(Y,paste0("means:", levels(study.test)[m])),nrow=length(which(study.test%in%levels(study.learn)[match.study[m]])),ncol=q,byrow=TRUE)
                }else{#if just one level in object$study, the normalisation attributes are not named the same
                    means.Y[which(study.test%in%levels(study.learn)[match.study[m]]),]=matrix(attr(Y,"scaled:center"),nrow=length(which(study.test%in%levels(study.learn)[match.study[m]])),ncol=q,byrow=TRUE)
                }
                
                
                if(scale==TRUE)
                {
                    # I want to build a vector with sigma.Y for each group
                    if(nlevels(object$study)>1)
                    {
                        sigma.Y[which(study.test%in%levels(study.learn)[match.study[m]]),]=matrix(attr(Y,paste0("sigma:", levels(study.test)[m])),nrow=length(which(study.test%in%levels(study.learn)[match.study[m]])),ncol=q,byrow=TRUE)
                    }else{#if just one level in object$study, the normalisation attributes are not named the same
                        sigma.Y[which(study.test%in%levels(study.learn)[match.study[m]]),]=matrix(attr(Y,"scaled:scale"),nrow=length(which(study.test%in%levels(study.learn)[match.study[m]])),ncol=q,byrow=TRUE)
                    }
                    
                    
                }
            }
        }
        
    }
    ### end if object is a meta analysis
    
    ### at this stage we have
    # X         # list of blocks
    # Y         # observation
    # newdata   #list of blocks for the prediction, same length as A, scaled
    
    
    # -----------------------
    #       prediction
    # -----------------------
    
    B.hat = t.pred = Y.hat = list() #= betay
    for (i in 1 : J)
    {
        Pmat = Cmat = Wmat = NULL
        
        ### Start estimation using formula Y = XW(P'W)C (+ Yr, residuals on Y) See page 136 La regression PLS Theorie et pratique Tenenhaus
        # Estimation matrix W, P and C
        Pmat = crossprod(X[[i]], variatesX[[i]])
        Cmat = crossprod(Y, variatesX[[i]])
        Wmat = loadingsX[[i]]
        
        # Prediction Y.hat, B.hat and t.pred
        Ypred = lapply(1 : ncomp[i], function(x){concat.newdata[[i]] %*% Wmat[, 1:x] %*% solve(t(Pmat[, 1:x]) %*% Wmat[, 1:x]) %*% t(Cmat)[1:x, ]})
        Ypred = sapply(Ypred, function(x){x*sigma.Y + means.Y}, simplify = "array")
        
        Y.hat[[i]] = Ypred
        
        t.pred[[i]] = concat.newdata[[i]] %*% Wmat %*% solve(t(Pmat) %*% Wmat)
        t.pred[[i]] = matrix(data = sapply(1:ncol(t.pred[[i]]),
        function(x) {t.pred[[i]][, x] * apply(variatesX[[i]], 2,
            function(y){(norm(y, type = "2"))^2})[x]}), nrow = nrow(concat.newdata[[i]]), ncol = ncol(t.pred[[i]]))
        
        B.hat[[i]] = sapply(1 : ncomp[i], function(x){Wmat[, 1:x] %*% solve(t(Pmat[, 1:x]) %*% Wmat[, 1:x]) %*% t(Cmat)[1:x, ]}, simplify = "array")
        ### End estimation using formula Y = XW(P'W)C (+ Yr, residuals on Y) See page 136 La regression PLS Theorie et pratique Tenenhaus
        
        rownames(t.pred[[i]]) = rownames(newdata[[i]])
        colnames(t.pred[[i]]) = paste("dim", c(1:ncomp[i]), sep = " ")
        rownames(Y.hat[[i]]) = rownames(newdata[[i]])
        colnames(Y.hat[[i]]) = colnames(Y)
        dimnames(Y.hat[[i]])[[3]]=paste("dim", c(1:ncomp[i]), sep = " ")
        rownames(B.hat[[i]]) = colnames(newdata[[i]])
        colnames(B.hat[[i]]) = colnames(Y)
        dimnames(B.hat[[i]])[[3]]=paste("dim", c(1:ncomp[i]), sep = " ")

    }
    

    #-- valeurs sortantes --#
    names(Y.hat)=names(t.pred)=names(B.hat)=names(object$X)

    
    # basic prediction results
    if(length(grep("block",class(object)))!=0 )
    {
        out=list(predict=Y.hat,variates=t.pred,B.hat=B.hat)
        out$newdata=concat.newdata
    }else{# not a block (pls/spls/plsda/splsda/meta...)
        out=list(predict=Y.hat[[1]],variates=t.pred[[1]],B.hat=B.hat[[1]])
        out$newdata=concat.newdata[[1]]

    }
    
    ### if the object is a DA analysis, we gives the class of each sample depending on 'method'
    if(length(grep("plsda",class(object)))>0) # a DA analysis (meta).(block).(s)plsda
    {
        Y.prim = unmap(object$Y)
        G = cls = list()
        for (i in 1 : J)
        {
            G[[i]] = sapply(1:q, function(x) {apply(as.matrix(variatesX[[i]][Y.prim[, x] == 1,,drop=FALSE]), 2, mean)})
            if (ncomp[i] == 1)
            G[[i]] = t(t(G[[i]]))
            else
            G[[i]] = t(G[[i]])
            colnames(G[[i]]) = paste("dim", c(1:ncomp[i]), sep = " ")

        }
        names(G)=names(object$X)
        
        ### Start: Maximum distance
        if (any(method == "all") || any(method == "max.dist"))
        {
            cls$max.dist = lapply(1:J, function(x){matrix(sapply(1:ncomp[x], ### List level
                function(y){apply(Y.hat[[x]][, , y, drop = FALSE], 1,  ### component level
                    function(z){which(z == max(z))}) ### matrix level
                }), nrow = nrow(newdata[[x]]), ncol = ncomp[x])
            })
            cls$max.dist = lapply(1:J, function(x){colnames(cls$max.dist[[x]]) = paste(rep("comp", ncomp[x]), 1 : ncomp[[x]], sep = " ");
                rownames(cls$max.dist[[x]]) = rownames(newdata[[x]]); return(cls$max.dist[[x]])})
            names(cls$max.dist)=names(object$X)
        }
        
        
        ### Start: Centroids distance
        if (any(method == "all") || any(method == "centroids.dist"))
        {
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
            
            for (i in 1 : J)
            {
                cl[[i]] = matrix(nrow = nrow(newdata[[1]]), ncol = ncomp[i])
                
                for (h in 1 : ncomp[[i]])
                {
                    cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 1, function(x) {centroids.fun(x = x, G = G, h = h, i = i)})
                    cl[[i]][, h] = cl.id
                }
            }
            
            cls$centroids.dist = lapply(1:J, function(x){colnames(cl[[x]]) = paste(rep("comp", ncomp[x]), 1 : ncomp[[x]], sep = " ");
                rownames(cl[[x]]) = rownames(newdata[[x]]); return(cl[[x]])})
            names(cls$centroids.dist)=names(object$X)
        }### End: Centroids distance
        
        
        ### Start: Mahalanobis distance
        if (any(method == "all") || any(method == "mahalanobis.dist"))
        {
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
            
            cls$mahalanobis.dist = lapply(1:J, function(x){colnames(cl[[x]]) = paste(rep("comp", ncomp[x]), 1 : ncomp[[x]], sep = " ");
                rownames(cl[[x]]) = rownames(newdata[[x]]);return(cl[[x]])})
            names(cls$mahalanobis.dist)=names(object$X)
        } ### End: Mahalanobis distance
        
        out$class=cls
        
        if(length(grep("block",class(object)))!=0 ) # a block
        {
            out$centroids = G
        }else{ #not a block
            out$centroids = G[[1]]
            out$class=lapply(out$class,function(x){x[[1]]})
        }
        if (any(method == "all")) method = "all"
        out$method = method
    }### End if discriminant analysis is performed
    
    
    
    out
    
}





