#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 24-08-2016
# last modified: 05-10-2017
#
# Copyright (C) 2015
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
# tune.splsda: chose the optimal number of parameters per component on a splsda method
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# multilevel: repeated measurement only
# validation: Mfold or loo cross validation
# folds: if validation = Mfold, how many folds?
# nrepeat: number of replication of the Mfold process
# ncomp: the number of components to include in the model. Default to 1.
# choice.keepX: a vector giving keepX on the components that were already tuned
# test.keepX: grid of keepX among which to chose the optimal one
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# dist: distance to classify samples. see predict
# auc:
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# progressBar: show progress,
# class.object
# cl: if parallel, the clusters
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# misdata: optional. any missing values in the data? list, misdata[[q]] for each data set
# is.na.A: optional. where are the missing values? list, is.na.A[[q]] for each data set (if misdata[[q]] == TRUE)
# parallel: logical.


MCVfold.spls = function(
X,
Y,
multilevel = NULL, # repeated measurement only
validation,
folds,
nrepeat = 1,
ncomp,
choice.keepX = NULL, # keepX chosen on the first components
test.keepX, # a vector of value(keepX) to test on the last component. There needs to be names(test.keepX)
measure = "MSE", # MAE (Mean Absolute Error: MSE without the square), Bias (average of the differences)
dist = "max.dist",
auc = FALSE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
progressBar = TRUE,
class.object,
cl,
scale,
misdata,
is.na.A,
parallel
)
{
    auc=FALSE
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    #-- set up a progress bar --#
    if (progressBar ==  TRUE)
    {
        pb = txtProgressBar(style = 3)
        nBar = 1
    } else {
        pb = FALSE
    }
    
    design = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)
    
    if(ncomp>1)
    {
        keepY = rep(ncol(Y), ncomp-1)
    } else {keepY = NULL}
    
    M = length(folds)
    features = features.j = NULL
    auc.all = prediction.comp = class.comp = list()
    class.comp[[1]] = array(0, c(nrow(X), nrepeat, length(test.keepX)))# prediction of all samples for each test.keepX and  nrep at comp fixed
    prediction.keepX = vector("list", length=length(test.keepX))
    
    for(i in 1:length(test.keepX))
        prediction.keepX[[i]] = array(0,c(nrow(X), ncol(Y), nrepeat), dimnames = list(rownames(X), colnames(Y), paste0("nrep.", 1:nrepeat)))

    folds.input = folds
    for(nrep in 1:nrepeat)
    {
        prediction.comp[[nrep]] = array(0, c(nrow(X), ncol(Y), length(test.keepX)), dimnames = list(rownames(X), colnames(Y), names(test.keepX)))
        rownames(prediction.comp[[nrep]]) = rownames(X)
        colnames(prediction.comp[[nrep]]) = colnames(Y)
        
        
        n = nrow(X)
        repeated.measure = 1:n
        if (!is.null(multilevel))
        {
            repeated.measure = multilevel[,1]
            n = length(unique(repeated.measure)) # unique observation: we put every observation of the same "sample" in the either the training or test set
        }
        
        
        #-- define the folds --#
        if (validation ==  "Mfold")
        {
            
            if (nrep > 1) # reinitialise the folds
            folds = folds.input
            
            if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
            {
                stop("Invalid number of folds.")
            } else {
                M = round(folds)
                if (is.null(multilevel))
                {
                    folds = suppressWarnings(split(sample(1:n), rep(1:n, length = M)))
                } else {
                    folds = split(sample(1:n), rep(1:M, length = n)) # needs to have all repeated samples in the same fold
                }
            }
        } else if (validation ==  "loo") {
            folds = split(1:n, rep(1:n, length = n))
            M = n
        }
        
        M = length(folds)
        
        error.sw = matrix(0,nrow = M, ncol = length(test.keepX))
        rownames(error.sw) = paste0("fold",1:M)
        colnames(error.sw) = names(test.keepX)
        # for the last keepX (i) tested, prediction combined for all M folds so as to extract the error rate per class
        # prediction.all = vector(length = nrow(X))
        # in case the test set only includes one sample, it is better to advise the user to
        # perform loocv
        stop.user = FALSE

        # function instead of a loop so we can use lapply and parLapply. Can't manage to put it outside without adding all the arguments
        
        #save(list=ls(),file="temp2.Rdata")
        #result.all=list()
        fonction.j.folds = function(j)#for (j in 1:M)
        {
            if (progressBar ==  TRUE)
            setTxtProgressBar(pb, (M*(nrep-1)+j-1)/(M*nrepeat))
            
            #print(j)
            #---------------------------------------#
            #-- set up leave out samples. ----------#

            omit = which(repeated.measure %in% folds[[j]] == TRUE)
            
            # get training and test set
            X.train = X[-omit, ]
            Y.train.mat = Y[-omit, , drop = FALSE]
            q = ncol(Y.train.mat)
            X.test = X[omit, , drop = FALSE]#matrix(X[omit, ], nrow = length(omit)) #removed to keep the colnames in X.test
            Y.test = Y[omit, , drop = FALSE]
            #-- set up leave out samples. ----------#
            #---------------------------------------#

            #------------------------------------------#
            #-- split the NA in training and testing --#
            if(any(misdata))
            {
                is.na.A.train = lapply(is.na.A, function(x){x[-omit,, drop=FALSE]})
                is.na.A.test = lapply(is.na.A, function(x){x[omit,, drop=FALSE]})
                
                ind.NA.train = lapply(is.na.A.train, function(x){which(apply(x, 1, sum) > 0)}) # calculated only once
                #ind.NA.test = which(apply(is.na.A.test, 1, sum) > 0) # calculated only once

                ind.NA.col.train = lapply(is.na.A.train, function(x){which(apply(x, 2, sum) > 0)}) # calculated only once
                #ind.NA.col.test = which(apply(is.na.A.test, 2, sum) > 0) # calculated only once

            } else {
                is.na.A.train = is.na.A.test =NULL
                ind.NA.train = NULL
                ind.NA.col.train = NULL
            }
            #-- split the NA in training and testing --#
            #------------------------------------------#

            #---------------------------------------#
            #-- near.zero.var ----------------------#
            
            # first remove variables with no variance
            var.train = apply(X.train, 2, var)
            ind.var = which(var.train == 0)
            if (length(ind.var) > 0)
            {
                X.train = X.train[, -c(ind.var),drop = FALSE]
                X.test = X.test[, -c(ind.var),drop = FALSE]
                
                # reduce choice.keepX and test.keepX if needed
                if (any(choice.keepX > ncol(X.train)))
                choice.keepX[which(choice.keepX>ncol(X.train))] = ncol(X.train)

                # reduce test.keepX if needed
                if (any(test.keepX > ncol(X.train)))
                test.keepX[which(test.keepX>ncol(X.train))] = ncol(X.train)
                    
            }
            
            if(near.zero.var == TRUE) # would need to be done for Y as well
            {
                remove.zero = nearZeroVar(X.train)$Position
                
                if (length(remove.zero) > 0)
                {
                    names.var = colnames(X.train)[remove.zero]
                    
                    X.train = X.train[, -c(remove.zero),drop = FALSE]
                    X.test = X.test[, -c(remove.zero),drop = FALSE]
                    
                    # reduce choice.keepX and test.keepX if needed
                    if (any(choice.keepX > ncol(X.train)))
                    choice.keepX[which(choice.keepX>ncol(X.train))] = ncol(X.train)
                    
                    # reduce test.keepX if needed
                    if (any(test.keepX > ncol(X.train)))
                    test.keepX[which(test.keepX>ncol(X.train))] = ncol(X.train)
                    
                    
                }
                #print(remove.zero)
            }
            
            #-- near.zero.var ----------------------#
            #---------------------------------------#
            
            
            prediction.comp.j = array(0, c(length(omit), ncol(Y), length(test.keepX)), dimnames = list(rownames(X.test), colnames(Y), names(test.keepX)))
            
            
            class.comp.j = list()
            class.comp.j[[1]] = matrix(0, nrow = length(omit), ncol = length(test.keepX))# prediction of all samples for each test.keepX and  nrep at comp fixed
            
            # shape input for `internal_mint.block' (keepA, test.keepA, etc)
            result = internal_wrapper.mint(X=X.train, Y=Y.train.mat, study=factor(rep(1,nrow(X.train))), ncomp=ncomp,
            keepX=choice.keepX, keepY=rep(ncol(Y.train.mat), ncomp-1), test.keepX=test.keepX, test.keepY=ncol(Y.train.mat),
            mode="regression", scale=scale, near.zero.var=near.zero.var,
            max.iter=max.iter, logratio="none", DA=TRUE, multilevel=NULL,
            misdata = misdata, is.na.A = is.na.A.train, ind.NA = ind.NA.train,
            ind.NA.col = ind.NA.col.train, all.outputs=FALSE)
            
            # `result' returns loadings and variates for all test.keepX on the ncomp component
            
            # need to find the best keepX/keepY among all the tested models
            
            #---------------------------------------#
            #-- scaling X.test ---------------------#
            
            # we prep the test set for the successive prediction: scale and is.na.newdata
            if (!is.null(attr(result$A[[1]], "scaled:center")))
            X.test = sweep(X.test, 2, STATS = attr(result$A[[1]], "scaled:center"))
            if (scale)
            X.test = sweep(X.test, 2, FUN = "/", STATS = attr(result$A[[1]], "scaled:scale"))
            
            means.Y = matrix(attr(result$A[[2]], "scaled:center"),nrow=nrow(X.test),ncol=q,byrow=TRUE);
            if (scale)
            {sigma.Y = matrix(attr(result$A[[2]], "scaled:scale"),nrow=nrow(X.test),ncol=q,byrow=TRUE)}else{sigma.Y=matrix(1,nrow=nrow(X.test),ncol=q)}
            
            #-- scaling X.test ---------------------#
            #---------------------------------------#


            #-----------------------------------------#
            #-- prediction on X.test for all models --#

            # record prediction results for each test.keepX
            keepA = result$keepA
            test.keepA = keepA[[ncomp]]
            
            for(i in 1:nrow(test.keepA))
            {
                #print(i)
                # creates temporary splsda object to use the predict function
                object.spls.temp = result
                # add the "splsda" class
                class(object.spls.temp) = c("spls")
                
                object.spls.temp$X = result$A$X
                object.spls.temp$Y =  result$A$Y
                
                # only pick the loadings and variates relevant to that test.keepX
                
                names.to.pick = NULL
                if(ncomp>1)
                names.to.pick = unlist(lapply(1:(ncomp-1), function(x){
                    paste(paste0("comp",x),apply(keepA[[x]],1,function(x) paste(x,collapse="_")), sep=":")
                    
                }))
                
                names.to.pick.ncomp = paste(paste0("comp",ncomp),paste(keepA[[ncomp]][i,],collapse="_"), sep=":")
                names.to.pick = c(names.to.pick, names.to.pick.ncomp)


                
                object.spls.temp$variates = lapply(result$variates, function(x){if(ncol(x)!=ncomp) {x[,colnames(x)%in%names.to.pick, drop=FALSE]}else{x}})
                object.spls.temp$loadings = lapply(result$loadings, function(x){if(ncol(x)!=ncomp) {x[,colnames(x)%in%names.to.pick, drop=FALSE]}else{x}})
                
                # added: record selected features
                if (any(class.object == "splsda") & length(test.keepX) ==  1) # only done if splsda and if only one test.keepX as not used if more so far
                # note: if plsda, 'features' includes everything: to optimise computational time, we don't evaluate for plsda object
                features.j = selectVar(object.spls.temp, comp = ncomp)$name
                
                # do the prediction, we are passing to the function some invisible parameters:
                # the scaled newdata and the missing values
                test.predict.sw <- predict(object.spls.temp, newdata.scale = X.test, dist = dist, misdata.all=any(misdata), is.na.X = list(X=is.na.A.train), is.na.newdata = list(X=is.na.A.test))
                prediction.comp.j[, , i] =  test.predict.sw$predict[, , ncomp]
                
                if(ncol(Y) == 1)
                class.comp.j[[1]][,i] = test.predict.sw$predict[, , ncomp]
                
            } # end i
            
            #-- prediction on X.test for all models --#
            #-----------------------------------------#


            return(list(class.comp.j = class.comp.j, prediction.comp.j = prediction.comp.j, features = features.j, omit = omit))
            #result.all[[j]] = list(class.comp.j = class.comp.j,  prediction.comp.j = prediction.comp.j, features = features.j, omit = omit)

        } # end fonction.j.folds
        
            
        if (parallel == TRUE)
        {
            clusterEvalQ(cl, library(mixOmics))
            clusterExport(cl, ls(), envir=environment())
            result.all = parLapply(cl, 1: M, fonction.j.folds)
        } else {
            result.all = lapply(1: M, fonction.j.folds)
            
        }
        #---------------------------#
        #--- combine the results ---#
        
        save(list=ls(),file="temp33.Rdata")

        for(j in 1:M)
        {
            omit = result.all[[j]]$omit
            prediction.comp.j = result.all[[j]]$prediction.comp.j
            class.comp.j = result.all[[j]]$class.comp.j

            prediction.comp[[nrep]][omit, , ] = prediction.comp.j
            
            if(ncol(Y) == 1)
            class.comp[[1]][omit, nrep, ] = class.comp.j[[1]]
            
            if (length(test.keepX) ==  1) # only done if splsda and if only one test.keepX as not used if more so far
            features = c(features, result.all[[j]]$features)

        }
        
        # create a array, for each keepX: n*q*nrep
        for(i in 1:length(test.keepX))
        {
            prediction.keepX[[i]][,,nrep] = prediction.comp[[nrep]][,,i, drop=FALSE]
            
        }
        
        #--- combine the results ---#
        #---------------------------#
        

        if (progressBar ==  TRUE)
        setTxtProgressBar(pb, (M*nrep)/(M*nrepeat))
        
        #---------------------------#
        #----- AUC on the test -----#
        if(auc)
        {
            data=list()
            for (i in 1:length(test.keepX))
            {
                data$outcome = Y
                data$data = prediction.comp[[nrep]][, , i]
                auc.all[[nrep]][, , i] = as.matrix(statauc(data))
            }
        }
        #----- AUC on the test -----#
        #---------------------------#

    } #end nrep 1:nrepeat

    names(prediction.comp) = paste0("nrep.", 1:nrepeat)
    # class.comp[[ijk]] is a matrix containing all prediction for test.keepX, all nrepeat and all distance, at comp fixed
    
    #-------------------------------------------------------------#
    #----- average AUC over the nrepeat, for each test.keepX -----#
    if(auc)
    {
        
        if(nlevels(Y)>2)
        {
            auc.mean.sd =  array(0, c(nlevels(Y),2, length(test.keepX)), dimnames = list(rownames(auc.all[[1]]), c("AUC.mean","AUC.sd"), names(test.keepX)))
        }else{
            auc.mean.sd =  array(0, c(1,2, length(test.keepX)), dimnames = list(rownames(auc.all[[1]]), c("AUC.mean","AUC.sd"), names(test.keepX)))
        }
        
        for(i in 1:length(test.keepX))
        {
            temp = NULL
            for(nrep in 1:nrepeat)
            {
                temp = cbind(temp, auc.all[[nrep]][, 1, i])
            }
            auc.mean.sd[, 1, i] = apply(temp,1,mean)
            auc.mean.sd[, 2, i] = apply(temp,1,sd)
        }
    } else {
        auc.mean.sd = auc.all = NULL
        
    }
    #----- average AUC over the nrepeat, for each test.keepX -----#
    #-------------------------------------------------------------#
    save(list=ls(),file="temp3.Rdata")

    result = list()
    error.mean = error.sd = error.per.class.keepX.opt.comp = keepX.opt = test.keepX.out = mat.error.final = choice.keepX.out = list()

    if (any(measure == "MSE"))
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",names(test.keepX))
        
        #finding the best keepX depending on the error measure: overall or BER
        # MSE error for each nrep and each test.keepX: summing over all samples
        varY=apply(Y,2,var)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
        {
            temp = (Y-x)^2
            temp2=t(t(temp)/varY)
            temp3=apply(temp2,2, sum)
            temp3
        })})
        mat.error.final[[ijk]] = lapply(error, function(x){x/nrow(Y)})  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)

        # we want to average the error per keepX over nrepeat and choose the minimum error
        error.mean[[ijk]] = sapply(error, mean)/nrow(Y) #length number of components
        if (!nrepeat ==  1)
        error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)/nrow(Y)}), mean)
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
        
        result$"MSE"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
        result$"MSE"$error.rate.sd = error.sd
        
        result$"MSE"$mat.error.rate = mat.error.final
        result$"MSE"$keepX.opt = test.keepX.out
    }
    
    if (any(measure == "MAE")) # MAE (Mean Absolute Error: MSE without the square)
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",names(test.keepX))
        
        #finding the best keepX depending on the error measure: overall or BER
        # MSE error for each nrep and each test.keepX: summing over all samples
        varY=apply(Y,2,var)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
            {
                temp = abs(Y-x)
                temp2=t(t(temp)/varY)
                temp3=apply(temp2,2, sum)
                temp3
            })})

        mat.error.final[[ijk]] = lapply(error, function(x){x/nrow(Y)})  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)

        # we want to average the error per keepX over nrepeat and choose the minimum error
        error.mean[[ijk]] = sapply(error, mean)/nrow(Y) #length number of components
        if (!nrepeat ==  1)
        error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)/nrow(Y)}), mean)
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
        
        result$"MAE"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
        result$"MAE"$error.rate.sd = error.sd
        
        result$"MAE"$mat.error.rate = mat.error.final
        result$"MAE"$keepX.opt = test.keepX.out
    }
    
    if (any(measure == "Bias")) # Bias (average of the differences)
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",names(test.keepX))
        
        #finding the best keepX depending on the error measure: overall or BER
        # MSE error for each nrep and each test.keepX: summing over all samples
        varY=apply(Y,2,var)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
            {
                temp = (Y-x)
                temp2=t(t(temp)/varY)
                temp3=abs(apply(temp2,2, sum)) # absolute value of the bias
                temp3
            })})

        mat.error.final[[ijk]] = lapply(error, function(x){x/nrow(Y)})  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)

        # we want to average the error per keepX over nrepeat and choose the minimum error
        error.mean[[ijk]] = sapply(error, mean)/nrow(Y) #length number of components
        if (!nrepeat ==  1)
        error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)/nrow(Y)}), mean)


        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
        
        result$"Bias"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
        result$"Bias"$error.rate.sd = error.sd
        
        result$"Bias"$mat.error.rate = mat.error.final
        result$"Bias"$keepX.opt = test.keepX.out
    }
    
    

 
    result$prediction.comp = prediction.comp
    result$class.comp = class.comp
    return(result)
}
