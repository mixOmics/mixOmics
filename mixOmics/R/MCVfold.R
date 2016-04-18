#############################################################################################################
# Author :
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2015
# last modified: 11-04-2016
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
# ncomp: the number of components to include in the model. Default to 1.
# test.keepX: grid of keepX among which to chose the optimal one
# already.tested.X: a vector giving keepX on the components that were already tuned
# validation: Mfold or loo cross validation
# folds: if validation=Mfold, how many folds?
# dist: distance to classify samples. see predict
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# progressBar: show progress,
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# nrepeat: number of replication of the Mfold process
# logratio = c('none','CLR'). see splsda
# verbose: if TRUE, shows component and nrepeat being tested.


get.confusion_matrix=function(Y.learn,Y.test,pred)
{
    ClassifResult=array(0,c(nlevels(factor(Y.learn)),nlevels(factor(Y.learn))))
    rownames(ClassifResult)=levels(factor(Y.learn))
    colnames(ClassifResult)=paste("predicted.as.",levels(factor(Y.learn)),sep="")
    #--------record of the classification accuracy for each level of Y
    for(i in 1:nlevels(factor(Y.learn)))
    {
        ind.i=which(Y.test==levels(factor(Y.learn))[i])
        for(ij in 1:nlevels(factor(Y.learn)))
        {
            ClassifResult[i,ij]=sum(pred[ind.i]==levels(Y.learn)[ij])
            
        }
    }
    ClassifResult
}

get.BER=function(X)
{
    if(!is.numeric(X)| !is.matrix(X) | length(dim(X)) != 2 | nrow(X)!=ncol(X))
    stop("'X' must be a square numeric matrix")
    
    nlev=nrow(X)
    #calculation of the BER
    ClassifResult.temp=X
    diag(ClassifResult.temp)=0
    BER=sum(apply(ClassifResult.temp,1,sum,na.rm=TRUE)/apply(X,1,sum,na.rm=TRUE),na.rm=TRUE)/nlev
    return(BER)
}



MCVfold.splsda = function(
X,
Y,
validation,
folds,
nrepeat = 1,
ncomp,
choice.keepX,
test.keepX,
measure = c("overall"), # one of c("overall","BER")
dist = "max.dist",
logratio = c('none','CLR'),
near.zero.var = FALSE,
progressBar = TRUE,
class.object = NULL
)
{    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    #-- set up a progress bar --#
    if (progressBar == TRUE)
    {
        pb = txtProgressBar(style = 3)
        nBar = 1
    } else {
        pb = FALSE
    }
    
    M = length(folds)
    features = NULL
    prediction.comp = list()
    for(ijk in dist)
    prediction.comp[[ijk]] = array(0, c(nrow(X), nrepeat, length(test.keepX)))# prediction of all samples for each test.keepX and  nrep at comp fixed
    for(nrep in 1:nrepeat)
    {
        n = nrow(X)
        #-- define the folds --#
        if (validation == "Mfold")
        {
            
            if (nrep > 1) # reinitialise the folds
            folds = M
            
            if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
            {
                stop("Invalid number of folds.")
            } else {
                M = round(folds)
                folds = split(sample(1:n), rep(1:M, length = n))
            }
        }else{
            folds = split(1:n, rep(1:n, length = n))
            M = n
        }
        
        
        error.sw = matrix(0,nrow = M, ncol = length(test.keepX))
        rownames(error.sw) = paste0("fold",1:M)
        colnames(error.sw) = test.keepX
        # for the last keepX (i) tested, prediction combined for all M folds so as to extract the error rate per class
        # prediction.all = vector(length = nrow(X))
        # in case the test set only includes one sample, it is better to advise the user to
        # perform loocv
        stop.user = FALSE
        for (j in 1:M)
        {
            if (progressBar == TRUE)
            setTxtProgressBar(pb, (M*(nrep-1)+j)/(M*nrepeat))
            
            #print(j)
            #set up leave out samples.
            omit = folds[[j]]
            
            # see below, we stop the user if there is only one sample drawn on the test set using MFold
            if(length(omit) == 1)
            stop.user = TRUE
            
            # get training and test set
            X.train = X[-omit, ]
            Y.train = Y[-omit]
            X.test = X[omit, , drop = FALSE]#matrix(X[omit, ], nrow = length(omit)) #removed to keep the colnames in X.test
            Y.test = Y[omit]
            
            #---------------------------------------#
            #-- logratio transformation of X.test --#
            # done on X.learn in the splsda function
            
            transfo=logratio.transfo(X=X.test,logratio=logratio)
            X.test=transfo$X
            
            #-- logratio transformation ------------#
            #---------------------------------------#
            
            #---------------------------------------#
            #-- near.zero.var ----------------------#
            if(near.zero.var==TRUE)
            {
                remove.zero = nearZeroVar(X.train)$Position
                
                if (length(remove.zero) > 0)
                {
                    X.train = X.train[, -c(remove.zero),drop=FALSE]
                    X.test = X.test[, -c(remove.zero),drop=FALSE]
                }
            }
            #-- near.zero.var ----------------------#
            #---------------------------------------#
            
            for (i in 1:length(test.keepX))
            {
                object.res = splsda(X.train, Y.train, ncomp = ncomp, keepX = c(choice.keepX, test.keepX[i]), logratio = logratio, near.zero.var = FALSE, mode = "regression")
                
                # added: record selected features
                if (any(class.object %in% c("splsda")) & length(test.keepX) == 1) # only done if splsda and if only one test.keepX as not used if more so far
                # note: if plsda, 'features' includes everything: to optimise computational time, we don't evaluate for plsda object
                features = c(features, selectVar(object.res, comp = ncomp)$name)
                
                test.predict.sw <- predict(object.res, newdata=X.test, method = dist)
                save(list=ls(),file="temp.Rdata")
                
                for(ijk in dist)
                prediction.comp[[ijk]][omit,nrep,i]= levels(Y)[test.predict.sw$class[[ijk]][, ncomp]]
                
            } # end i
            
        } # end j 1:M (M folds)
        
        
    } #end nrep 1:nrepeat
    
    # prediction.comp[[ijk]] is a matrix containing all prediction for test.keepX, all nrepeat and all distance, at comp fixed
    
    
    result=list()
    
    
    error.mean = error.sd = error.per.class.keepX.opt.comp = keepX.opt = test.keepX.out = mat.error.final = choice.keepX = list()
    
    if (any(measure=="overall"))
    {
        for(ijk in dist)
        {
            rownames(prediction.comp[[ijk]]) = rownames(X)
            colnames(prediction.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            dimnames(prediction.comp[[ijk]])[[3]] = paste0("test.keepX.",test.keepX)
            
            #finding the best keepX depending on the error measure: overall or BER
            # classification error for each nrep and each test.keepX: summing over all samples
            error = apply(prediction.comp[[ijk]],c(3,2),function(x)
            {
                sum(as.character(Y) != x)
            })
            rownames(error) = test.keepX
            colnames(error) = paste0("nrep.",1:nrepeat)
            
            # we want to average the error per keepX over nrepeat and choose the minimum error
            error.mean[[ijk]] = apply(error,1,mean)/length(Y)
            if (!nrepeat == 1)
            error.sd[[ijk]] = apply(error,1,sd)/length(Y)
            
            mat.error.final[[ijk]] = error/length(Y)  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
            
            keepX.opt[[ijk]] = which(error.mean[[ijk]] == min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt.comp[[ijk]] = apply(prediction.comp[[ijk]][, , keepX.opt[[ijk]], drop = FALSE], 2, function(x)
            {
                conf = get.confusion_matrix(Y.learn = factor(Y), Y.test = factor(Y), pred = x)
                out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
            })
            
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            colnames(error.per.class.keepX.opt.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            
            
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            choice.keepX[[ijk]] = c(choice.keepX[[ijk]], test.keepX.out)
            
            
            result$"overall"$error.rate.mean = error.mean
            if (!nrepeat == 1)
            result$"overall"$error.rate.sd = error.sd
            
            result$"overall"$confusion = error.per.class.keepX.opt.comp
            result$"overall"$mat.error.rate = mat.error.final
            result$"overall"$keepX.opt = test.keepX.out
        }
    }
    
    if (any(measure == "BER"))
    {
        for(ijk in dist)
        {
            rownames(prediction.comp[[ijk]]) = rownames(X)
            colnames(prediction.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            dimnames(prediction.comp[[ijk]])[[3]] = paste0("test.keepX.",test.keepX)
            
            error = apply(prediction.comp[[ijk]],c(3,2),function(x)
            {
                conf=get.confusion_matrix(Y.learn=factor(Y),Y.test=factor(Y),pred=x)
                get.BER(conf)
            })
            rownames(error) = test.keepX
            colnames(error) = paste0("nrep.",1:nrepeat)
            
            # average BER over the nrepeat
            error.mean[[ijk]] = apply(error,1,mean)
            if (!nrepeat == 1)
            error.sd[[ijk]] = apply(error,1,sd)
            
            mat.error.final[[ijk]] = error  # BER for each test.keepX (rows) and each nrepeat (columns)
            
            keepX.opt[[ijk]] = which(error.mean[[ijk]] == min(error.mean[[ijk]]))[1]
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt.comp[[ijk]] = apply(prediction.comp[[ijk]][, , keepX.opt[[ijk]], drop = FALSE], 2, function(x)
            {
                conf = get.confusion_matrix(Y.learn = factor(Y), Y.test = factor(Y), pred = x)
                out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
            })
            
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            colnames(error.per.class.keepX.opt.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
            
            
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            choice.keepX[[ijk]] = c(choice.keepX[[ijk]], test.keepX.out)
            
            result$"BER"$error.rate.mean = error.mean
            if (!nrepeat == 1)
            result$"BER"$error.rate.sd = error.sd
            
            result$"BER"$confusion = error.per.class.keepX.opt.comp
            result$"BER"$mat.error.rate = mat.error.final
            result$"BER"$keepX.opt = test.keepX.out
            
        }
        
        
    }
    
    
    result$prediction.comp = prediction.comp
    result$features$stable = sort(table(as.factor(features))/M/nrepeat, decreasing = TRUE)
    return(result)
}



