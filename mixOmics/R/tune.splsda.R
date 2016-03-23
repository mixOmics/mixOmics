#############################################################################################################
# Author :
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2015
# last modified: 15-03-2016
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


confusion_matrix=function(Y.learn,Y.test,pred)
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

get.BER=function(out,nlev)
{
    if(missing(nlev)) nlev=2
    #calculation of the BER
    ClassifResult.temp=out
    diag(ClassifResult.temp)=0
    BER_test=sum(apply(ClassifResult.temp,1,sum,na.rm=TRUE)/apply(out,1,sum,na.rm=TRUE),na.rm=TRUE)/nlev
    return(BER_test)
}



tune.splsda <- function (X, Y, ncomp = 1, test.keepX = c(5, 10, 15),
already.tested.X = NULL, validation = "Mfold",folds=10,
dist ="max.dist",
measure=c("overall","BER"),
progressBar = TRUE,
near.zero.var = FALSE,nrepeat=1,
logratio = c('none','CLR'),
verbose=TRUE)
{    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #------------------#
    #-- check entries --#
    X = as.matrix(X)
    
    if (length(dim(X)) != 2 || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
    
    
    # Testing the input Y
    if (is.null(dim(Y)))
    {
        Y = as.factor(Y)
        ind.mat = unmap(as.numeric(Y))
    }else {
        stop("'Y' should be a factor or a class vector.")
    }
    
    
    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    
    #-- validation
    choices = c("Mfold", "loo")
    validation = choices[pmatch(validation, choices)]
    
    if(validation=='loo')
    {
        nrepeat=1
    }
    
    #if ((!is.null(already.tested.X)) && (length(already.tested.X) != (ncomp - 1)) )
    #stop("The number of already tested parameters should be NULL or ", ncomp - 1, " since you set ncomp = ", ncomp)
    
    if(length(already.tested.X)>=ncomp)
    stop("'ncomp' needs to be higher than the number of components already tuned ('length(already.tested.X)')", call. = FALSE)

    if ((!is.null(already.tested.X)) && (!is.numeric(already.tested.X)))
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)
    
    if (!is.null(already.tested.X))
    cat("Number of variables selected on the first ", length(already.tested.X), "component(s) was ", already.tested.X,"\n")
    


    
    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    # -------------------------------------
    # added: first check for near zero var on the whole data set
    
    # and then we start from the X data set with the nzv removed
    
    #-- cross-validation approach  ---------------------------------------------#
    #---------------------------------------------------------------------------#
    
    
    test.keepX=sort(test.keepX) #sort test.keepX so as to be sure to chose the smallest in case of several minimum
    
    # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
    if((!is.null(already.tested.X)))
    {
        comp.real=(length(already.tested.X)+1):ncomp
    }else{
        comp.real=1:ncomp
    }
    
    
    dist = match.arg(dist, several.ok = TRUE)
    
    mat.error = matrix(nrow = length(test.keepX), ncol = nrepeat,
    dimnames = list(test.keepX,c(paste('repeat', 1:nrepeat))))
    rownames(mat.error) = test.keepX
    
    mat.error.final = list()
    error.per.class = list()
    final=list()
    choice.keepX = already.tested.X

    mat.sd.error = matrix(0,nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real))))
    mat.mean.error = matrix(nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real))))

    error.per.class.mean = matrix(nrow = nlevels(Y), ncol = ncomp-length(already.tested.X),
        dimnames = list(c(levels(Y)), c(paste('comp', comp.real))))
    error.per.class.sd = matrix(0,nrow = nlevels(Y), ncol = ncomp-length(already.tested.X),
        dimnames = list(c(levels(Y)), c(paste('comp', comp.real))))
        
    #-- set up a progress bar --#
    if (progressBar == TRUE)
    {
        pb = txtProgressBar(style = 3)
        nBar = 1
    }
    
    prediction.all=list()
    error.per.class.keepX.opt=list()
    # successively tune the components until ncomp: comp1, then comp2, ...
    for(comp in 1:length(comp.real))
    {
        if(verbose)
        cat("\ncomp",comp.real[comp])
        
        prediction.comp=array(0,c(nrow(X), nrepeat,length(test.keepX)))# prediction of all samples for each test.keepX and  nrep at comp fixed
        for(nrep in 1:nrepeat)
        {
            if(verbose)
            cat("\nnrepeat",nrep,"\n")

            #print(nrep)
            if(near.zero.var==TRUE)
            {
                nzv = nearZeroVar(X)
                if (length(nzv$Position > 0))
                {
                    warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
                    X = X[, -nzv$Position, drop=TRUE]
                    
                    if(ncol(X)==0)
                    {
                        stop("No more predictors after Near Zero Var has been applied!")
                    }
                    
                }
            }
            
            n = nrow(X)
            
            res = list()
            #-- define the folds --#
            if (validation == "Mfold")
            {
                
                if(nrep>1 || comp>1)
                folds=M
                
                if (is.null(folds) || !is.numeric(folds) || folds <
                2 || folds > n)
                {
                    stop("Invalid number of folds.")
                }else{
                    M = round(folds)
                    folds = split(sample(1:n), rep(1:M, length = n))
                }
            }else{
                folds = split(1:n, rep(1:n, length = n))
                M = n
            }

            error.sw = matrix(0,nrow = M, ncol = length(test.keepX))
            rownames(error.sw)=paste0("fold",1:M)
            colnames(error.sw)=test.keepX
            # for the last keepX (i) tested, prediction combined for all M folds so as to extract the error rate per class
            # prediction.all = vector(length = nrow(X))
            # in case the test set only includes one sample, it is better to advise the user to
            # perform loocv
            stop.user = FALSE
            for (j in 1:M)
            {
                #print(j)
                #set up leave out samples.
                omit = folds[[j]]
                
                # see below, we stop the user if there is only one sample drawn on the test set using MFold
                if(length(omit) == 1) stop.user = TRUE
                
                # the training set is NOT scaled
                X.train = X[-omit, ]
                Y.train = Y[-omit]
                X.test = X[omit, ,drop=FALSE]#matrix(X[omit, ], nrow = length(omit)) #removed to keep the colnames in X.test
                Y.test=Y[omit]
                
                #---------------------------------------#
                #-- logratio transformation of X.test --#
                # done on X.learn in the splsda function
                
                transfo=logratio.transfo(X=X.test,logratio=logratio)
                X.test=transfo$X
                
                #-- logratio transformation ------------#
                #---------------------------------------#
                
                if(near.zero.var==TRUE)
                {
                    remove.zero = nearZeroVar(X.train)$Position
                    
                    if (length(remove.zero) > 0)
                    {
                        X.train = X.train[, -c(remove.zero),drop=FALSE]
                        X.test = X.test[, -c(remove.zero),drop=FALSE]
                    }
                }
                
                for (i in 1:length(test.keepX))
                {
                    object.res = splsda(X.train, Y.train, ncomp=comp+length(already.tested.X), keepX=c(choice.keepX[1:(comp+length(already.tested.X)-1)], test.keepX[i]),logratio=logratio,near.zero.var=FALSE,mode="regression")
                    test.predict.sw <- predict(object.res, X.test, method = dist)
                    
                    #Prediction.sw <-
                    #error.sw[j, i] <- sum(as.character(Y.test) != Prediction.sw)
                   
                    prediction.comp[omit,nrep,i]= levels(Y)[test.predict.sw$class[[dist]][, comp.real[comp]]]
                    
                    
                } # end i
                
                
                ## for the last keepX (i) tested, prediction combined for all M folds:
                # prediction.all[omit] = Prediction.sw
                
            } # end j 1:M (M folds)
            
            
            if (progressBar == TRUE) {setTxtProgressBar(pb, (nBar/(length(comp.real)*nrepeat)))
                nBar=nBar+1}
            
        } #end nrep 1:nrepeat
        
        # prediction.comp is a matrix containing all prediction for test.keepX and all nrepeat, at comp fixed
        
        
        
        #finding the best keepX depending on the error measure: overall or BER
        if(measure=="overall")
        {
            # classification error for each nrep and each test.keepX: summing over all samples
            error=apply(prediction.comp,c(3,2),function(x){sum(as.character(Y) != x)})
            rownames(error)=test.keepX
            colnames(error)=paste0("nrep",1:nrepeat)

            # we want to average the error per keepX over nrepeat and choose the minimum error
            error.mean=apply(error,1,mean)/length(Y)
            if(!nrepeat==1)
                error.sd=apply(error,1,sd)/length(Y)
            
            mat.error.final[[comp]] = error/length(Y)  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
            
            keepX.opt=which(error.mean==min(error.mean))[1] # chose the lowest keepX if several minimum
            
            
        }else if(measure=="BER")
        {
            error=apply(prediction.comp,c(3,2),function(x){conf=confusion_matrix(Y.learn=factor(Y),Y.test=factor(Y),pred=x); get.BER(conf,nlev=nlevels(Y))})
            rownames(error)=test.keepX
            colnames(error)=paste0("nrep",1:nrepeat)
            
            # average BER over the nrepeat
            error.mean=apply(error,1,mean)
            if(!nrepeat==1)
                error.sd=apply(error,1,sd)

            mat.error.final[[comp]] = error  # BER for each test.keepX (rows) and each nrepeat (columns)

            keepX.opt=which(error.mean==min(error.mean))[1]
        }
        
        
        mat.sd.error[,comp]=error.sd
        mat.mean.error[,comp]=error.mean
        
        # confusion matrix for keepX.opt
        error.per.class.keepX.opt.comp=apply(prediction.comp[,,keepX.opt,drop=FALSE],2,function(x){conf=confusion_matrix(Y.learn=factor(Y),Y.test=factor(Y),pred=x); out=(apply(conf,1,sum)-diag(conf))/summary(Y)})
        
        rownames(error.per.class.keepX.opt.comp)=levels(Y)
        colnames(error.per.class.keepX.opt.comp)=paste0("nrep",1:nrepeat)

        choice.keepX = c(choice.keepX, test.keepX[keepX.opt])
        
        error.per.class.keepX.opt[[comp]]=error.per.class.keepX.opt.comp
        
    } # end comp
       
    names(mat.error.final)=c(paste('comp', comp.real))
    names(error.per.class.keepX.opt)=c(paste('comp', comp.real))
    
    if (progressBar == TRUE) cat('\n')
    return(list(
    mat.mean.error = mat.mean.error,
    mat.sd.error=mat.sd.error,
    mat.error.final=mat.error.final,
    choice.keepX =choice.keepX ,
    # error per class for last keepX tested
    #error.per.class.mean=error.per.class.mean,
    #error.per.class.sd=error.per.class.sd,
    error.per.class=error.per.class.keepX.opt))
    
}