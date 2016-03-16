# Copyright (C) 2015

# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France

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


tune.splsda <- function (X, Y, ncomp = 1, test.keepX = c(5, 10, 15),
already.tested.X = NULL, validation = "Mfold",folds=10,
dist ="max.dist",
progressBar = TRUE,
near.zero.var = FALSE,nrepeat=1,
logratio = c('none','CLR'))
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
    
    if ((!is.null(already.tested.X)) && (length(already.tested.X) != (ncomp - 1)) )
    stop("The number of already tested parameters should be NULL or ", ncomp - 1, " since you set ncomp = ", ncomp)
    
    if ((!is.null(already.tested.X)) && (!is.numeric(already.tested.X)))
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)
    
    if (!is.null(already.tested.X))
    cat("Number of variables selected on the first ", ncomp - 1, "component(s) was ", already.tested.X)
    
    
    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    # -------------------------------------
    # added: first check for near zero var on the whole data set
    
    # and then we start from the X data set with the nzv removed
    
    #-- cross-validation approach  ---------------------------------------------#
    #---------------------------------------------------------------------------#
    
    if((!is.null(already.tested.X)))
    {
        comp.real=ncomp
        ncomp=1
    }else{
        comp.real=1:ncomp
    }
    
    
    dist = match.arg(dist, several.ok = TRUE)
    
    mat.error = matrix(nrow = length(test.keepX), ncol = nrepeat,
    dimnames = list(test.keepX,c(paste('repeat', 1:nrepeat))))
    rownames(mat.error) = test.keepX
    
    mat.error.final = list()
    
    if(!nrepeat==1)
    {
        mat.sd.error = matrix(nrow = length(test.keepX), ncol = ncomp,
        dimnames = list(c(test.keepX), c(paste('comp', comp.real))))
    }else{
        mat.sd.error = matrix(0,nrow = length(test.keepX), ncol = ncomp,
        dimnames = list(c(test.keepX), c(paste('comp', comp.real))))
    }
    
    mat.mean.error = matrix(nrow = length(test.keepX), ncol = ncomp,
    dimnames = list(c(test.keepX), c(paste('comp', comp.real))))
    
    choice.keepX = already.tested.X
    
    error.per.class = list()
    
    prediction.all=list()
    
    final=list()
    error.per.class.mean = matrix(nrow = nlevels(Y), ncol = ncomp,
        dimnames = list(c(levels(Y)), c(paste('comp', comp.real))))
    if(!nrepeat==1)
    {
        error.per.class.sd = matrix(nrow = nlevels(Y), ncol = ncomp,
            dimnames = list(c(levels(Y)), c(paste('comp', comp.real))))
    }else{
        error.per.class.sd = matrix(rep(0,(nlevels(Y)*ncomp)),nrow = nlevels(Y), ncol = ncomp,
            dimnames = list(c(levels(Y)), c(paste('comp', comp.real))))
    }
    #-- set up a progress bar --#
    if (progressBar == TRUE)
    {
        pb = txtProgressBar(style = 3)
        nBar = 1
    }
    
    # successively tune the components until ncomp: comp1, then comp2, ...
    for(comp in 1:ncomp)
    {
        
        for(nrep in 1:nrepeat)
        {
            
            print(nrep)
            
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
            
            if((is.null(already.tested.X)))
            {
                comp.2=comp
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
                print(j)
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
                        X.train = X.train[, -c(remove.zero)]
                        X.test = X.test[, -c(remove.zero)]
                    }
                }
                
                for (i in 1:length(test.keepX))
                {
                    object.res = splsda(X.train, Y.train, ncomp=comp.2, keepX=c(choice.keepX[1:(comp.2-1)], test.keepX[i]),logratio=logratio,near.zero.var=FALSE)
                    test.predict.sw <- predict.mint.splsda(object.res, X.test, method = dist)
                    
                    Prediction.sw <- levels(Y)[test.predict.sw$class[[dist]][, comp]]
                    error.sw[j, i] <- sum(as.character(Y.test) != Prediction.sw)
                    if(j==1)
                    {
                        prediction.all[[i]]=vector(length = nrow(X))
                    }
                    prediction.all[[i]][omit]=Prediction.sw
                    
                    
                } # end i
                
                
                ## for the last keepX (i) tested, prediction combined for all M folds:
                # prediction.all[omit] = Prediction.sw
                
            } # end J (M folds)
            
            #prediction.all[[i]] is the prediction for all folds for test.keepX[i]

            for(len in 1:length(prediction.all))
            {
                if(nrep==1)
                {
                    error.per.class[[len]] = matrix(nrow = nlevels(Y), ncol = nrepeat,
                    dimnames = list(levels(Y),c(paste('repeat', 1:nrepeat))))
                    
                }
                tab.pred =  table(Y, factor(prediction.all[[len]], levels = levels(Y)))
                error.per.class[[len]][ ,nrep] = (rowSums(tab.pred) - diag(tab.pred))/summary(Y)
                
            }
            result <- colSums(error.sw)/length(Y)
            names(result) = paste("var",test.keepX, sep = "")
            mat.error[,nrep]=result
            
            
            if (progressBar == TRUE) {setTxtProgressBar(pb, (nBar/(ncomp*nrepeat)))
                nBar=nBar+1}
            
        } #end nrep 1:nrepeat
        
        
        
        mat.mean.error[, comp] = rowMeans(mat.error)
        
        if(!nrepeat==1)
            mat.sd.error[, comp] = apply(mat.error, 1, sd)
        # get the keepX for the minimal error rate, first keepX value
        which.min = which(mat.mean.error[, comp] == min(mat.mean.error[, comp]))[1]
        mat.error.final[[comp]] = mat.error
        
        error.per.class.mean[, comp] = rowMeans(error.per.class[[which.min]])
        
        if(!nrepeat==1)
            error.per.class.sd[, comp] = apply(error.per.class[[which.min]], 1, sd)
        
        
        final[[comp]]=error.per.class[[which.min]]
        choice.keepX = c(choice.keepX, test.keepX[which.min])
        
    } # end comp
       
    names(mat.error.final)=c(paste('comp', comp.real))
    names(final)=c(paste('comp', comp.real))
    
    if (progressBar == TRUE) cat('\n')
    return(list(
    mat.mean.error = mat.mean.error,
    mat.sd.error=mat.sd.error,
    mat.error.final=mat.error.final,
    choice.keepX =choice.keepX ,
    # error per class for last keepX tested
    error.per.class.mean=error.per.class.mean,
    error.per.class.sd=error.per.class.sd,
    error.per.class=final))
    
}