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


tune.splsda <- function (X, Y, ncomp = 1, test.keepX = c(5, 10, 15),
already.tested.X = NULL, validation = "Mfold",folds=10,
dist ="max.dist",
measure=c("overall"), # one of c("overall","BER")
progressBar = TRUE,
near.zero.var = FALSE,nrepeat=1,
logratio = c('none','CLR')
)
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
    
    mat.error.rate = list()
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
        

    
    # first check for near zero var on the whole data set
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
    

    prediction.all=list()
    error.per.class.keepX.opt=list()
    # successively tune the components until ncomp: comp1, then comp2, ...
    for(comp in 1:length(comp.real))
    {

        if (progressBar == TRUE)
        cat("\ncomp",comp.real[comp], "\n")
        

        result = MCVfold.splsda (X, Y, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = comp + length(choice.keepX), choice.keepX = choice.keepX,
        test.keepX = test.keepX, measure = measure, dist = dist, logratio = logratio , near.zero.var = near.zero.var, progressBar = progressBar, class.object = "splsda")
        
        save(list=ls(),file="temp.Rdata")
        # in the following, there is [[1]] because 'tune' is working with only 1 distance and 'MCVfold.splsda' can work with multiple distances
        
        if (!is.null(result[[measure]]$error.rate.sd[[1]]))
        mat.sd.error[, comp]=result[[measure]]$error.rate.sd[[1]]
        mat.mean.error[, comp]=result[[measure]]$error.rate.mean[[1]]
        
        # confusion matrix for keepX.opt
        error.per.class.keepX.opt[[comp]]=result[[measure]]$confusion[[1]]

        # best keepX
        choice.keepX = c(choice.keepX, result[[measure]]$keepX.opt[[1]])
        
        mat.error.rate[[comp]] = result[[measure]]$mat.error.rate[[1]]
        
        #prediction of each samples for each fold and each repeat, on each comp
        prediction.all[[comp]] = result$prediction.comp[[1]]
        
    } # end comp

    names(mat.error.rate)=c(paste('comp', comp.real))
    names(error.per.class.keepX.opt)=c(paste('comp', comp.real))
    names(prediction.all)=c(paste('comp', comp.real))
    
    if (progressBar == TRUE)
    cat('\n')
    
    return(list(
    mat.mean.error = mat.mean.error,
    mat.sd.error = mat.sd.error,
    mat.error.rate = mat.error.rate,
    choice.keepX = choice.keepX ,
    # error per class for last keepX tested
    #error.per.class.mean=error.per.class.mean,
    #error.per.class.sd=error.per.class.sd,
    error.per.class = error.per.class.keepX.opt,
    prediction.all = prediction.all))
    
}