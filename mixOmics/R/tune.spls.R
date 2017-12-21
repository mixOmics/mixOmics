#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#
# created: 2013
# last modified: 05-10-2017
#
# Copyright (C) 2013
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
# tune.spls: chose the optimal number of parameters per component on a spls method, based on MSE
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a numeric matrix
# ncomp: the number of components to include in the model. Default to 1.
# test.keepX: grid of keepX among which to chose the optimal one
# already.tested.X: a vector giving keepX on the components that were already tuned
# validation: Mfold or loo cross validation
# folds: if validation=Mfold, how many folds?
# dist: distance to classify samples. see predict
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# scale: boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# auc: calculate AUC
# progressBar: show progress,
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# nrepeat: number of replication of the Mfold process
# logratio: one of "none", "CLR"
# multilevel: repeated measurement. `multilevel' is passed to multilevel(design = ) in withinVariation. Y is ommited and shouldbe included in `multilevel'
# light.output: if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
# cpus: number of cpus to use. default to no parallel



# I start with no selection on Y. Otherwise we need to be able to tell which submodel of Y is better. I'm afraid a sum(MSE_Yi) is gonna lead to very sparse model (only 1Y), to test.

tune.spls = function (X, Y,
ncomp = 1,
test.keepX = c(5, 10, 15),
already.tested.X,
validation = "Mfold",
folds = 10,
dist = "max.dist",
measure = "MSE", # MAE (Mean Absolute Error: MSE without the square), Bias (average of the differences), MAPE (average of the absolute errors, as a percentage of the actual values)
# can do a R2 per Y (correlation, linear regression R2),
# need to call MSEP (see perf.spls)
scale = TRUE,
progressBar = TRUE,
tol = 1e-06,
max.iter = 100,
near.zero.var = FALSE,
nrepeat = 1,
logratio = c('none','CLR'),
multilevel = NULL,
light.output = TRUE,
cpus
)
{    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #------------------#
    #-- check entries --#
    X = as.matrix(X)
    
    if (length(dim(X)) != 2 || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
    
    Y = as.matrix(Y)

    if (length(dim(Y)) != 2 || !is.numeric(Y))
    stop("'Y' must be a numeric matrix.")
    
    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop("invalid number of variates, 'ncomp'.")
    
    
    #-- validation
    choices = c("Mfold", "loo")
    validation = choices[pmatch(validation, choices)]
    if (is.na(validation))
    stop("'validation' must be either 'Mfold' or 'loo'")
    
    if (validation == 'loo')
    {
        if (nrepeat != 1)
        warning("Leave-One-Out validation does not need to be repeated: 'nrepeat' is set to '1'.")
        nrepeat = 1
    }
    
    
    if (missing(already.tested.X))
    {
        already.tested.X = NULL
    } else {
        if(is.null(already.tested.X) | length(already.tested.X)==0)
        stop("''already.tested.X' must be a vector of keepX values")

        if(is.list(already.tested.X))
        stop("''already.tested.X' must be a vector of keepX values")

        message(paste("Number of variables selected on the first", length(already.tested.X), "component(s):", paste(already.tested.X,collapse = " ")))
    }
    
    if(length(already.tested.X) >= ncomp)
    stop("'ncomp' needs to be higher than the number of components already tuned, which is length(already.tested.X)=",length(already.tested.X) , call. = FALSE)
    
    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- test.keepX
    if (is.null(test.keepX) | length(test.keepX) == 1 | !is.numeric(test.keepX))
    stop("'test.keepX' must be a numeric vector with more than two entries", call. = FALSE)
    
    if(!missing(cpus))
    {
        if(!is.numeric(cpus) | length(cpus)!=1)
        stop("'cpus' must be a numerical value")
        
        parallel = TRUE
        cl = makeCluster(cpus, type = "SOCK")
        clusterExport(cl, c("splsda","selectVar"))
        
        if(progressBar == TRUE)
        message(paste("As code is running in parallel, the progressBar will only show 100% upon completion of each nrepeat/ component.",sep=""))

    } else {
        parallel = FALSE
        cl = NULL
    }
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:ncol(X), sep = "")
        dimnames(X)[[2]] = X.names
    }
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = 1:nrow(X)
        rownames(X)  = ind.names
    }
    
    if (length(unique(rownames(X))) != nrow(X))
    stop("samples should have a unique identifier/rowname")
    if (length(unique(X.names)) != ncol(X))
    stop("Unique indentifier is needed for the columns of X")


    #-- end checking --#
    #------------------#
    
   
    
    #---------------------------------------------------------------------------#
    #-- NA calculation      ----------------------------------------------------#
    
    misdata = c(X=anyNA(X), Y=anyNA(Y)) # Detection of missing data. we assume no missing values in the factor Y
    
    if (any(misdata))
    {
        is.na.A.X = is.na(X)
        is.na.A.Y = is.na(Y)
        is.na.A = list(X=is.na.A.X, Y=is.na.A.Y)
    } else {
        is.na.A = NULL
    }
    #-- NA calculation      ----------------------------------------------------#
    #---------------------------------------------------------------------------#



    test.keepX = sort(test.keepX) #sort test.keepX so as to be sure to chose the smallest in case of several minimum
    names(test.keepX) = test.keepX
    # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
    if ((!is.null(already.tested.X)))
    {
        comp.real = (length(already.tested.X) + 1):ncomp
    } else {
        comp.real = 1:ncomp
    }
    
    
    choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist")
    dist = match.arg(dist, choices, several.ok = TRUE)
    
    
    mat.error.rate = list()
    error.per.class = list()

    mat.sd.error = matrix(0,nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real, sep=''))))
    mat.mean.error = matrix(nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
    dimnames = list(c(test.keepX), c(paste('comp', comp.real, sep=''))))        
   
    # first: near zero var on the whole data set
    if(near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
            X = X[, -nzv$Position, drop=TRUE]
            
            if(ncol(X)==0)
            stop("No more predictors after Near Zero Var has been applied!")
            
        }
    }

    if (is.null(multilevel) | (!is.null(multilevel) && ncol(multilevel) == 2))
    {
        if(light.output == FALSE)
        prediction.all = class.all = list()
        

        
        error.per.class.keepX.opt = list()
        error.per.class.keepX.opt.mean = matrix(0, nrow = nlevels(Y), ncol = length(comp.real),
        dimnames = list(c(levels(Y)), c(paste('comp', comp.real, sep=''))))
        # successively tune the components until ncomp: comp1, then comp2, ...
        for(comp in 1:length(comp.real))
        {

            if (progressBar == TRUE)
            cat("\ncomp",comp.real[comp], "\n")
            
            result = MCVfold.spls (X, Y, multilevel = multilevel, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = 1 + length(already.tested.X),
            choice.keepX = already.tested.X,
            test.keepX = test.keepX, measure = measure, dist = dist, scale=scale,
            near.zero.var = near.zero.var, progressBar = progressBar, tol = tol, max.iter = max.iter,
            cl = cl, parallel = parallel,
            misdata = misdata, is.na.A = is.na.A, class.object="spls")
            
            save(list=ls(),file="temp.Rdata")
            # in the following, there is [[1]] because 'tune' is working with only 1 distance and 'MCVfold.splsda' can work with multiple distances
            mat.error.rate[[comp]] = result[[measure]]$mat.error.rate[[1]]
            mat.mean.error[, comp]=result[[measure]]$error.rate.mean[[1]]
            if (!is.null(result[[measure]]$error.rate.sd[[1]]))
            mat.sd.error[, comp]=result[[measure]]$error.rate.sd[[1]]
            

            # best keepX
            already.tested.X = c(already.tested.X, result[[measure]]$keepX.opt[[1]])
             
            if(light.output == FALSE)
            {
                #prediction of each samples for each fold and each repeat, on each comp
                class.all[[comp]] = result$class.comp[[1]]
                prediction.all[[comp]] = result$prediction.comp
            }
            

            
        } # end comp
        if (parallel == TRUE)
        stopCluster(cl)
        
        names(mat.error.rate) = c(paste('comp', comp.real, sep=''))
        names(already.tested.X) = c(paste('comp', 1:ncomp, sep=''))
        
        if (progressBar == TRUE)
        cat('\n')
        
        # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
        if(nrepeat > 2 & length(comp.real) >1)
        {
            keepX = already.tested.X
            error.keepX = NULL
            for(comp in 1:length(comp.real))
            {
                ind.row = match(keepX[[comp.real[comp]]],test.keepX)
                error.keepX = cbind(error.keepX, apply(matrix(mat.error.rate[[comp]][[ind.row]],nrow=ncol(Y)),2,mean)) # average MSE for all Y, per nrepeat
            }
            colnames(error.keepX) = c(paste('comp', comp.real, sep=''))
            
            opt = t.test.process(error.keepX)
            
            ncomp_opt = comp.real[opt]
        } else {
            ncomp_opt = error.keepX = NULL
        }
        
        result = list(
        error.rate = mat.mean.error,
        error.rate.sd = mat.sd.error,
        error.rate.all = mat.error.rate,
        choice.keepX = already.tested.X,
        choice.ncomp = list(ncomp = ncomp_opt, values = error.keepX))
        
        if(light.output == FALSE)
        {
            names(class.all) = names(prediction.all) = c(paste('comp', comp.real, sep=''))
            result$predict = prediction.all
            result$class = class.all
        }

        result$measure = measure
        result$call = match.call()

        class(result) = "tune.spls"

        return(result)
    } else {
        # if multilevel with 2 factors, we can not do as before because withinvariation depends on the factors, we maximase a correlation
        message("For a two-factor analysis, the tuning criterion is based on the maximisation of the correlation between the components on the whole data set")

        cor.value = vector(length = length(test.keepX))
        names(cor.value) = test.keepX

        for (i in 1:length(test.keepX))
        {
            spls.train = mixOmics::splsda(X, Y, ncomp = ncomp, keepX = c(already.tested.X, test.keepX[i]), logratio = logratio, near.zero.var = FALSE, mode = "regression")
            
            # Note: this is performed on the full data set
            # (could be done with resampling (bootstrap) (option 1) and/or prediction (option 2))
            cor.value[i] = cor(spls.train$variates$X[, ncomp], spls.train$variates$Y[, ncomp])
            #
        }
        return(list(cor.value = cor.value))

    }
}



plot.tune.spls = #plot.spca <- plot.ipca <- plot.sipca <-
function(x, optimal = TRUE, sd = TRUE, legend.position = "topright", col, ...)
{
    
    if (!is.logical(optimal))
    stop("'optimal' must be logical.", call. = FALSE)
    
    
    error <- x$error.rate
    if(sd & !is.null(x$error.rate.sd))
    {
        error.rate.sd = x$error.rate.sd
        ylim = range(c(error + error.rate.sd), c(error - error.rate.sd))
    } else {
        error.rate.sd = NULL
        ylim = range(error)
    }
    
    select.keepX <- x$choice.keepX[colnames(error)]
    comp.tuned = length(select.keepX)
    
    legend=NULL
    measure = x$measure
    
    if (length(select.keepX) < 10)
    {
        #only 10 colors in color.mixo
        if(missing(col))
        col = color.mixo(1:comp.tuned)
    } else {
        #use color.jet
        if(missing(col))
        col = color.jet(comp.tuned)
    }
    if(length(col) != comp.tuned)
    stop("'col' should be a vector of length ", comp.tuned,".")
    
    if(measure == "overall")
    {
        ylab = "Classification error rate"
    } else if (measure == "BER")
    {
        ylab = "Balanced error rate"
    }else if (measure == "MSE"){
        ylab = "MSE"
    }else if (measure == "MAE"){
        ylab = "MAE"
    }else if (measure == "Bias"){
        ylab = "Bias"
    }
    
    #save(list=ls(),file="temp.Rdata")
    
    matplot(rownames(error),error, type = "l", axes = TRUE, lwd = 2, lty = 1, log = "x",
    xlab = "Number of selected features", ylab = ylab,
    col = col, ylim = ylim)
    
    if(optimal)
    {
        for(i in 1:comp.tuned)
        {
            # store coordinates of chosen keepX
            index = which(rownames(error) == select.keepX[i])
            # choseen keepX:
            points(rownames(error)[index], error[index,i], col = col[i], lwd=2, cex=3, pch = 18)
        }
    }
    
    if(!is.null(error.rate.sd))
    {
        for(j in 1:ncol(error))
        plot_error_bar(x = as.numeric(names(error[, j])), y =error[, j] , uiw=error.rate.sd[, j], add=T, col = rep(col[j],each=nrow(error)))#, ...)
    }
    
    
    
    if(length(x$choice.keepX) == 1) #only first comp tuned
    {
        legend = "comp1"
    } else if(length(x$choice.keepX) == comp.tuned) # all components have been tuned
    {
        legend = c("comp1", paste("comp1 to", colnames(error)[-1]))
    } else { #first component was not tuned
        legend = paste("comp1 to", colnames(error))
    }
    
    legend(legend.position, lty = 1, lwd = 2, horiz = FALSE, col = col,
    legend = legend)
    
}


