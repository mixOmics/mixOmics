#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2013
# last modified: 24-08-2016
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


tune.block.splsda = function (
X,
Y,
indY,
ncomp = 2,
test.keepX,
already.tested.X,
constraint = FALSE, #if TRUE, expect a list in already.tested.X, otherwise a number(keepX)
validation = "Mfold",
folds = 10,
dist = "max.dist",
measure = "BER", # one of c("overall","BER")
auc = FALSE,
progressBar = TRUE,
near.zero.var = FALSE,
nrepeat = 1,
max.iter = 500,
design,
scheme,
mode,
scale = TRUE,
bias,
init ,
tol = 1e-06,
verbose,
light.output = TRUE # if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #------------------#
    #-- check entries --#
    
    
    # check inpuy 'Y' and transformation in a dummy matrix
    if(!missing(Y))
    {
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
        stop("'Y' should be a factor with more than one level")
        
    } else if(!missing(indY)) {
        Y = X[[indY]]
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(temp) == 1)
        stop("'X[[indY]]' should be a factor with more than one level")
        
        X = X[-indY] #remove Y from X to pass the arguments simpler to block.splsda
        
    } else if(missing(indY)) {
        stop("Either 'Y' or 'indY' is needed")
        
    }

    
    #-- progressBar
    if (!is.logical(progressBar))
    stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
    
    #-- ncomp
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
        warnings("Leave-One-Out validation does not need to be repeated: 'nrepeat' is set to '1'.")
        nrepeat = 1
    }
    
    
    #-- measure
    measure.input = measure
    if(measure == "BER")
    {
        measure = "Overall.BER"
    } else if (measure == "overall"){
        measure = "Overall.ER"
    } else {
        stop("'measure' must be 'overall' or 'BER'")
    }


    #-- already.tested.X
    if (missing(already.tested.X))
    {
        if(constraint == TRUE)
        {
            already.tested.X = list()
        } else {
            already.tested.X = NULL
        }
    } else {
        if(is.null(already.tested.X) | length(already.tested.X)==0)
        stop("''already.tested.X' must be a vector of keepX values (if 'constraint'= FALSE) or a list (if'constraint'= TRUE) ")
        
        if(constraint == TRUE)
        {
            if(!is.list(already.tested.X))
            stop("''already.tested.X' must be a list since 'constraint' is set to TRUE")
            
            print(paste("A total of",paste(lapply(already.tested.X,length),collapse=" "),"specific variables ('already.tested.X') were selected on the first ", length(already.tested.X), "component(s)"))
        } else {
            if(is.list(already.tested.X))
            stop("''already.tested.X' must be a vector of keepX values since 'constraint' is set to FALSE")
            
            print(paste("Number of variables selected on the first", length(already.tested.X), "component(s):", paste(already.tested.X,collapse = " ")))
        }
    }
    if(length(already.tested.X) >= ncomp)
    stop("'ncomp' needs to be higher than the number of components already tuned, which is length(already.tested.X)=",length(already.tested.X) , call. = FALSE)
    
    if (any(is.na(validation)) || length(validation) > 1)
    stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- test.keepX
    if(missing(test.keepX))
    {
        test.keepX = lapply(1:length(X),function(x){c(5,10,15)[which(c(5,10,15)<ncol(X[[x]]))]})
        names(test.keepX) = names(X)
        
    } else {
        if(length(test.keepX) != length(X))
        stop(paste("test.keepX should be a list of length ", length(X),", corresponding to the blocks: ", paste(names(X),collapse=", "), sep=""))
    }
    
    message(paste("You have provided a sequence of keepX of length: ", paste(lapply(test.keepX,length),collapse=", "), " for the blocks: ",paste(names(test.keepX),collapse=", "),", respectively.\nThis results in ",prod(sapply(test.keepX,length)), " models being fitted, it may take some time...",sep=""))
    
    
    #-- end checking --#
    #------------------#
    
    test.keepX = lapply(test.keepX,sort) #sort test.keepX so as to be sure to chose the smallest in case of several minimum

    grid = expand.grid (test.keepX[length(test.keepX):1])[length(test.keepX):1] # each row is to be tested, the reordering is just a personal preference, works without it


    # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
    if ((!is.null(already.tested.X)))
    {
        comp.real = (length(already.tested.X) + 1):ncomp
        #check and match already.tested.X to X
        if(constraint == TRUE & length(already.tested.X) >0)
        {
            already.tested.X = get.keepA.and.keepA.constraint (X = X, keepX.constraint = already.tested.X, ncomp = length(already.tested.X))$keepA.constraint$X
            #transform already.tested.X to characters
            already.tested.X = relist(colnames(X), skeleton = already.tested.X)
        }
        
    } else {
        comp.real = 1:ncomp
    }


    keepX = error.rate = mat.sd.error = NULL
    mat.error.rate = list()
    error.per.class.keepX.opt=list()
    # successively tune the components until ncomp: comp1, then comp2, ...
    for(comp in 1:length(comp.real))
    {
        if (progressBar == TRUE)
        cat("\ncomp",comp.real[comp], "\n")
        
        result.comp = matrix(nrow=nrow(grid),ncol=1)
        error.mean = error.sd = mat.error.rate.keepX = NULL
        error.per.class = array(0,c(nlevels(Y),nrepeat,nrow(grid)),dimnames = list(levels(Y), paste("nrep",1:nrepeat,sep=".")))
        for(indice.grid in 1 : nrow(grid))
        {
            # keepX to be tested for each block on component "comp"
            test.keepX.comp = grid[indice.grid,]
            
            #print(paste("keepX:",paste(keepX.temp.comp, collapse = " ")))
            
            
            # test.keepX.comp: keepX for each block on component "comp.real[comp]"
            # already.tested.X: either keepX (constraint=FALSE) or keepX.constraint.temp (constraint=TRUE) for all block on all component 1:(comp.real[comp]-1)
            # keepX.temp: keepX for all block on all component 1:comp.real[comp], only used if constraint=FALSE
            
            if(constraint == FALSE)
            {
                keepX.temp = lapply(1:length(X), function(x){c(already.tested.X[[x]],test.keepX.comp[[x]])})
                names(keepX.temp) = names(X)

            }

            #print(keepX.temp)
            
            
            # run block.splsda
            model = block.splsda(X = X, Y = Y, ncomp=comp,
            keepX.constraint = if(constraint){already.tested.X}else{NULL},
            keepX = if(constraint){test.keepX.comp}else{keepX.temp})

#           model = block.splsda(X = X, Y = Y, ncomp=ncomp,
#           keepX.constraint=keepX.constraint.temp, keepX=keepX.temp, design=design, scheme=scheme, mode=mode, scale=scale,
#           bias=bias, init=init, tol=tol, verbose=verbose, max.iter=max.iter, near.zero.var=near.zero.var)
            
            # run perf on the model
            cvPerf = lapply(1 : nrepeat, function(u){perf(model, validation = validation, folds = folds, dist = dist)})

            # record results
            ## Majority Vote
            cvPerf2 = lapply(1 : nrepeat, function(x){unlist(cvPerf[[x]][names(cvPerf[[x]]) == "MajorityClass.error.rate"], recursive = FALSE)})
            mat.error.rate.temp = simplify2array( lapply(cvPerf2, function(x) x$"MajorityClass.error.rate"[measure,comp])) # error over the nrepeat
            mat.error.rate.keepX = rbind(mat.error.rate.keepX, mat.error.rate.temp)
            
            error.mean = c(error.mean, mean(mat.error.rate.temp))#apply(simplify2array(lapply(cvPerf2, function(x) x$"MajorityClass.error.rate")), c(1,2), mean)[measure,comp]
            error.per.class[,,indice.grid] = simplify2array(lapply(cvPerf2, function(x){ x$"MajorityClass.error.rate"[1:nlevels(Y),comp]}))
            
            if(nrepeat > 1)
            error.sd = c(error.sd, apply(simplify2array(lapply(cvPerf2, function(x) x$"MajorityClass.error.rate")), c(1,2), sd)[measure,comp])
            
            
                #-- set up a progress bar --#
                if (progressBar ==  TRUE & comp == 1)
                {
                    pb = txtProgressBar(style = 3)
                    nBar = 1
                }
            
            if (progressBar ==  TRUE)
            setTxtProgressBar(pb, (indice.grid)/nrow(grid))

        }
        
        names(error.mean) = apply(grid,1,function(x){paste(x, collapse = "_")})
        if(nrepeat > 1)
        names(error.sd) = names(error.mean)
        
        min.keepX = names(which.min(error.mean))
        a=as.numeric(strsplit(min.keepX,"_")[[1]]) # vector of each optimal keepX for all block on component comp.real[comp]
        
        # best keepX
        opt.keepX.comp = as.list(a)
        names(opt.keepX.comp) = names(X)
        
        error.per.class.keepX.opt[[comp]] = error.per.class[, , which.min(error.mean)]
        error.rate = cbind(error.rate, error.mean)
        mat.sd.error = cbind(mat.sd.error, error.sd)
        mat.error.rate [[comp]] = mat.error.rate.keepX
        
        if(!constraint)
        {
            already.tested.X = lapply(1:length(X), function(x){c(already.tested.X[[x]],opt.keepX.comp[[x]])})

        } else {
            fit = block.splsda(X = X, Y = Y, ncomp=comp,
            keepX.constraint = already.tested.X,
            keepX = opt.keepX.comp)
            
            varselect = selectVar(fit, comp = comp.real[comp])
            varselect = varselect[which(names(varselect) %in% names(X))]
            
            save(list=ls(),file = "temp.Rdata")
            if(length(already.tested.X) == 0)
            {
                already.tested.X = lapply(1:length(X), function(x){already.tested.X[[x]] = list()})
                already.tested.X = lapply(1:length(X),function(x){already.tested.X[[x]][[comp.real[comp]]] = list("comp1" = selectVar(fit, comp = 1)[[x]]$"name")})

            } else {
                already.tested.X = lapply(1:length(X),function(x){already.tested.X[[x]] = c(already.tested.X[[x]], list(selectVar(fit, comp = comp.real[comp])[[x]]$"name")); names(already.tested.X[[x]]) = paste("comp",1:comp.real[comp],sep=""); return(already.tested.X[[x]])})

                #paste("comp",comp.real[comp],sep="") =
            # already.tested needs to be [dataset][comp]
            
            }
        }
        names(already.tested.X) = names(X)


#print(already.tested.X)
    }
    cat("\n")
   
    #    print(already.tested.X)
    colnames(error.rate) = paste("comp", comp.real)
    names(mat.error.rate) = c(paste('comp', comp.real))
    mat.error.rate = lapply(mat.error.rate, function(x) {colnames(x) = paste("nrep",1:nrepeat,sep="."); rownames(x) = rownames(error.rate);x})
    
    names(error.per.class.keepX.opt) = c(paste('comp', comp.real))

    result = list(
    error.rate = error.rate,
    error.rate.sd = mat.sd.error,
    error.rate.all = mat.error.rate,
    choice.keepX = if(constraint){lapply(already.tested.X, function(x){sapply(x,length)})}else{already.tested.X},
    choice.keepX.constraint = if(constraint){already.tested.X}else{NULL},
    error.rate.class = error.per.class.keepX.opt)
    
    result$measure = measure.input
    result$call = match.call()

    class(result) = "tune.block.splsda"
    
    return(result)
    
}

