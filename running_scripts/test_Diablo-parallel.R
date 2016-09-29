#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

## sGCC-DA
# -------------
library(mixOmics)
#source("mixOmics/R/plotIndiv.R")
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- wrapper.sgccda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,#c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE)

set.seed(43)
a=perf(nutrimouse.sgccda)

set.seed(43)
a2=perf(nutrimouse.sgccda,cpus=4)




#test perf diablo
#source("mixOmics/R/tune.diablo.R")

library(parallel)

set.seed(45)
X = data
Y = Y
design = design
ncomp = 2
scheme = "centroid"
verbose = FALSE
bias = FALSE
nrepeat = 11
constraint = FALSE
validation = "Mfold"
folds = 10
dist = "max.dist"
measure = "BER"
progressBar = TRUE
max.iter = 100
near.zero.var = FALSE
scale = TRUE
tol = 1e-06
light.output = TRUE
parallel = TRUE

cpus=2
scheme = "centroid"
mode = "regression"
bias = FALSE
init = "svd"
verbose = FALSE

measure = "Overall.BER"


if (is.null(dim(Y)))
{
    Y = factor(Y)
} else {
    stop("'Y' should be a factor or a class vector.")
}

if (nlevels(Y) == 1)
stop("'Y' should be a factor with more than one level")

#-- already.tested.X
already.tested.X = NULL
test.keepX = lapply(1:length(X),function(x){c(5,10,15)[which(c(5,10,15)<ncol(X[[x]]))]})
names(test.keepX) = names(X)



l = sapply(test.keepX,length)
n = names(test.keepX)
temp = data.frame(l, n)


message(paste("You have provided a sequence of keepX of length: ", paste(apply(temp, 1, function(x) paste(x,collapse=" for block ")), collapse= " and "), ".\nThis results in ",prod(sapply(test.keepX,length)), " models being fitted for each component and each nrepeat, this may take some time to run, be patient!",sep=""))
test.keepX = lapply(test.keepX,sort) #sort test.keepX so as to be sure to chose the smallest in case of several minimum

grid = expand.grid (test.keepX[length(test.keepX):1])[length(test.keepX):1] # each row is to be tested, the reordering is just a personal preference, works without it

# if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
if ((!is.null(already.tested.X)) & length(already.tested.X) > 0)
{
    comp.real = (length(already.tested.X[[1]]) + 1):ncomp
    #check and match already.tested.X to X
    if(constraint == TRUE & length(already.tested.X[[1]]) >0)
    {
        already.tested.X = get.keepA.and.keepA.constraint (X = X, keepX.constraint = already.tested.X, ncomp = rep(length(already.tested.X[[1]]), length(X)))$keepA.constraint
        
        # to get characters in already.tested.X
        already.tested.X = lapply(1:length(already.tested.X), function(y){temp = lapply(1:length(already.tested.X[[y]]), function(x){colnames(X[[y]])[already.tested.X[[y]][[x]]]}); names(temp) = paste("comp", 1:length(already.tested.X[[1]]), sep=""); temp})
        names(already.tested.X) = names(X)
        
    } else if(length(already.tested.X[[1]]) >0)
    {
        if(length(unique(names(already.tested.X)))!=length(already.tested.X) | sum(is.na(match(names(already.tested.X),names(X)))) > 0)
        stop("Each entry of 'already.tested.X' must have a unique name corresponding to a block of 'X'")
        
    }
    
} else {
    comp.real = 1:ncomp
}

if (parallel == TRUE)
{
    cl <- makeCluster(cpus, type = "SOCK")
}

keepX = error.rate = mat.sd.error = NULL
mat.error.rate = list()
error.per.class.keepX.opt=list()


comp = 1
if (progressBar == TRUE)
cat("\ncomp",comp.real[comp], "\n")

#-- set up a progress bar --#
if (progressBar ==  TRUE & comp == 1)
{
    pb = txtProgressBar(style = 3)
    nBar = 1
}

result.comp = matrix(nrow=nrow(grid),ncol=1)
error.mean = error.sd = mat.error.rate.keepX = NULL
error.per.class = array(0,c(nlevels(Y),nrepeat,nrow(grid)),dimnames = list(levels(Y), paste("nrep",1:nrepeat,sep=".")))















design2 = matrix(runif(9),nrow=3,ncol=3)

source("mixOmics/R/tune.diablo.R")

# classic tune
tune = tune.block.splsda(
X = data,
Y = Y,
design = design,
ncomp = 2,#c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE,
nrepeat = 20,
)
tune

system.time(tune.block.splsda(
X = data,
Y = Y,
design = design,
ncomp = 2,#c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE,
nrepeat = 11,
))


source("mixOmics/R/tune.diablo.R")
tune2 = tune.block.splsda(
X = data,
Y = Y,
design = design2,
ncomp = 2,#c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE,
nrepeat = 11,
cpus=2
)

system.time(tune.block.splsda(
X = data,
Y = Y,
design = design,
ncomp = 2,#c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE,
nrepeat = 11,
cpus=10
))


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    cat("no additional tests")










    # classic tune, with test.keepX as input
    tune = tune.block.splsda(
    X = data,
    Y = Y,
    design = design,
    ncomp = 2,#c(2, 2),
    test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),
    scheme = "centroid",
    verbose = FALSE,
    bias = FALSE,
    progressBar=TRUE
    )
    tune

    # tune with constraint
    tune = tune.block.splsda(
    X=data,
    Y = Y,
    design = design,
    constraint=TRUE,
    nrepeat=4,
    ncomp = 2,#c(2, 2),
    scheme = "centroid",
    verbose = FALSE,
    bias = FALSE
    )
    tune
    
    # tune with constraint, with test.keepX as input
    tune = tune.block.splsda(
    X = data,
    Y = Y,
    design = design,
    ncomp = 2,#c(2, 2),
    test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),
    constraint = TRUE,
    scheme = "centroid",
    verbose = FALSE,
    bias = FALSE
    )
    tune

    # tune without constraint, but only component 2
    tune = tune.block.splsda(
    X=data,
    Y = Y,
    design = design,
    constraint = FALSE,
    already.tested.X = list(gene=c(10), lipid=c(15)),
    nrepeat=4,
    ncomp = 2,#c(2, 2),
    scheme = "centroid",
    verbose = FALSE,
    bias = FALSE
    )
    tune


    # tune with constraint, and only component 2 and 3
    tune = tune.block.splsda(
    X=data,
    Y = Y,
    design = design,
    constraint=TRUE,
    already.tested.X = list(gene=list(comp1=sample(colnames(data$gene),10)), lipid=list(comp1=sample(colnames(data$lipid),5))),
    nrepeat=4,
    ncomp = 3,#c(2, 2),
    scheme = "centroid",
    verbose = FALSE,
    bias = FALSE
    )
    tune


}

par(opar)
.





time3 = proc.time()
for(indice.grid in 1 : nrow(grid))
{
    # test.keepX.comp: keepX for each block on component "comp.real[comp]"
    # already.tested.X: either keepX (constraint=FALSE) or keepX.constraint.temp (constraint=TRUE) for all block on components 1:(comp.real[comp]-1)
    # keepX.temp: keepX for all block on all component 1:comp.real[comp], only used if constraint=FALSE
    set.seed(43)
    test.keepX.comp = grid[indice.grid,]
    
    if(constraint == FALSE)
    {
        keepX.temp = lapply(1:length(X), function(x){c(already.tested.X[[x]],test.keepX.comp[[x]])})
        names(keepX.temp) = names(X)
    }
    
    # run block.splsda
    model = suppressMessages(block.splsda(X = X, Y = Y, ncomp=comp.real[comp],
    keepX.constraint = if(constraint){already.tested.X}else{NULL},
    keepX = if(constraint){test.keepX.comp}else{keepX.temp},
    design=design, scheme=scheme, mode=mode, scale=scale,
    bias=bias, init=init, tol=tol, verbose=verbose, max.iter=max.iter, near.zero.var=near.zero.var))
    
    
    # run perf on the model
    cvPerf = lapply(1 : nrepeat, function(u){out = suppressMessages(perf(model, validation = validation, folds = folds, dist = dist));
        if (progressBar ==  TRUE)
        setTxtProgressBar(pb, ((indice.grid-1)*nrepeat+u)/(nrow(grid)*nrepeat))
        out
    })
    
    # record results
    ## Majority Vote
    cvPerf2 = lapply(1 : nrepeat, function(x){unlist(cvPerf[[x]][names(cvPerf[[x]]) == "MajorityClass.error.rate"], recursive = FALSE)})
    mat.error.rate.temp = simplify2array( lapply(cvPerf2, function(x) x$"MajorityClass.error.rate"[measure,comp])) # error over the nrepeat
    mat.error.rate.keepX = rbind(mat.error.rate.keepX, mat.error.rate.temp)
    
    error.mean = c(error.mean, mean(mat.error.rate.temp))
    error.per.class[, , indice.grid] = simplify2array(lapply(cvPerf2, function(x){ x$"MajorityClass.error.rate"[1:nlevels(Y),comp]}))
    
    if(nrepeat > 1)
    error.sd = c(error.sd, apply(simplify2array(lapply(cvPerf2, function(x) x$"MajorityClass.error.rate")), c(1,2), sd)[measure,comp])
    
    if (progressBar ==  TRUE)
    setTxtProgressBar(pb, (indice.grid)/nrow(grid))
}
time4 = proc.time()
time4-time3


