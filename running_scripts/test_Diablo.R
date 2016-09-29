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
bias = FALSE,
tol=1e-30)


nutrimouse.sgccda2 <- wrapper.sgccda(X=c(data[-1],data[1]),
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,#c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE,
tol=1e-30)


a=perf(nutrimouse.sgccda)
plot(a)

nutrimouse.sgccda$design

plotIndiv(nutrimouse.sgccda,style="graphics")
plotIndiv(nutrimouse.sgccda,style="graphics",add.legend=FALSE,ind.names=FALSE,abline=TRUE)

plotLoadings(nutrimouse.sgccda)
plotLoadings(nutrimouse.sgccda,contrib="max")

par(mfrow=c(1,1))
circosPlot(nutrimouse.sgccda,cutoff=0.7,ncol=2,size.legend=1.1)
circosPlot(nutrimouse.sgccda,cutoff=0.9,ncol=2,size.legend=1.1)

circosPlot(nutrimouse.sgccda,cutoff=0.5,ncol=2,size.legend=1.1, size.variables = 1.5, comp =2)


circosPlot(nutrimouse.sgccda,cutoff=0.7,ncol=2,size.legend=1.1,size.labels=2)

#compatibility with par(mfrow)
par(mfrow=c(2,2))
circosPlot(nutrimouse.sgccda,cutoff=0.6,ncol=1,size.legend=0.6)
circosPlot(nutrimouse.sgccda,cutoff=0.7,ncol=1,size.legend=0.6)
circosPlot(nutrimouse.sgccda,cutoff=0.8,ncol=1,size.legend=0.6,showIntraLinks=TRUE)
plot(1:10)

plotIndiv(nutrimouse.sgccda,style="graphics",add.legend=FALSE,layout=c(2,2)) # Amrit function
plot(1:10)

par(mfrow=c(2,2))
plotDiablo(nutrimouse.sgccda)


# "something above" changes the plot as the legend of cimDiablo is not where it supposed to, the following command line is basically to open a new plot window
# the problem doesn't appear all the time, so I'm not sure where it comes from
#par(opar)
cimDiablo(nutrimouse.sgccda)
cimDiablo(nutrimouse.sgccda,size.legend=0.3)


plotIndiv(nutrimouse.sgccda)
#plotIndiv(nutrimouse.sgccda, blocks = c(1,2), group = nutrimouse$diet,
#ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgccda)

# variable plot on the selected variables
plotVar(nutrimouse.sgccda, col = color.mixo(1:2), cex = c(2,2))


# perf

perf=perf(nutrimouse.sgccda)
perf$features$stable

perf=perf(nutrimouse.sgccda,validation="loo")
perf$features$stable


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    cat("no additional tests")
    
    #test perf diablo
    #source("mixOmics/R/tune.diablo.R")
    
    set.seed(45)
    # classic tune
    tune = tune.block.splsda(
    X = data,
    Y = Y,
    design = design,
    ncomp = 2,#c(2, 2),
    scheme = "centroid",
    verbose = FALSE,
    bias = FALSE,
    nrepeat = 11
    )
    tune

    # classic tune, with test.keepX as input
    tune = tune.block.splsda(
    X = data,
    Y = Y,
    design = design,
    ncomp = 2,#c(2, 2),
    test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),
    scheme = "centroid",
    verbose = FALSE,
    bias = FALSE
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

