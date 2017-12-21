#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

## sGCC-DA
# -------------
#library(rARPACK)
#library(matrixStats)
#library(mixOmicsDD)
#source("mixOmics/R/plotIndiv.R")
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, outcome=Y, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- block.splsda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,#c(2, 2),
scheme = "centroid",
tol=1e-30)


head(nutrimouse.sgccda$variates$gene)


if(FALSE)
{source("mixOmicsDD/R/internal_wrapper.mint.R")
source("mixOmicsDD/R/internal_wrapper.mint.block.R")
source("mixOmicsDD/R/internal_mint.block.R")
source("mixOmicsDD/R/internal_mint.block_helpers.R")
source("mixOmicsDD/R/check_entry.R")

source("mixOmicsDD/R/tune.splsda.R")
source("mixOmicsDD/R/tune.block.splsda.R")
source("mixOmicsDD/R/block.splsda.R")
source("mixOmicsDD/R/block.spls.R")
source("mixOmicsDD/R/splsda.R")
source("mixOmicsDD/R/spls.R")
source("mixOmicsDD/R/pls.R")
source("mixOmicsDD/R/plsda.R")
source("mixOmicsDD/R/MCV.block.splsda.R")
source("mixOmicsDD/R/MCV.splsda.R")
source("mixOmicsDD/R/perf.R")

nutrimouse.sgccda <- block.splsda(X=data,
Y = Y,
design = design,
#keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,#c(2, 2),
scheme = "centroid",
tol=1e-30)
}
#source("mixOmicsDD/R/predict.mint.block.pls.R")
#source("mixOmicsDD/R/internal_predict.DA.R")
#source("mixOmicsDD/R/MCVfold.R")




pred=predict(nutrimouse.sgccda, data)
pred2=predict(nutrimouse.sgccda, data,weight=c(1,2))


#nutrimouse.sgccda2 <- wrapper.sgccda(X=c(data[-1],data[1]),
#Y = Y,
#design = design,
#keepX = list(gene=c(10,10), lipid=c(15,15)),
#ncomp = 2,#c(2, 2),
#scheme = "centroid",
#tol=1e-30)

set.seed(43)
a=perf(nutrimouse.sgccda,nrepeat=3, progressBar = FALSE)
plot(a)


set.seed(43)
a=perf(nutrimouse.sgccda, progressBar = FALSE)
plot(a)


set.seed(43)
a2=perf(nutrimouse.sgccda,cpus=4, progressBar = FALSE)


res = all.equal(a[-which(names(a)=="call")],a2[-which(names(a2)=="call")])
if(!isTRUE(res))
stop("problem parallel diablo")


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
plotIndiv(nutrimouse.sgccda,legend=TRUE,ind.names=FALSE, centroid=FALSE)
plotIndiv(nutrimouse.sgccda,legend=TRUE,ind.names=FALSE, centroid=TRUE)
#plotIndiv(nutrimouse.sgccda, blocks = c(1,2), group = nutrimouse$diet,
#ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgccda)

# variable plot on the selected variables
plotVar(nutrimouse.sgccda, col = color.mixo(1:2), cex = c(2,2))


# perf

perf=perf(nutrimouse.sgccda, progressBar = FALSE)
perf$features$stable

perf=perf(nutrimouse.sgccda,validation="loo", progressBar = FALSE)
perf$features$stable


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
   
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
    nrepeat = 11
    )
    tune
    
    plot(tune)

    # classic tune, with test.keepX as input
    set.seed(12)
    tune2 = tune.block.splsda(
    X = data,
    Y = Y,
    design = design,
    ncomp = 2,#c(2, 2),
    test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),
    scheme = "centroid",
    weighted=FALSE
    )
    tune2
    
    
    
    set.seed(12)
    a = tune.block.splsda(
    X = data,
    Y = Y,
    design = design,
    ncomp = 2,#c(2, 2),
    test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),
    scheme = "centroid",
    dist="centroids.dist",
    weighted=TRUE,
    nrepeat=2
    )
    a
    
    set.seed(12)
    a2 = tune.block.splsda(
    X = data,
    Y = Y,
    design = design,
    ncomp = 2,#c(2, 2),
    test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),
    scheme = "centroid",
    dist="centroids.dist",
    weighted=TRUE,
    nrepeat=2,
    cpus=2
    )
    a2
    
    res = all.equal(a[-which(names(a)=="call")],a2[-which(names(a2)=="call")])
    if(!isTRUE(res))
    stop("problem parallel diablo")
    



    # tune
    tune = tune.block.splsda(
    X=data,
    Y = Y,
    design = design,
    nrepeat=4,
    ncomp = 2,#c(2, 2),
    scheme = "centroid"
    )
    tune
    
    # tune with test.keepX as input
    tune = tune.block.splsda(
    X = data,
    Y = Y,
    design = design,
    ncomp = 2,#c(2, 2),
    test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),
    scheme = "centroid"
    )
    tune

    # tune, but only component 2
    tune = tune.block.splsda(
    X=data,
    Y = Y,
    design = design,
    already.tested.X = list(gene=c(10), lipid=c(15)),
    nrepeat=4,
    ncomp = 2,#c(2, 2),
    scheme = "centroid"
    )
    tune


 tune = tune.block.splsda(
    X = data,
    Y = Y,
    design = design,
    ncomp = 2,#c(2, 2),
    scheme = "centroid",
    nrepeat = 2,
    name.save="try",
    cpus=4
    )


}

par(opar)

