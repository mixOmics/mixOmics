#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

## sGCC-DA
# -------------
#library(mixOmicsv6)
#source("mixOmics/R/plotIndiv.R")
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- wrapper.sgccda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = c(2, 3),
scheme = "centroid",
verbose = FALSE,
bias = FALSE)



plotIndiv(nutrimouse.sgccda,style="graphics")
plotIndiv(nutrimouse.sgccda,style="graphics",add.legend=FALSE,ind.names=FALSE,abline.line=TRUE)

plotLoadings(nutrimouse.sgccda)
plotLoadings(nutrimouse.sgccda,contrib="max")

par(mfrow=c(1,1))
circosPlot(nutrimouse.sgccda,corThreshold=0.7,ncol=2,cex.legend=1.1)

circosPlot(nutrimouse.sgccda,corThreshold=0.7,ncol=2,cex.legend=1.1,cex.labels=2)

#compatibility with par(mfrow)
par(mfrow=c(2,2))
circosPlot(nutrimouse.sgccda,corThreshold=0.6,ncol=1,cex.legend=0.6)
circosPlot(nutrimouse.sgccda,corThreshold=0.7,ncol=1,cex.legend=0.6)
circosPlot(nutrimouse.sgccda,corThreshold=0.8,ncol=1,cex.legend=0.6,showIntraLinks=TRUE)
plot(1:10)

plotIndiv(nutrimouse.sgccda,style="graphics",add.legend=FALSE,layout=c(2,2)) # Amrit function
plot(1:10)

par(mfrow=c(2,2))
plotDiablo(nutrimouse.sgccda)


# "something above" changes the plot as the legend of cimDiablo is not where it supposed to, the following command line is basically to open a new plot window
# the problem doesn't appear all the time, so I'm not sure where it comes from
#par(opar)
cimDiablo(nutrimouse.sgccda)
cimDiablo(nutrimouse.sgccda,pos.legend="right")
cimDiablo(nutrimouse.sgccda,pos.legend="right",cex.legend=0.3)


plotIndiv(nutrimouse.sgccda)
#plotIndiv(nutrimouse.sgccda, blocks = c(1,2), group = nutrimouse$diet,
#plot.ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgccda)

# variable plot on the selected variables
plotVar(nutrimouse.sgccda, col = color.mixo(1:2), cex = c(2,2))


# perf

perf=perf(nutrimouse.sgccda)
perf$features$stable

perf=perf(nutrimouse.sgccda,validation="loo",parallel=TRUE,cpus=4)
perf$features$stable


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    cat("no additional tests")
    
}

par(opar)

