#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

toxicity.spls = spls(X, Y, ncomp = 2, keepX = c(50, 50),
keepY = c(10, 10))

plotLoadings(toxicity.spls)

# with xlim
xlim = matrix(c(-0.1,0.3, -0.4,0.6), nrow=2, byrow=T)
plotLoadings(toxicity.spls, xlim = xlim)


data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))

# contribution on comp 1, based on the median.
# Colors indicate the group in which the median expression is maximal
plotLoadings(splsda.liver, comp = 1, method = 'median')
plotLoadings(splsda.liver, comp = 1, method = 'median', contrib = "max")

# contribution on comp 2, based on median.
#Colors indicate the group in which the median expression is maximal
plotLoadings(splsda.liver, comp = 2, method = 'median', contrib = "max")

# contribution on comp 2, based on median.
# Colors indicate the group in which the median expression is minimal
plotLoadings(splsda.liver, comp = 2, method = 'median', contrib = 'min')

# changing the name to gene names
# if the user input a name.var but names(name.var) is NULL,
# then a warning will be output and assign names of name.var to colnames(X)
# this is to make sure we can match the name of the selected variables to the contribution plot.
name.var = liver.toxicity$gene.ID[, 'geneBank']
length(name.var)
plotLoadings(splsda.liver, comp = 2, method = 'median', name.var = name.var, title = "Liver data", contrib = "max")

# if names are provided: ok, even when NAs
name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)
plotLoadings(splsda.liver, comp = 2, method = 'median',
name.var = name.var, cex.name = 0.5, contrib = "max")

#missing names of some genes? complete with the original names
plotLoadings(splsda.liver, comp = 2, method = 'median',
name.var = name.var, cex.name = 0.5,complete.name.var=TRUE, contrib = "max")

# look at the contribution (median) for each variable
plot.contrib = plotLoadings(splsda.liver, comp = 2, method = 'median', plot = FALSE, contrib = "max")
head(plot.contrib)
# change the title of the legend and title name
plotLoadings(splsda.liver, comp = 2, method = 'median', legend.title = 'Time',
title = 'Contribution plot', contrib = "max")

# no legend
plotLoadings(splsda.liver, comp = 2, method = 'median', legend = FALSE, contrib = "max")

# change the color of the legend
plotLoadings(splsda.liver, comp = 2, method = 'median', legend.color = c(1:4), contrib = "max")


data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample,
stimul = vac18$stimulation)

# multilevel sPLS-DA model
res.1level = splsda(X, Y = design[,2], ncomp = 3, multilevel = design[,1],
 keepX = c(30, 137, 123))


name.var = vac18$tab.prob.gene[, 'Gene']
names(name.var) = colnames(X)

plotLoadings(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
name.var = name.var, cex.name = 0.2, contrib = "max")

# too many transcripts? only output the top ones
plotLoadings(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
name.var = name.var, cex.name = 0.5, ndisplay = 60, contrib = "max")

# breast tumors
# ---
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2,near.zero.var=TRUE)

name.var = as.character(breast.tumors$genes$name)
names(name.var) = colnames(X)

# with gene IDs, showing the top 60
plotLoadings(plsda.breast, contrib = 'max', comp = 1, method = 'median',
ndisplay = 60,
name.var = name.var,
cex.name = 0.6,
legend.color = color.mixo(1:2))

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$treatment[, 4]

plsda.liver <- plsda(X, Y, ncomp = 2)
plotIndiv(plsda.liver, ind.names = Y, ellipse = TRUE)


name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)

plotLoadings(plsda.liver, contrib = 'max', comp = 1, method = 'median', ndisplay = 100,
name.var = name.var, cex.name = 0.4,
legend.color = color.mixo(1:4))


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################
### splsda, plsda

if(additional.test==TRUE)
{
    data(liver.toxicity)
    X <- as.matrix(liver.toxicity$gene)
    Y <- as.factor(liver.toxicity$treatment[, 4])
    
    splsda.liver <- splsda(X, Y, ncomp = 3, keepX = c(20, 20, 20))
    
    #test contrib
    plotLoadings(splsda.liver,contrib='max')
    plotLoadings(splsda.liver,contrib='min')
    
    #test method
    plotLoadings(splsda.liver,method='mean', contrib = "max")
    plotLoadings(splsda.liver,method='median', contrib = "max")
    
    #test comp
    plotLoadings(splsda.liver,method='mean',comp=1, contrib = "max")
    plotLoadings(splsda.liver,method='mean',comp=2, contrib = "max")
    plotLoadings(splsda.liver,method='mean',comp=3, contrib = "max")
    plotLoadings(splsda.liver,method='median',comp=1, contrib = "max")
    plotLoadings(splsda.liver,method='median',comp=2, contrib = "max")
    plotLoadings(splsda.liver,method='median',comp=3, contrib = "max")
    
    #test ties
    plotLoadings(splsda.liver,show.ties=FALSE, contrib = "max")
    plotLoadings(splsda.liver,show.ties=TRUE, contrib = "max")
    
    #test cex.name
    plotLoadings(splsda.liver,cex.name=0.5, contrib = "max")
    plotLoadings(splsda.liver,cex.name=1.5, contrib = "max")
    
    #test cex.legend
    plotLoadings(splsda.liver,cex.legend=0.5, contrib = "max")
    plotLoadings(splsda.liver,cex.legend=1, contrib = "max")
    
    #test name.var
    plotLoadings(splsda.liver,name.var=liver.toxicity$gene.ID[, 'gene.title'], contrib = "max")
    plotLoadings(splsda.liver,name.var=1:3116, contrib = "max")
    
    #test legend
    plotLoadings(splsda.liver,legend=TRUE, contrib = "max")
    plotLoadings(splsda.liver,legend=FALSE, contrib = "max")
    
    #test legend.color
    plotLoadings(splsda.liver,legend.color=c("red","green","turquoise","blue"), contrib = "max")
    plotLoadings(splsda.liver,legend.color=1:4, contrib = "max")
    
    #test legend.title
    plotLoadings(splsda.liver,legend.title="legend", contrib = "max")
    plotLoadings(splsda.liver,legend.title="legend",legend=FALSE, contrib = "max")
    
    #test title
    plotLoadings(splsda.liver,title="Plot", contrib = "max")
    
    #test plot
    plotLoadings(splsda.liver,plot=FALSE,title="plot", contrib = "max")
    
    
    
    ##### multilevel
    data(vac18)
    X <- vac18$genes
    Y <- vac18$stimulation
    
    design <- data.frame(sample = vac18$sample,
    stimul = vac18$stimulation)
    
    res.1level <- splsda(X, Y = design[,2], ncomp = 3, multilevel = design[,1],
        keepX = c(30, 137, 123))
    
    
    plotLoadings(res.1level,contrib='max')
    plotLoadings(res.1level,contrib='min')
    
    #test method
    plotLoadings(res.1level,method='mean', contrib = "max")
    plotLoadings(res.1level,method='median', contrib = "max")
    
    #test comp
    plotLoadings(res.1level,method='mean',comp=1, contrib = "max")
    plotLoadings(res.1level,method='mean',comp=2, contrib = "max")
    plotLoadings(res.1level,method='mean',comp=3, contrib = "max")
    plotLoadings(res.1level,method='median',comp=1, contrib = "max")
    plotLoadings(res.1level,method='median',comp=2, contrib = "max")
    plotLoadings(res.1level,method='median',comp=3, contrib = "max")
    
    #test ties
    plotLoadings(res.1level,show.ties=FALSE, contrib = "max")
    plotLoadings(res.1level,show.ties=TRUE, contrib = "max")
    
    #test cex.name
    plotLoadings(res.1level,cex.name=0.5, contrib = "max")
    plotLoadings(res.1level,cex.name=1.5, contrib = "max")
    
    #test cex.legend
    plotLoadings(res.1level,cex.legend=0.5, contrib = "max")
    plotLoadings(res.1level,cex.legend=1, contrib = "max")
    
    #test name.var
    plotLoadings(res.1level,name.var=substr(liver.toxicity$gene.ID[1:ncol(X), 'gene.title'],1,20))
    plotLoadings(res.1level,name.var=substr(liver.toxicity$gene.ID[1:ncol(X), 'gene.title'],1,20), contrib = "max")
    plotLoadings(res.1level,name.var=1:ncol(X), contrib = "max")
    
    #test legend
    plotLoadings(res.1level,legend=TRUE, contrib = "max")
    plotLoadings(res.1level,legend=FALSE, contrib = "max")
    
    #test legend.color
    plotLoadings(res.1level,legend.color=c("red","green","turquoise","blue"), contrib = "max")
    plotLoadings(res.1level,legend.color=1:4, contrib = "max")
    
    #test legend.title
    plotLoadings(res.1level,legend.title="legend", contrib = "max")
    plotLoadings(res.1level,legend.title="legend",legend=FALSE, contrib = "max")
    
    #test title
    plotLoadings(res.1level,title="Plot", contrib = "max")
    
    #test plot
    a=plotLoadings(res.1level,plot=FALSE,title="plot", contrib = "max")
    head(a)
    
    
    ##### sgccda
    data(nutrimouse)
    Y = nutrimouse$diet
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    nutrimouse.sgccda <- wrapper.sgccda(X = data,
    Y = Y,
    design = design,
    keepX = list(gene = c(10,10), lipid = c(15,15)),
    ncomp = 2,
    scheme = "centroid",
    verbose = FALSE,
    bias = FALSE)

    plotLoadings(nutrimouse.sgccda,block=1, legend = FALSE)
    plotLoadings(nutrimouse.sgccda,block=2, legend = FALSE)
    plotLoadings(nutrimouse.sgccda,block=1:2, legend = FALSE)
    plotLoadings(nutrimouse.sgccda,block=1:2, legend = FALSE, title = "YAY")
    plotLoadings(nutrimouse.sgccda,block=1:2, legend = FALSE, title = "YAY", subtitle = c("a","Flo"))
    plotLoadings(nutrimouse.sgccda,block=1:2, contrib = "max", legend = FALSE)
    plotLoadings(nutrimouse.sgccda,block=1:2, contrib = "max")
    plotLoadings(nutrimouse.sgccda,block=2, contrib = "max")
    
    
    plotLoadings(nutrimouse.sgccda,contrib='max')
    plotLoadings(nutrimouse.sgccda,contrib='min')
    
    #test method
    plotLoadings(nutrimouse.sgccda,method='mean', contrib = "max")
    plotLoadings(nutrimouse.sgccda,method='median', contrib = "max")
    
    #test comp
    plotLoadings(nutrimouse.sgccda,method='mean',comp=1, contrib = "max")
    plotLoadings(nutrimouse.sgccda,method='mean',comp=2, contrib = "max")
    #try(plotLoadings(nutrimouse.sgccda,method='mean',comp=3)) # error comp>ncomp
    plotLoadings(nutrimouse.sgccda,method='median',comp=1, contrib = "max")
    plotLoadings(nutrimouse.sgccda,method='median',comp=2, contrib = "max")
    #try(plotLoadings(nutrimouse.sgccda,method='median',comp=3)) # error comp>ncomp
    
    #test ties
    plotLoadings(nutrimouse.sgccda,show.ties=FALSE, contrib = "max")
    plotLoadings(nutrimouse.sgccda,show.ties=TRUE, contrib = "max")
    
    #test cex.name
    plotLoadings(nutrimouse.sgccda,cex.name=0.5, contrib = "max")
    plotLoadings(nutrimouse.sgccda,cex.name=1.5, contrib = "max")
    
    #test cex.legend
    plotLoadings(nutrimouse.sgccda,cex.legend=0.5, contrib = "max")
    plotLoadings(nutrimouse.sgccda,cex.legend=1, contrib = "max")
    
    #test name.var
    plotLoadings(nutrimouse.sgccda,name.var=liver.toxicity$gene.ID[1:ncol(nutrimouse$gene), 'gene.title'], contrib = "max",block=1)
    plotLoadings(nutrimouse.sgccda,name.var=liver.toxicity$gene.ID[1:ncol(nutrimouse$gene), 'gene.title'], contrib = "max",block=1)

    plotLoadings(nutrimouse.sgccda,name.var=1:ncol(nutrimouse$gene), contrib = "max",block=1)
    
    #test legend
    plotLoadings(nutrimouse.sgccda,legend=TRUE, contrib = "max")
    plotLoadings(nutrimouse.sgccda,legend=FALSE, contrib = "max")
    
    #test legend.color
    plotLoadings(nutrimouse.sgccda,legend.color=c("red","green","turquoise","blue","black"), contrib = "max")
    plotLoadings(nutrimouse.sgccda,legend.color=1:5, contrib = "max")
    
    #test legend.title
    plotLoadings(nutrimouse.sgccda,legend.title="legend", contrib = "max")
    plotLoadings(nutrimouse.sgccda,legend.title="legend",legend=FALSE, contrib = "max")
    
    #test title
    plotLoadings(nutrimouse.sgccda,title="Plot", contrib = "max")
    
    #test plot
    plotLoadings(nutrimouse.sgccda,plot=FALSE,title="plot", contrib = "max")
    
}
par(opar)

