#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))

# contribution on comp 1, based on the median.
# Colors indicate the group in which the median expression is maximal
plotContrib(splsda.liver, comp = 1, method = 'median')

# contribution on comp 2, based on median.
#Colors indicate the group in which the median expression is maximal
plotContrib(splsda.liver, comp = 2, method = 'median')

# contribution on comp 2, based on median.
# Colors indicate the group in which the median expression is minimal
plotContrib(splsda.liver, comp = 2, method = 'median', contrib = 'min')

# changing the name to gene names
# if the user input a name.var but names(name.var) is NULL,
# then a warning will be output and assign names of name.var to colnames(X)
# this is to make sure we can match the name of the selected variables to the contribution plot.
name.var = liver.toxicity$gene.ID[, 'geneBank']
length(name.var)
plotContrib(splsda.liver, comp = 2, method = 'median', name.var = name.var,main="Liver data")

# if names are provided: ok, even when NAs
name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)
plotContrib(splsda.liver, comp = 2, method = 'median',
name.var = name.var, cex.name = 0.5)

#missing names of some genes? complete with the original names
plotContrib(splsda.liver, comp = 2, method = 'median',
name.var = name.var, cex.name = 0.5,complete.name.var=TRUE)

# look at the contribution (median) for each variable
plot.contrib = plotContrib(splsda.liver, comp = 2, method = 'median', plot = FALSE)
head(plot.contrib$contrib)
# change the title of the legend and title name
plotContrib(splsda.liver, comp = 2, method = 'median', legend.title = 'Time',
main = 'Contribution plot')

# no legend
plotContrib(splsda.liver, comp = 2, method = 'median', legend = FALSE)

# change the color of the legend
plotContrib(splsda.liver, comp = 2, method = 'median', legend.color = c(1:4))


data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample,
stimul = vac18$stimulation)

# multilevel sPLS-DA model
res.1level = splsda(X, ncomp = 3, multilevel = design,
 keepX = c(30, 137, 123))


name.var = vac18$tab.prob.gene[, 'Gene']
names(name.var) = colnames(X)

plotContrib(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
name.var = name.var, cex.name = 0.2)

# too many transcripts? only output the top ones
plotContrib(res.1level, comp = 2, method = 'median', legend.title = 'Stimu',
name.var = name.var, cex.name = 0.5, ndisplay = 60)

# breast tumors
# ---
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2,near.zero.var=TRUE)

name.var = as.character(breast.tumors$genes$name)
names(name.var) = colnames(X)

# with gene IDs, showing the top 60
plotContrib(plsda.breast, contrib = 'max', comp = 1, method = 'median',
ndisplay = 60,
name.var = name.var,
cex.name = 0.6,
legend.color = color.mixo(1:2))

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$treatment[, 4]

plsda.liver <- plsda(X, Y, ncomp = 2)
plotIndiv(plsda.liver, ind.names = Y, plot.ellipse = TRUE)


name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)

plotContrib(plsda.liver, contrib = 'max', comp = 1, method = 'median', ndisplay = 100,
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
    plotContrib(splsda.liver,contrib='max')
    plotContrib(splsda.liver,contrib='min')
    
    #test method
    plotContrib(splsda.liver,method='mean')
    plotContrib(splsda.liver,method='median')
    
    #test comp
    plotContrib(splsda.liver,method='mean',comp=1)
    plotContrib(splsda.liver,method='mean',comp=2)
    plotContrib(splsda.liver,method='mean',comp=3)
    plotContrib(splsda.liver,method='median',comp=1)
    plotContrib(splsda.liver,method='median',comp=2)
    plotContrib(splsda.liver,method='median',comp=3)
    
    #test ties
    plotContrib(splsda.liver,show.ties=FALSE)
    plotContrib(splsda.liver,show.ties=TRUE)
    
    #test cex.name
    plotContrib(splsda.liver,cex.name=0.5)
    plotContrib(splsda.liver,cex.name=1.5)
    
    #test cex.legend
    plotContrib(splsda.liver,cex.legend=0.5)
    plotContrib(splsda.liver,cex.legend=1)
    
    #test name.var
    plotContrib(splsda.liver,name.var=liver.toxicity$gene.ID[, 'gene.title'])
    plotContrib(splsda.liver,name.var=1:3116)
    
    #test legend
    plotContrib(splsda.liver,legend=TRUE)
    plotContrib(splsda.liver,legend=FALSE)
    
    #test legend.color
    plotContrib(splsda.liver,legend.color=c("red","green","turquoise","blue"))
    plotContrib(splsda.liver,legend.color=1:4)
    
    #test legend.title
    plotContrib(splsda.liver,legend.title="legend")
    plotContrib(splsda.liver,legend.title="legend",legend=FALSE)
    
    #test main
    plotContrib(splsda.liver,main="Plot")
    
    #test plot
    plotContrib(splsda.liver,plot=FALSE,main="plot")
    
    
    
    ##### multilevel
    data(vac18)
    X <- vac18$genes
    Y <- vac18$stimulation
    
    design <- data.frame(sample = vac18$sample,
    stimul = vac18$stimulation)
    
    res.1level <- splsda(X, ncomp = 3, multilevel = design,
        keepX = c(30, 137, 123))
    
    
    plotContrib(res.1level,contrib='max')
    plotContrib(res.1level,contrib='min')
    
    #test method
    plotContrib(res.1level,method='mean')
    plotContrib(res.1level,method='median')
    
    #test comp
    plotContrib(res.1level,method='mean',comp=1)
    plotContrib(res.1level,method='mean',comp=2)
    plotContrib(res.1level,method='mean',comp=3)
    plotContrib(res.1level,method='median',comp=1)
    plotContrib(res.1level,method='median',comp=2)
    plotContrib(res.1level,method='median',comp=3)
    
    #test ties
    plotContrib(res.1level,show.ties=FALSE)
    plotContrib(res.1level,show.ties=TRUE)
    
    #test cex.name
    plotContrib(res.1level,cex.name=0.5)
    plotContrib(res.1level,cex.name=1.5)
    
    #test cex.legend
    plotContrib(res.1level,cex.legend=0.5)
    plotContrib(res.1level,cex.legend=1)
    
    #test name.var
    plotContrib(res.1level,name.var=liver.toxicity$gene.ID[1:ncol(X), 'gene.title'])
    plotContrib(res.1level,name.var=1:ncol(X))
    
    #test legend
    plotContrib(res.1level,legend=TRUE)
    plotContrib(res.1level,legend=FALSE)
    
    #test legend.color
    plotContrib(res.1level,legend.color=c("red","green","turquoise","blue"))
    plotContrib(res.1level,legend.color=1:4)
    
    #test legend.title
    plotContrib(res.1level,legend.title="legend")
    plotContrib(res.1level,legend.title="legend",legend=FALSE)
    
    #test main
    plotContrib(res.1level,main="Plot")
    
    #test plot
    plotContrib(res.1level,plot=FALSE,main="plot")
    
    
    
    ##### sgccda
    data(nutrimouse)
    Y = nutrimouse$diet
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    nutrimouse.sgccda <- wrapper.sgccda(X = data,
    Y = Y,
    design = design,
    keepX = list(c(10,10), c(15,15)),
    ncomp = c(2, 2),
    scheme = "centroid",
    verbose = FALSE,
    bias = FALSE)
    
    plotContrib(nutrimouse.sgccda,block=2)
    
    
    plotContrib(nutrimouse.sgccda,contrib='max')
    plotContrib(nutrimouse.sgccda,contrib='min')
    
    #test method
    plotContrib(nutrimouse.sgccda,method='mean')
    plotContrib(nutrimouse.sgccda,method='median')
    
    #test comp
    plotContrib(nutrimouse.sgccda,method='mean',comp=1)
    plotContrib(nutrimouse.sgccda,method='mean',comp=2)
    #try(plotContrib(nutrimouse.sgccda,method='mean',comp=3)) # error comp>ncomp
    plotContrib(nutrimouse.sgccda,method='median',comp=1)
    plotContrib(nutrimouse.sgccda,method='median',comp=2)
    #try(plotContrib(nutrimouse.sgccda,method='median',comp=3)) # error comp>ncomp
    
    #test ties
    plotContrib(nutrimouse.sgccda,show.ties=FALSE)
    plotContrib(nutrimouse.sgccda,show.ties=TRUE)
    
    #test cex.name
    plotContrib(nutrimouse.sgccda,cex.name=0.5)
    plotContrib(nutrimouse.sgccda,cex.name=1.5)
    
    #test cex.legend
    plotContrib(nutrimouse.sgccda,cex.legend=0.5)
    plotContrib(nutrimouse.sgccda,cex.legend=1)
    
    #test name.var
    plotContrib(nutrimouse.sgccda,name.var=liver.toxicity$gene.ID[1:ncol(nutrimouse$gene), 'gene.title'])
    plotContrib(nutrimouse.sgccda,name.var=1:ncol(nutrimouse$gene))
    
    #test legend
    plotContrib(nutrimouse.sgccda,legend=TRUE)
    plotContrib(nutrimouse.sgccda,legend=FALSE)
    
    #test legend.color
    plotContrib(nutrimouse.sgccda,legend.color=c("red","green","turquoise","blue","black"))
    plotContrib(nutrimouse.sgccda,legend.color=1:5)
    
    #test legend.title
    plotContrib(nutrimouse.sgccda,legend.title="legend")
    plotContrib(nutrimouse.sgccda,legend.title="legend",legend=FALSE)
    
    #test main
    plotContrib(nutrimouse.sgccda,main="Plot")
    
    #test plot
    plotContrib(nutrimouse.sgccda,plot=FALSE,main="plot")
    
}

