### splsda, plsda

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
plotContrib(splsda.liver,ties=FALSE)
plotContrib(splsda.liver,ties=TRUE)

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

res.1level <- multilevel(X, ncomp = 3, design = design,
                         method = "splsda", keepX = c(30, 137, 123))


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
plotContrib(res.1level,ties=FALSE)
plotContrib(res.1level,ties=TRUE)

#test cex.name
plotContrib(res.1level,cex.name=0.5)
plotContrib(res.1level,cex.name=1.5)

#test cex.legend
plotContrib(res.1level,cex.legend=0.5)
plotContrib(res.1level,cex.legend=1)

#test name.var
plotContrib(res.1level,name.var=liver.toxicity$gene.ID[, 'gene.title'])
plotContrib(res.1level,name.var=1:3116)

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
