
##### spls

# plotIndiv

X <- pentosys$gene
Y <- pentosys$metabolom
pento.spls <- spls(X, Y, ncomp = 3,
                      keepX = c(50, 50, 50), keepY = c(10, 10, 10))


col <- c("blue", "red", "darkgreen", "darkviolet", "turquoise")
pch <- c(15, 16, 17, 18, 19)

cex <- c(1, 1.2, 1.4, 1.6, 1.8)

groupe=pentosys$conditions$condition

plotIndiv(pento.spls, comp = c(1,2), ind.names = rownames(pentosys$conditions),
          group = groupe)

plotIndiv(pento.spls, comp = c(1,3), ind.names = FALSE,
          group = groupe)



plotIndiv(pento.spls,rep.space = "Y-variate",cex=cex,
          ind.names=FALSE,style="lattice", plot.ellipse=TRUE,
          group=groupe)


plotIndiv(pento.spls, ind.names= FALSE, add.legend=TRUE,abline.line = TRUE,col=col,group=groupe)
plotIndiv(pento.spls,add.legend=TRUE,ind.names=TRUE,style="lattice",group=groupe)
plotIndiv(pento.spls,add.legend=TRUE,ind.names=TRUE,style="graphics",group=groupe)

plotIndiv(pento.spls,col=col,
          add.legend=list(space="right",title="Groupe",text=list(unique(as.character(groupe))),point=list(col=col)),
          ind.names=FALSE,style="lattice",group=groupe)


# plot.perf
pento.perf <- perf(pento.spls, validation = "Mfold")

plot(pento.perf,criterion = "R2")
plot(pento.perf,criterion = "Q2")
plot(pento.perf,criterion = "Q2.total")
plot(pento.perf,criterion = "MSEP") 

# cim
cim(pento.spls)

cim(pento.spls,mapping='Y',sample.sideColors=rep(col,each=3),sample.names = FALSE,
    legend=list( legend = unique(pentosys$conditions$condition), col = col,title = "Condition"))

# select_var
selectVar(pento.spls, comp = 2)


##### rcc

rcc(pentosys$gene,pentosys$metabolom) ## pbm





##### plsda

X <- pentosys$gene
Y <- pentosys$conditions$condition
pento.plsda=plsda(X, Y)

# plotIndiv

plotIndiv(pento.plsda, ind.names = TRUE,style="lattice")
plotIndiv(pento.plsda, ind.names = FALSE,ellipse.level=0.5,plot.ellipse=TRUE)
plotIndiv(pento.plsda,rep.space="Y-variate", ind.names = F,style="graphics",add.legend=TRUE)

# plotperf
pento.perf = perf(pento.plsda,validation="loo",method.predict="all")

plot(pento.perf,criterion = "R2")
plot(pento.perf,criterion = "Q2")
plot(pento.perf,criterion = "Q2.total")
plot(pento.perf,criterion = "MSEP") 


# plotContrib (#mauvaise version de plotContrib)
plotContrib(pento.plsda,contrib='max',method='mean',comp=2)
plotContrib(pento.plsda,contrib='min'method='median',comp=2)

# cim
cim(pento.plsda,mapping='Y')
