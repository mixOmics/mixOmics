data(multidrug)

#test ncomp
multi.pca=pca(multidrug$ABC.trans,ncomp=1)
plotIndiv(multi.pca)####pbm avec plotIndiv
plot.perf(multi.pca)###pbm avec plotperf
plot.perf(multi.pca,criterion="R2")### plotperf
plot.pca(multi.pca)

multi.pca=pca(multidrug$ABC.trans,ncomp=2)
plotIndiv(multi.pca)
plot.pca(multi.pca)
plotContrib(multi.pca)
multi.pca=pca(multidrug$ABC.trans,ncomp=48)
plotIndiv(multi.pca,comp=c(1,48))
plot.pca(multi.pca)

#test center ??? aucune difference
multi.pca=pca(multidrug$ABC.trans,center=TRUE)
plotIndiv(multi.pca)
biplot(multi.pca)

multi.pca=pca(multidrug$ABC.trans,center=FALSE)
plotIndiv(multi.pca)
biplot(multi.pca)

#test scale
multi.pca=pca(multidrug$ABC.trans,scale=TRUE)
plotIndiv(multi.pca)
biplot(multi.pca)

multi.pca=pca(multidrug$ABC.trans,scale=FALSE)
plotIndiv(multi.pca)
biplot(multi.pca)

#test comp.tol
multi.pca=pca(multidrug$ABC.trans,ncomp=10,tol=1e-1)
plotIndiv(multi.pca)
plot.pca(multi.pca)
biplot(multi.pca)

multi.pca=pca(multidrug$ABC.trans,ncomp=10,tol=1e-5)
plotIndiv(multi.pca)
plot.pca(multi.pca)
biplot(multi.pca)

multi.pca=pca(multidrug$ABC.trans,ncomp=10,tol=10)
plotIndiv(multi.pca)
plot.pca(multi.pca)
biplot(multi.pca)

#test max.iter
multi.pca=pca(multidrug$ABC.trans,ncomp=10,max.iter=10) 
plotIndiv(multi.pca)
plot.pca(multi.pca)
biplot(multi.pca)

multi.pca=pca(multidrug$ABC.trans,ncomp=10,max.iter=5000) 
plotIndiv(multi.pca)
plot.pca(multi.pca)
biplot(multi.pca)

multi.pca=pca(multidrug$ABC.trans,ncomp=10,max.iter=500) 
plotIndiv(multi.pca)
plot.pca(multi.pca)
biplot(multi.pca)

