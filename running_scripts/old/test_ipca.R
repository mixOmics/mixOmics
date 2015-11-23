data(multidrug)
data(liver.toxicity)

#test ncomp
liver.ipca=ipca(multidrug$ABC.trans)###
liver.ipca=ipca(liver.toxicity$gene)

liver.ipca=pca(liver.toxicity,ncomp=2)####
liver.ipca=pca(liver.toxicity$gene,ncomp=2)
plotIndiv(liver.ipca)
plot.pca(liver.ipca)
liver.ipca=pca(liver.toxicity$gene,ncomp=64)
plotIndiv(liver.ipca,comp=c(63,64))
plot.pca(liver.ipca)

#test mode
liver.ipca=ipca(liver.toxicity$gene,mode="deflation")
plotIndiv(liver.ipca)

liver.ipca=ipca(liver.toxicity$gene,mode="parallel") #??
plotIndiv(liver.ipca)

#test fun
liver.ipca=ipca(liver.toxicity$gene,fun="logcosh")
plotIndiv(liver.ipca)

liver.ipca=ipca(liver.toxicity$gene,fun="exp")
plotIndiv(liver.ipca)

#test w.init
mat=matrix(1:9,3,3)
liver.ipca=ipca(liver.toxicity$gene,w.init=mat)
plotIndiv(liver.ipca)

mat=matrix(9:17,3,3)
liver.ipca=ipca(liver.toxicity$gene,w.init=mat)
plotIndiv(liver.ipca)

mat=matrix(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),3,3)
liver.ipca=ipca(liver.toxicity$gene,w.init=mat)
plotIndiv(liver.ipca)

mat=matrix(rep(0,9),3,3)
liver.ipca=ipca(liver.toxicity$gene,w.init=mat) ###
plotIndiv(liver.ipca)

mat=matrix(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),3,3)
mat[2,1]=NA
liver.ipca=ipca(liver.toxicity$gene,w.init=mat)###
plotIndiv(liver.ipca)