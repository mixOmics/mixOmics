###############################
############# PLS #############
###############################
#### Difference because of not the same starting point between the 2 algo "pls" and "spls" in mixOmics package when tol = 1e-06
rm(list=ls())
setwd("/Users/florian/Work/git/package-mixOmics/")
library(mixOmics)


load("test_scripts/meta.block_test/Fibro-ESC-iPSC.6exp.167samples.light.Rdata") #load data, type.id, exp among others



source("mixOmics/R/check_entry.R")
source("mixOmics/R/helpers.R")
source("mixOmics/R/meta.block.spls.R")
source("mixOmics/R/meta.spls.hybrid.R")
source("mixOmics/R/mixOmics.R")
source("mixOmics/R/pls.R")
source("mixOmics/R/plsda.R")
source("mixOmics/R/spls.R")
source("mixOmics/R/splsda.R")
source("mixOmics/R/meta.pls.R")
source("mixOmics/R/meta.plsda.R")
source("mixOmics/R/meta.spls.R")
source("mixOmics/R/meta.splsda.R")



#res.sgcca=meta.spls.hybrid.sgcca(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,5),ncomp=3,near.zero.var=FALSE,tol=1e-25)
Y.mat=unmap(type.id)
rownames(Y.mat)=rownames(data)

#wraper.meta.spls.hybrid
res=wrapper.meta.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,tol=1e-25)

## ======================================================================
## =========================        data for one study
## ======================================================================
ind=which(exp=="Briggs")
study.light=exp[ind]
type.id.light=type.id[ind]
data.light=data[ind,]
Y.mat.light=unmap(type.id.light)
rownames(Y.mat.light)=rownames(data.light)

#wraper.pls
res=wrapper.pls(X=data.light,Y=Y.mat.light,ncomp=3,near.zero.var=FALSE)


#wraper.plsda
res=wrapper.plsda(X=data.light,Y=type.id.light,ncomp=3,near.zero.var=FALSE)


#wraper.spls
res=wrapper.spls(X=data.light,Y=Y.mat.light,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15))


#wraper.splsda
res=wrapper.splsda(X=data.light,Y=type.id.light,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15))





res=mixOmics(data.light,type.id.light,ncomp=3) #plsda
res=mixOmics(data.light,type.id.light,ncomp=3,keepX=c(10,5,15))#splsda
res=mixOmics(data.light,Y.mat.light,ncomp=3) #pls
res=mixOmics(data.light,Y.mat.light,ncomp=3,keepX=c(10,5,15)) #spls


res=mixOmics(data,type.id,ncomp=3,study=exp) #pls
res=mixOmics(X=data,Y=type.id,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,15),ncomp=3)

## ======================================================================
## =========================        checks for wrong inputs
## ======================================================================

#Y too long
res=mixOmics(data,type.id.light,ncomp=3,study=exp) #pls
# bad keepX.constraint
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576",100)),keepX=c(10,15),ncomp=3)
# bad ncomp
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000006576")),keepX=c(10,15),ncomp=2)


#with gene names in keepX.constraint
#res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)
which(res.spls.hybrid$loadings$X[,1]!=0)
which(res.spls.hybrid$loadings$X[,2]!=0)

#with numbers in keepX.constraint
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),keepX=c(100,50),ncomp=3)
which(res.spls.hybrid$loadings$X[,1]!=0)
which(res.spls.hybrid$loadings$X[,2]!=0)




system.time(wrapper.meta.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE,tol=1e-25))
system.time(mixOmics(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))


system.time(splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))
system.time(mixOmics(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))

system.time(splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE))
system.time(mixOmics(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE))
system.time(wrapper.splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE))

