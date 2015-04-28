###############################
############# PLS #############
###############################
#### Difference because of not the same starting point between the 2 algo "pls" and "spls" in mixOmics package when tol = 1e-06
#source("/Users/florian/Work/git/package-mixomics/test_scripts/meta.block_test/chack_mixOmics.R")
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

Y.mat=unmap(type.id)
rownames(Y.mat)=rownames(data)


#res.sgcca=meta.spls.hybrid.sgcca(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,5),ncomp=3,near.zero.var=FALSE,tol=1e-25)

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



res=mixOmics(data,type.id,ncomp=3,study=exp) #meta.plsda
res=mixOmics(X=data,Y=type.id,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,15),ncomp=3)#meta.splsda

res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,ncomp=3)#meta.pls

#with gene names in keepX.constraint
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)#meta.spls
which(res.spls.hybrid$loadings$X[,1]!=0)
which(res.spls.hybrid$loadings$X[,2]!=0)

#with numbers in keepX.constraint
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),keepX=c(100,50),ncomp=3)#meta.spls
which(res.spls.hybrid$loadings$X[,1]!=0)
which(res.spls.hybrid$loadings$X[,2]!=0)




system.time(wrapper.meta.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE,tol=1e-25))
system.time(mixOmics(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))


system.time(splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))
system.time(mixOmics(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))

system.time(splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE))
system.time(wrapper.splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE))
system.time(mixOmics(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,scale=TRUE))
#note: difference of time due to different number of iterations (due to the convergence that we have now):
sp1=splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE)
sp2=mixOmics(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,scale=TRUE)
sp1$iter;sp2$iter

## ======================================================================
## =========================        checks for wrong inputs
## ======================================================================

#Y too long
res=mixOmics(data,type.id.light,ncomp=3,study=exp) #pls
# bad keepX.constraint
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576",100)),keepX=c(10,15),ncomp=3)
# bad ncomp
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000006576")),keepX=c(10,15),ncomp=2)
#extra parameters
res=mixOmics(data.light,Y.mat.light,ncomp=3,keepX=c(10,5,15),scheme="centroid") #spls





## ======================================================================
## =========================        data is a list
## ======================================================================


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

A=list(X=data,Y=Y.mat)

res=mixOmics(X=A,indY=2,ncomp=c(2,2))

res=mixOmics(X=A,indY=2,keepX.constraint=list(X=list(1:10)),ncomp=c(2,2))
res=mixOmics(X=A,indY=2)

#OK keep variable 5 and 10 on the 2first comp of block1, keep variable 2 on comp1 for block 2. completed by keeping all variable on comp2 for block2
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX.constraint=list(c(5,10),2),ncomp=c(2,2))

#OK keep variable 5 on comp1 and 10 on comp2 of block1; keep variable 2 on comp1 for block 2. completed by keeping all variable on comp3 for block1 and comp2/3 for block2
res1=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX.constraint=list(c(5,10),2),ncomp=c(3,3))

#OK keep variable 5 and 10 on the first comp of block1, variable 3 on the second comp of block1; keep variable 2 on comp1 for block 2. completed by keeping all variable on comp3 for block1 and comp2/3 for block2
res2=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX.constraint=list(list(c(5,10),3),list(2)),ncomp=c(3,3))

which(res1$loadings$X[,1]!=0)
which(res1$loadings$X[,2]!=0)
sum(res1$loadings$X[,3]!=0)

which(res2$loadings$X[,1]!=0)
which(res2$loadings$X[,2]!=0)
sum(res2$loadings$X[,3]!=0)


## ======================================================================
## =========================        checks for wrong inputs
## ======================================================================
res=mixOmics(X=A) #missing stuff

res=mixOmics(X=Y.mat,Y=A)# bad Y

res=mixOmics(X=A,indY=2,keepX.constraint=list(X=list(1:10,1:5),list(1:4)))#bad keepX.constraint

#check RGCCA
res=mixOmics(X=A,tau=1,mode="canonical") #error length tau
res=mixOmics(X=A,tau=1,mode="regression") # error length tau (and mode)
res=mixOmics(X=A,tau=c(1,1),mode="regression") # error mode
res=mixOmics(X=A,tau=c(1,1),mode="canonical",init="svd") #error svd
res=mixOmics(X=A,tau=c(1,1),mode="canonical",scheme="bla") #error scheme
res=mixOmics(X=A,tau=c(1,1),mode="canonical",study=rep(1,nrow(data))) #error study
res=mixOmics(X=A,tau=c(1,1),mode="canonical",study=factor(rep(1,nrow(data)))) #OK because study is factor with only one level


res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX=list(list(1:5),list(1:10)),ncomp=c(1,1)) #error keepX
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX=list(c(5,10)),ncomp=c(1,1)) #error keepX (keepX=c(5,10) for block1, but ncomp=1
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX=list(5,10),ncomp=c(1,1)) #error keepX, keepY=10 for Y, should be less than 3
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX.constraint=list(list(c(5,10),3)),ncomp=c(1,1)) #error keepX.constraint. Keep variables 5,10 on comp1 and 3 on comp2, but ncomp=1

res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX=list(5,2),keepX.constraint=list(list(c(5,10))),ncomp=c(2,2)) #error keepX+keepX.constraint !=ncomp


res=mixOmics(X=A,indY=2,keepX.constraint=list(X=list(1:10)),ncomp=c(3,1)) # OK
res=mixOmics(X=A,indY=2,tau=c(1,1),ncomp=c(3,1)) #OK





## ======================================================================
## =========================        to do
## ======================================================================

# what if rgcca with one of the blocks being a factor?
# what if X is a list and Y a factor? it should be a meta.block.splsda, not a the moment: error Y must be a numeric matrix

res=mixOmics(X=list(X=data,Y=type.id),tau=1) #A[[2]] must be numeric matrix
res=mixOmics(X=A,Y=type.id) #Y must be numeric matrix
res=mixOmics(X=list(data),Y=type.id) #Y must be numeric matrix
res=mixOmics(X=list(data),Y=unmap(type.id)) #OK


