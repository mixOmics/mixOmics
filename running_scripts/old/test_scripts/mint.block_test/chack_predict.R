###############################
############# PLS #############
###############################
#### Difference because of not the same starting point between the 2 algo "pls" and "spls" in mixOmics package when tol = 1e-06
#source("/Users/florian/Work/git/package-mixomics/test_scripts/mint.block_test/chack_mixOmics.R")
rm(list=ls())
setwd("/Users/florian/Work/git/package-mixOmics/")
library(mixOmics)


load("test_scripts/mint.block_test/Fibro-ESC-iPSC.6exp.167samples.light.Rdata") #load data, type.id, exp among others



source("mixOmics/R/check_entry.R")
source("mixOmics/R/helpers.R")
source("mixOmics/R/sparse.mint.block.R")

source("mixOmics/R/mint.spls.hybrid.R")
source("mixOmics/R/pls.R")
source("mixOmics/R/plsda.R")
source("mixOmics/R/spls.R")
source("mixOmics/R/splsda.R")
source("mixOmics/R/mint.pls.R")
source("mixOmics/R/mint.plsda.R")
source("mixOmics/R/mint.spls.R")
source("mixOmics/R/mint.splsda.R")

source("mixOmics/R/wrapper.sparse.mint.block.R")
source("mixOmics/R/block.pls.R")
source("mixOmics/R/block.spls.R")
source("mixOmics/R/block.plsda.R")
source("mixOmics/R/block.splsda.R")
source("mixOmics/R/mint.block.pls.R")
source("mixOmics/R/mint.block.spls.R")
source("mixOmics/R/mint.block.plsda.R")
source("mixOmics/R/mint.block.splsda.R")


source("mixOmics/R/wrapper.rgcca.R")
source("mixOmics/R/wrapper.sparse.rgcca.R")
source("mixOmics/R/wrapper.sgcca.R") # need to be checked or deleted, should be same as block.pls
source("mixOmics/R/wrapper.sgccda.R") #same as block.splsda


source("mixOmics/R/mixOmics.R")




source("mixOmics/R/predict.mint.block.pls.R")






Y.mat=unmap(type.id)
rownames(Y.mat)=rownames(data)
ind=which(exp=="Briggs")
study.light=exp[ind]
type.id.light=type.id[ind]
data.light=data[ind,]
Y.mat.light=unmap(type.id.light)
rownames(Y.mat.light)=rownames(data.light)

A=list(X=data,Y=Y.mat)
A.light=list(X=data.light,Y=Y.mat.light)



## =========================================================================================================
## =========================================================================================================
## ===========================                  test_wrappers               ==============================##
## =========================================================================================================
## =========================================================================================================


#wraper.mint.spls.hybrid, function not available to user, so not a proper check to conduct, no need to do a predict function either
res=wrapper.mint.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,tol=1e-25)
res=wrapper.mint.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,5),ncomp=3,near.zero.var=FALSE,tol=1e-25)


## ======================================================================
## ==      data for one study
## ======================================================================



#wraper.pls
res=wrapper.pls(X=data.light,Y=Y.mat.light,ncomp=5)
pred=predict(res,newdata=data.light)

#wraper.spls
res=wrapper.spls(X=data.light,Y=Y.mat.light,ncomp=3,keepX=c(10,5,15))
pred=predict(res,newdata=data.light)

#wraper.plsda
res=wrapper.plsda(X=data.light,Y=type.id.light,ncomp=3)
pred=predict(res,newdata=data.light)

#wraper.splsda
res=wrapper.splsda(X=data.light,Y=type.id.light,ncomp=3,keepX=c(10,5,15))
pred=predict(res,newdata=data.light)



#wraper.block.pls
res=wrapper.block.pls(X=A.light,indY=2)
pred=predict(res,newdata=list(A.light[[1]]))
res=wrapper.block.pls(list(data),Y=Y.mat,ncomp=2)
pred=predict(res,newdata=list(data))

#wraper.block.spls
res=wrapper.block.spls(X=A.light,indY=2,keepX=list(block1=c(10,5,15),block2=c(3,2)),ncomp=c(3,3))
pred=predict(res,newdata=list(A.light[[1]]))
res=wrapper.block.spls(list(data),Y=Y.mat,keepX=list(block1=c(10,5,15),block2=c(3,2)),ncomp=3)
pred=predict(res,newdata=list(data))

#wraper.block.plsda
res=wrapper.block.plsda(X=list(X=data,Y=type.id),indY=2,ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.block.plsda(list(data),Y=type.id,ncomp=2)
pred=predict(res,newdata=list(data))

#wraper.block.splsda
res=wrapper.block.splsda(X=list(X=data,Y=type.id),keepX=list(block1=c(10,5)),indY=2,ncomp=c(3,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.block.splsda(X=list(X=data,Y=type.id),keepX=list(block1=c(10,5)),indY=2,ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))


# same results for sgccda and block.splsda. outputs are different though
res=wrapper.sgccda(X=data,Y=type.id,keepA=c(10,5,10),ncomp=3)
res2=wrapper.block.splsda(X=list(X=data),Y=type.id,keepX=c(10,5,10),ncomp=3,mode="canonical")
#no predict function?


## ======================================================================
## ==      data for several study
## ======================================================================

#wraper.mint.pls
res=wrapper.mint.pls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE,study=exp)
sample=sample(1:nrow(data))
pred=predict(res,newdata=data,study.test=exp)
pred2=predict(res,newdata=data[sample,],study.test=exp[sample])

all.equal(pred$predict[sample,,],pred2$predict)
all.equal(pred$variates[sample,],pred2$variates)
all.equal(pred$newdata[sample,],pred2$newdata)
## all GOOD

#wraper.mint.plsda
res=wrapper.mint.plsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,study=exp)
pred=predict(res,newdata=data,study.test=exp)

#wraper.mint.spls
res=wrapper.mint.spls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)
pred=predict(res,newdata=data,study.test=exp)

#wraper.mint.splsda
res=wrapper.mint.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)
pred=predict(res,newdata=data,study.test=exp)



#wraper.mint.block.pls
res=wrapper.mint.block.pls(X=list(X=data,Y=Y.mat),indY=2,ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.mint.block.pls(list(data),Y=Y.mat,ncomp=2)
pred=predict(res,newdata=list(X=data))

#wraper.mint.block.spls
res=wrapper.mint.block.spls(X=list(X=data,Y=Y.mat),indY=2,keepX=list(block1=c(10,5)),ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.mint.block.spls(list(data),Y=Y.mat,ncomp=2)
pred=predict(res,newdata=list(X=data))

#wraper.mint.block.plsda
res=wrapper.mint.block.plsda(X=list(X=data,Y=type.id),indY=2,ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.mint.block.plsda(list(data),Y=type.id,ncomp=2)
pred=predict(res,newdata=list(X=data))

#wraper.mint.block.splsda
res=wrapper.mint.block.splsda(X=list(X=data,Y=type.id),indY=2,keepX=list(block1=c(10,5)),ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.mint.block.splsda(list(data),Y=type.id,ncomp=2)
pred=predict(res,newdata=list(X=data))





## =========================================================================================================
## =========================================================================================================
## ===========================                  test_mixOmics               ==============================##
## =========================================================================================================
## =========================================================================================================



## ======================================================================
## ==      data for one study
## ======================================================================

res=wrapper.pls(data.light,Y.mat.light,ncomp=3) #pls
pred=predict(res,newdata=data.light)
res2=mixOmics(data.light,Y.mat.light,ncomp=3) #pls
pred2=predict(res2,newdata=data.light)

all.equal(res,res2)

res=mixOmics(data.light,Y.mat.light,ncomp=3,keepX=c(10,5,15)) #spls
pred=predict(res,newdata=data.light)
res=mixOmics(data.light,type.id.light,ncomp=3) #plsda
pred=predict(res,newdata=data.light)
res=mixOmics(data.light,type.id.light,ncomp=3,keepX=c(10,5,15))#splsda
pred=predict(res,newdata=data.light)


#block.pls
res=mixOmics(X=list(data),Y=unmap(type.id))
pred=predict(res,newdata=list(data))
res=mixOmics(X=A,indY=2,ncomp=c(2,2))
pred=predict(res,newdata=A) #error too many blocks
pred=predict(res,newdata=list(A[[1]]))
res=mixOmics(X=A,indY=2)
pred=predict(res,newdata=list(A[[1]]))


#block.spls
res=mixOmics(X=list(data),Y=unmap(type.id),keepX=c(10,5,15))
pred=predict(res,newdata=list(data))
res=mixOmics(X=list(data=data,Y=Y.mat),indY=2,keepX=list(c(10,5,15)),ncomp=c(3,3))
pred=predict(res,newdata=list(data))
res=mixOmics(X=A,indY=2,keepX.constraint=list(X=list(1:10)),ncomp=c(2,2))
pred=predict(res,newdata=list(A[[1]]))

# block.plsda
res=mixOmics(X=A,Y=type.id)
pred=predict(res,newdata=A,method=c("max.dist", "centroids.dist"))
res=mixOmics(X=list(data=data,Y=type.id),indY=2,ncomp=c(3,3))
pred=predict(res,newdata=list(data))
res=mixOmics(X=A,Y=type.id)
pred=predict(res,newdata=A,method=c("max.dist", "centroids.dist"))
res=mixOmics(X=list(data),Y=type.id)
pred=predict(res,newdata=list(data))

# block.splsda
res=mixOmics(X=list(data=data,Y=type.id),indY=2,keepX=list(c(10,5,15)),ncomp=c(3,3))
pred=predict(res,newdata=list(data))


system.time(splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))
system.time(mixOmics(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))

system.time(splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE))
system.time(wrapper.splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,scale=TRUE))
system.time(mixOmics(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,scale=TRUE))
#note: difference of time due to different number of iterations (due to the convergence that we have now) and maybe other things:
sp1=splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE)
sp2=mixOmics(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,scale=TRUE)
sp1$iter;sp2$iter



## ======================================================================
## ==      data for several study
## ======================================================================

#mint.pls
res=mixOmics(data,Y=Y.mat,ncomp=3,study=exp)
pred=predict(res,newdata=data,study=exp)

#mint.spls
res=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)
res1=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)#with gene names in keepX.constraint
res2=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),keepX=c(100,50),ncomp=3)#with numbers in keepX.constraint
all.equal(res1,res2)

pred=predict(res,newdata=data,study=exp)
pred=predict(res1,newdata=data,study=exp)
pred=predict(res2,newdata=data,study=exp)

res=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),ncomp=3)#mint.spls missing keepX is completed by pls-like
res$keepX
pred=predict(res,newdata=data,study=exp)

#mint.plsda
res=mixOmics(data,type.id,ncomp=3,study=exp)
pred=predict(res,newdata=data,study=exp)

#mint.splsda
res=mixOmics(X=data,Y=type.id,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,15),ncomp=3)
pred=predict(res,newdata=data,study=exp)



source("mixOmics/R/predict.mint.block.pls.R")

#block.pls
res=mixOmics(X=list(data),Y=unmap(type.id),study=exp)
pred=predict(res,newdata=list(data),study=exp)

#block.spls
res=mixOmics(X=list(data),Y=unmap(type.id),keepX=c(10,5,15),study=exp)
pred=predict(res,newdata=list(data),study=exp)
res=mixOmics(X=list(data=data,Y=Y.mat),indY=2,keepX=list(c(10,5,15)),ncomp=c(3,3),study=exp)
pred=predict(res,newdata=list(data),study=exp)
res=mixOmics(X=A,indY=2,keepX.constraint=list(X=list(1:10)),ncomp=c(3,1)) # OK
pred=predict(res,newdata=res$X)

# block.plsda
res=wrapper.mint.block.plsda(X=A,Y=type.id,study=exp,ncomp=c(2,2))
res=mixOmics(X=A,Y=type.id,study=exp,ncomp=c(2,2))
pred=predict(res,newdata=res$X,study.test=exp)

#block.plsda
res=mixOmics(X=list(data=data,Y=type.id),indY=2,ncomp=c(3,3),study=exp)
pred=predict(res,newdata=res$X,study=exp)

#block.splsda
res=mixOmics(X=list(data=data,Y=type.id),indY=2,keepX=list(c(10,5,15)),ncomp=c(3,3),study=exp)
pred=predict(res,newdata=res$X,study=exp)



## ======================================================================
## ==      RGCCA/sparse.RGCCA
## ======================================================================

#RGCCA
res=mixOmics(X=A,tau=c(1,1))
pred=predict(res,newdata=res$X,study=exp)




## =========================================================================================================
## =========================================================================================================
## ===========================                  wrong_inputs                ==============================##
## =========================================================================================================
## =========================================================================================================


## ======================================================================
## ==     X is a matrix
## ======================================================================


#Y too long
res=mixOmics(data,type.id.light,ncomp=3,study=exp) #plsda
# bad keepX.constraint
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576",100)),keepX=c(10,15),ncomp=3)
# bad ncomp
res.spls.hybrid=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000006576")),keepX=c(10,15),ncomp=2)
#extra parameters
res=mixOmics(data.light,Y.mat.light,ncomp=3,keepX=c(10,5,15),scheme="centroid") #spls


## ======================================================================
## ==     X is a list
## ======================================================================




#block.pls



#block.spls
res=mixOmics(X=A,indY=2,keepX.constraint=list(X=list(1:10,1:5),list(1:4)))#bad keepX.constraint




#block.plsda




#block.splsda





