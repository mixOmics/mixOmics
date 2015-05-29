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
source("mixOmics/R/sparse.meta.block.R")

source("mixOmics/R/meta.spls.hybrid.R")
source("mixOmics/R/pls.R")
source("mixOmics/R/plsda.R")
source("mixOmics/R/spls.R")
source("mixOmics/R/splsda.R")
source("mixOmics/R/meta.pls.R")
source("mixOmics/R/meta.plsda.R")
source("mixOmics/R/meta.spls.R")
source("mixOmics/R/meta.splsda.R")

source("mixOmics/R/wrapper.sparse.meta.block.R")
source("mixOmics/R/block.pls.R")
source("mixOmics/R/block.spls.R")
source("mixOmics/R/block.plsda.R")
source("mixOmics/R/block.splsda.R")
source("mixOmics/R/meta.block.pls.R")
source("mixOmics/R/meta.block.spls.R")
source("mixOmics/R/meta.block.plsda.R")
source("mixOmics/R/meta.block.splsda.R")


source("mixOmics/R/wrapper.rgcca.R")
source("mixOmics/R/wrapper.sparse.rgcca.R")
source("mixOmics/R/wrapper.sgcca.R") # need to be checked or deleted, should be same as block.pls
source("mixOmics/R/wrapper.sgccda.R") #same as block.splsda


source("mixOmics/R/mixOmics.R")




source("mixOmics/R/predict.meta.block.pls.R")






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




res=wrapper.meta.block.pls(X=list(X=data,Y=Y.mat),indY=2,ncomp=c(2,2))
res=wrapper.meta.pls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE,study=exp,scale=TRUE)

object=res
newdata=data
study.test=exp



## ======================================================================
# from sgcca
## ======================================================================
source("mixOmics/R/predict.meta.block.pls.R")
object=wrapper.pls(X=data.light,Y=Y.mat.light,ncomp=3)
object=wrapper.meta.pls(X=data.light,Y=Y.mat.light,ncomp=3,study=c(rep("Flo",20),rep("haha",31)),scale=TRUE)

source("mixOmics/R/predict.meta.block.pls.R")
object2=wrapper.meta.pls(X=data.light[c(1:13),],Y=Y.mat.light[c(1:13),],ncomp=3,study=c(rep("Flo",13)),scale=TRUE)
pred2=predict.meta.block.pls(object2,newdata=data.light[c(1:13),],study.test=factor(object2$study))

object=wrapper.meta.pls(X=rbind(data.light[c(1:13),],data.light[c(1:13),]),Y=rbind(Y.mat.light[c(1:13),],Y.mat.light[c(1:13),]),ncomp=3,study=c(rep("Flo",13),rep("haha",13)),scale=TRUE)

newdata=rbind(data.light[c(1:13),],data.light[c(1:13),])
rownames(newdata)=NULL
pred=predict.meta.block.pls(object,newdata=newdata,study.test=factor(object$study))

newdata=rbind(data.light[c(1:13),],data.light[c(1:13),])
study.test=factor(object$study)




object3=wrapper.meta.plsda(X=data.light[c(1:13),],Y=type.id.light[c(1:13)],ncomp=3,study=c(rep("Flo",13)),scale=TRUE)
pred3=predict.meta.block.pls(object3,newdata=data.light[c(1:13),],study.test=factor(object3$study))


source("mixOmics/R/predict.R")
object4=plsda(X=data.light[c(1:13),],Y=type.id.light[c(1:13)],ncomp=3,near.zero.var=FALSE)
pred4=predict(object4,newdata=data.light[c(1:13),])


#pred3$Y.hat[[1]][,,1] and pred4$predict[,,1] should be the same
all.equal(pred3$Y.hat[[1]][,,1],pred4$predict[,,1])
#it's not// pretty close though. all good because pls is not converging anyway!


object=wrapper.meta.pls(X=data.light[c(1:13),],Y=Y.mat.light[c(1:13),],ncomp=3,study=c(rep("Flo",13)),scale=TRUE)
newdata=list(data.light[c(1:13),])




## =========================================================================================================
## =========================================================================================================
## ===========================                  test_wrappers               ==============================##
## =========================================================================================================
## =========================================================================================================


#wraper.meta.spls.hybrid, function not available to user, so not a proper check to conduct
res=wrapper.meta.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,tol=1e-25)
res=wrapper.meta.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,5),ncomp=3,near.zero.var=FALSE,tol=1e-25)






## ======================================================================
## ==      data for one study
## ======================================================================



#wraper.pls
res=wrapper.pls(X=data.light,Y=Y.mat.light,ncomp=3)
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
pred=predict(res,newdata=A.light)
res=wrapper.block.pls(list(data),Y=Y.mat,ncomp=2)
pred=predict(res,newdata=list(data))

#wraper.block.spls
res=wrapper.block.spls(X=A.light,indY=2,keepX=list(block1=c(10,5,15),block2=c(3,2)),ncomp=c(3,3))
pred=predict(res,newdata=A.light)
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

#wraper.meta.pls
res=wrapper.meta.pls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE,study=exp)
pred=predict(res,newdata=data,study.test=exp)

#wraper.meta.plsda
res=wrapper.meta.plsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,study=exp)
pred=predict(res,newdata=data,study.test=exp)

#wraper.meta.spls
res=wrapper.meta.spls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)
pred=predict(res,newdata=data,study.test=exp)

#wraper.meta.splsda
res=wrapper.meta.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)
pred=predict(res,newdata=data,study.test=exp)



#wraper.meta.block.pls
res=wrapper.meta.block.pls(X=list(X=data,Y=Y.mat),indY=2,ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.meta.block.pls(list(data),Y=Y.mat,ncomp=2)
pred=predict(res,newdata=list(X=data))

#wraper.meta.block.spls
res=wrapper.meta.block.spls(X=list(X=data,Y=Y.mat),indY=2,keepX=list(block1=c(10,5)),ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.meta.block.spls(list(data),Y=Y.mat,ncomp=2)
pred=predict(res,newdata=list(X=data))

#wraper.meta.block.plsda
res=wrapper.meta.block.plsda(X=list(X=data,Y=type.id),indY=2,ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.meta.block.plsda(list(data),Y=type.id,ncomp=2)
pred=predict(res,newdata=list(X=data))

#wraper.meta.block.splsda
res=wrapper.meta.block.splsda(X=list(X=data,Y=type.id),indY=2,keepX=list(block1=c(10,5)),ncomp=c(2,2))
pred=predict(res,newdata=list(X=data))
res=wrapper.meta.block.splsda(list(data),Y=type.id,ncomp=2)
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

#meta.pls
res=mixOmics(data,Y=Y.mat,ncomp=3,study=exp)
pred=predict(res,newdata=data,study=exp)

#meta.spls
res=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)
res1=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)#with gene names in keepX.constraint
res2=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),keepX=c(100,50),ncomp=3)#with numbers in keepX.constraint
all.equal(res1,res2)

pred=predict(res,newdata=data,study=exp)
pred=predict(res1,newdata=data,study=exp)
pred=predict(res2,newdata=data,study=exp)

res=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),ncomp=3)#meta.spls missing keepX is completed by pls-like
res$keepX
pred=predict(res,newdata=data,study=exp)

#meta.plsda
res=mixOmics(data,type.id,ncomp=3,study=exp)
pred=predict(res,newdata=data,study=exp)

#meta.splsda
res=mixOmics(X=data,Y=type.id,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,15),ncomp=3)
pred=predict(res,newdata=data,study=exp)



source("mixOmics/R/predict.meta.block.pls.R")

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
res=mixOmics(X=A,Y=type.id,study=exp)
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

res=mixOmics(X=A,tau="optimal")
res=mixOmics(X=A,indY=2,tau=c(1,1),ncomp=c(3,1)) #OK

#sparse RGCCA
# keep variable 5 and 10 on the 2first comp of block1, keep variable 2 on comp1 for block 2. completed by keeping all variable on comp2 for block2
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX.constraint=list(c(5,10),2),ncomp=c(2,2))
res1=mixOmics(X=A,tau=c(1,1),keepX.constraint=list(list(5,10),list(2)),ncomp=c(2,2))#same

# keep variable 5 on comp1 and 10 on comp2 of block1; keep variable 2 on comp1 for block 2. completed by keeping all variable on comp3 for block1 and comp2/3 for block2
res2=mixOmics(X=A,tau=c(1,1),keepX.constraint=list(c(5,10),2),ncomp=c(3,3))

# keep variable 5 and 10 on the first comp of block1, variable 3 on the second comp of block1; keep variable 2 on comp1 for block 2. completed by keeping all variable on comp3 for block1 and comp2/3 for block2
res3=mixOmics(X=A,tau=c(1,1),keepX.constraint=list(list(c(5,10),3),list(2)),ncomp=c(3,3))


all.equal(res,res3)

which(res$loadings$X[,1]!=0)
which(res$loadings$X[,2]!=0)

which(res1$loadings$X[,1]!=0)
which(res1$loadings$X[,2]!=0)

which(res2$loadings$X[,1]!=0)
which(res2$loadings$X[,2]!=0)
sum(res2$loadings$X[,3]!=0)

which(res3$loadings$X[,1]!=0)
which(res3$loadings$X[,2]!=0)
sum(res3$loadings$X[,3]!=0)

#OK keep variable 5 and 10 on the first comp of block1, variable 3 on the second comp of block1. completed by keeping all variable on comp3 for block1 and comp1/2/3 for block2
res=mixOmics(X=A,tau=c(1,1),keepX.constraint=list(list(c(5,10),3)),ncomp=c(2,2)) #OK

#OK keep variable 5 and 10 on the first comp of block1, variable 3 on the second comp of block1. keep the 10 most important variables on comp3 of block1. completed by keeping all variable on comp3 for block1 and comp1/2/3 for block2
res=mixOmics(X=A,tau=c(1,1),keepX.constraint=list(list(c(5,10),3)),keepX=list(10),ncomp=c(3,3)) #OK


res=mixOmics(X=A,tau=c(1,1),keepX.constraint=list(),keepX=list(c(10,10)),ncomp=c(2,2))#OK
res=mixOmics(X=A,tau=c(1,1),keepX.constraint=list(),keepX=list(c(10,10)),ncomp=c(2,1))#OK
res=mixOmics(X=A,tau=c(1,1),keepX=list(c(10,10)),ncomp=c(3,3))#OK
res=mixOmics(X=A,tau=c(1,1),keepX=list(5,2),keepX.constraint=list(list(c(5,10))),ncomp=c(2,2)) #OK




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



#RGCCA

#check RGCCA
res=mixOmics(X=A,tau=1) #error length tau
res=mixOmics(X=A,tau=1,mode="regression") # error length tau (and mode)
res=mixOmics(X=A,tau=c(1,1),mode="regression") # error mode
res=mixOmics(X=A,tau=c(1,1),init="svd") #error svd
res=mixOmics(X=A,tau=c(1,1),scheme="bla") #error scheme
res=mixOmics(X=A,tau=c(1,1),study=rep(1,nrow(data))) #warnings study
res=mixOmics(X=A,tau=c(1,1),study=factor(rep(1,nrow(data)))) #warnings study
res=mixOmics(X=list(X=data,Y=type.id),tau=1) #A[[2]] must be numeric matrix




#sparse.RGCCA
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX=list(list(1:5),list(1:10)),ncomp=c(1,1)) #error keepX
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX=list(c(5,10)),ncomp=c(1,1)) #error keepX (keepX=c(5,10) for block1, but ncomp=1
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX=list(5,10),ncomp=c(1,1)) #error keepX, keepY=10 for Y, should be less than 3
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX.constraint=list(list(c(5,10),3)),ncomp=c(1,1)) #error keepX.constraint. Keep variables 5,10 on comp1 and 3 on comp2, but ncomp=1
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX.constraint=list(list(c(5,10),3)),keepX=list(10),ncomp=c(2,2)) #should give an  error
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX=list(5,2),keepX.constraint=list(c(5,10)),ncomp=c(2,2)) #error keepX+keepX.constraint !=ncomp
res=mixOmics(X=A,tau=c(1,1),mode="canonical",keepX=list(5,2),keepX.constraint=list(list(c(5,10))),ncomp=c(1,2)) #error keepX+keepX.constraint !=ncomp



