###############################
############# PLS #############
###############################
#### Difference because of not the same starting point between the 2 algo "pls" and "spls" in mixOmics package when tol = 1e-06
#source("/Users/florian/Work/git/package-mixomics/test_scripts/mint.block_test/chack_mixOmics.R")
rm(list=ls())
setwd("/Users/florian/Work/git/package-mixOmics/")
library(mixOmicsv6)


load("test_scripts/mint.block_test/Fibro-ESC-iPSC.6exp.167samples.light.Rdata") #load data, type.id, exp among others

if(FALSE)
{
    sourceDir <- function(path, trace = TRUE, ...) {
        for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
            if(trace) cat(nm,":")
            source(file.path(path, nm), ...)
            if(trace) cat("\n")
        }
    }
    sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE)
    source("/Users/florian/Work/git/package-mixOmics/mixOmics/R/predict.mint.block.pls.R")

source("test_scripts/mint.block_test/test_all.R")
}



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


#add useless column for near.zero.var=TRUE tests

data.near.zero=cbind(data,matrix(0,nrow=nrow(data),ncol=10))
colnames(data.near.zero)=c(colnames(data),paste0("nearzero",1:10))




## =========================================================================================================
## =========================================================================================================
## ===========================                  test_wrappers               ==============================##
## =========================================================================================================
## =========================================================================================================


#wraper.mint.spls.hybrid, function not available to user, so not a proper check to conduct
#res=wrapper.mint.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,tol=1e-25)
#res=wrapper.mint.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,5),ncomp=3,near.zero.var=FALSE,tol=1e-25)






## ======================================================================
## ==      data for one study
## ======================================================================



#######  wraper.pls
res=pls(X=data.light,Y=Y.mat.light)
res=pls(X=data.light,Y=Y.mat.light,ncomp=3)
res=pls(X=data.near.zero,Y=Y.mat,ncomp=3)
res=pls(X=data.near.zero,Y=Y.mat,near.zero.var=TRUE)



#######  wraper.spls
res=spls(X=data.light,Y=Y.mat.light)
res=spls(X=data.light,Y=Y.mat.light,ncomp=3,keepX=c(10,5,15))
res=spls(X=data.light,Y=Y.mat.light,ncomp=3,keepX=c(10,5)) #complete keepX
res=spls(X=data.light,Y=Y.mat.light,ncomp=3,keepX=c(10),keepX.constraint=list(comp1=c(100,1,3),comp2=c(10)))
res=spls(X=data.light,Y=Y.mat.light,ncomp=3,keepX=c(1),keepX.constraint=list(comp1=c(100),comp2=c(10)))
res=spls(X=data.light,Y=Y.mat.light,ncomp=3,keepX.constraint=list(comp1=c(100),comp2=c(10)))
res=spls(X=data.light,Y=Y.mat.light,ncomp=3,keepY=c(3),keepY.constraint=list(comp1=c(1),comp2=c(2)))

#######  wraper.plsda
res=plsda(X=data.light,Y=type.id.light,ncomp=3)


#######  wraper.splsda
res=splsda(X=data.light,Y=type.id.light)
res=splsda(X=data.light,Y=type.id.light,ncomp=3,keepX=c(10,5,15))
res=splsda(X=data.light,Y=type.id.light,ncomp=3,keepX=c(10,5)) #complete keepX
res=splsda(X=data.light,Y=type.id.light,ncomp=3,keepX=c(10),keepX.constraint=list(comp1=c(100,1,3),comp2=c(10)))
res=splsda(X=data.light,Y=type.id.light,ncomp=3,keepX=c(1),keepX.constraint=list(comp1=c(100),comp2=c(10)))
res=splsda(X=data.light,Y=type.id.light,ncomp=3,keepX.constraint=list(comp1=c(100),comp2=c(10)))



#######  wraper.block.pls
res=block.pls(X=A.light,indY=2)
res=block.pls(list(data),Y=Y.mat,ncomp=2)
res2=block.pls(list(data),Y=Y.mat,ncomp=2,bias=TRUE)
res3=block.pls(list(data),Y=Y.mat,ncomp=2,scale=FALSE)
all.equal(res,res2)
all.equal(res,res3)
#res=block.pls(list(data.light),Y=Y.mat.light,ncomp=2,verbose=TRUE,tol=1e-30)
res=block.pls(list(data.near.zero),Y=Y.mat,ncomp=2,near.zero.var=TRUE)




#######  wraper.block.spls
res=block.spls(X=A.light,indY=2,keepX=list(block1=c(10,5,15),block2=c(3,2)),ncomp=c(3,3))
res=block.spls(list(data),Y=Y.mat,keepX=list(block1=c(10,5,15),block2=c(3,2)),ncomp=3)



#######  wraper.block.plsda
res=block.plsda(X=list(X=data,Y=type.id),indY=2,ncomp=c(2,2))
res=block.plsda(list(data),Y=type.id,ncomp=2)



#######  wraper.block.splsda
res=block.splsda(X=list(X=data,Y=type.id),keepX=list(block1=c(10,5)),indY=2,ncomp=c(3,2))
res=block.splsda(X=list(X=data,Y=type.id),keepX=list(block1=c(10,5)),indY=2,ncomp=c(2,2))
res=block.splsda(X=list(X=data),Y=type.id,ncomp=3,keepX=list(c(100)),
keepX.constraint=list(list(comp1=c("ENSG00000001084","ENSG00000001461"),comp2=c("ENSG00000000938"))))





# same results for sgccda and block.splsda. outputs are different though
res=wrapper.sgccda(X=data,Y=type.id,keepA=list(c(10,5,10)),ncomp=3)
res2=block.splsda(X=list(X=data),Y=type.id,keepX=list(c(10,5,10)),ncomp=3,mode="canonical")



## ======================================================================
## ==      data for several study
## ======================================================================


#######  wraper.mint.pls
res=mint.pls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE,study=exp)
res=mint.pls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE,study=as.character(exp))
res=mint.pls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE)


#######  wraper.mint.plsda
res=mint.plsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,study=exp)


#######  wraper.mint.spls
res=mint.spls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)

#######  wraper.mint.splsda
res=mint.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)



#######  wraper.mint.block.pls
res=mint.block.pls(X=list(X=data,Y=Y.mat),indY=2,ncomp=c(2,2))
res=mint.block.pls(list(data),Y=Y.mat,ncomp=2)

#######  wraper.mint.block.spls
res=mint.block.spls(X=list(X=data,Y=Y.mat),indY=2,keepX=list(block1=c(10,5)),ncomp=c(2,2))
res=mint.block.spls(list(data),Y=Y.mat,ncomp=2)

#######  wraper.mint.block.plsda
res=mint.block.plsda(X=list(X=data,Y=type.id),indY=2,ncomp=c(2,2))
res=mint.block.plsda(list(data),Y=type.id,ncomp=2)

#######  wraper.mint.block.splsda
res=mint.block.splsda(X=list(X=data,Y=type.id),indY=2,keepX=list(block1=c(10,5)),ncomp=c(2,2))
res=mint.block.splsda(list(data),Y=type.id,ncomp=2)





## =========================================================================================================
## =========================================================================================================
## ===========================                  test_mixOmics               ==============================##
## =========================================================================================================
## =========================================================================================================


type.id.light=factor(type.id.light)
type.id=factor(type.id)

## ======================================================================
## ==      data for one study
## ======================================================================

res=mixOmics(data.light,Y.mat.light,ncomp=3) #pls
res=mixOmics(data.light,Y.mat.light,ncomp=3,keepX=c(10,5,15)) #spls
res=mixOmics(data.light,type.id.light,ncomp=3) #plsda
res=mixOmics(data.light,type.id.light,ncomp=3,keepX=c(10,5,15))#splsda


#block.pls
res=mixOmics(X=list(data),Y=unmap(type.id))
res=mixOmics(X=A,indY=2,ncomp=c(2,2))
res=mixOmics(X=A,indY=2)


#block.spls
res=mixOmics(X=list(data),Y=unmap(type.id),keepX=list(c(10,5,15)),ncomp=3)
res=mixOmics(X=list(data=data,Y=Y.mat),indY=2,keepX=list(c(10,5,15)),ncomp=c(3,3))
res=mixOmics(X=A,indY=2,keepX.constraint=list(X=list(1:10)),ncomp=c(2,2))

# block.plsda
res=mixOmics(X=A,Y=type.id)
res=mixOmics(X=list(data=data,Y=type.id),indY=2,ncomp=c(3,3))
res=mixOmics(X=A,Y=type.id)
res=mixOmics(X=list(data),Y=type.id)

# block.splsda
res=mixOmics(X=list(data=data,Y=type.id),indY=2,keepX=list(c(10,5,15)),ncomp=c(3,3))


system.time(splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))
system.time(mixOmics(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=TRUE))

system.time(splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE))
system.time(splsda(X=data,Y=type.id,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,scale=TRUE))
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

#mint.spls
res=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)
res1=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(100,50),ncomp=3)#with gene names in keepX.constraint
res2=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),keepX=c(100,50),ncomp=3)#with numbers in keepX.constraint
all.equal(res1,res2)
res=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),ncomp=3)#mint.spls missing keepX is completed by pls-like
res$keepX

#mint.plsda
res=mixOmics(data,type.id,ncomp=3,study=exp)

#mint.splsda
res=mixOmics(X=data,Y=type.id,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,15),ncomp=3)




#block.pls
res=mixOmics(X=list(data),Y=unmap(type.id),study=exp)

#block.spls
res=mixOmics(X=list(data),Y=unmap(type.id),keepX=list(10,5,15),study=exp)
res=mixOmics(X=list(data=data,Y=Y.mat),indY=2,keepX=list(c(10,5,15)),ncomp=c(3,3),study=exp)
res=mixOmics(X=A,indY=2,keepX.constraint=list(X=list(1:10)),ncomp=c(3,1)) # OK

# block.plsda
res=mixOmics(X=A,Y=type.id,study=exp)

#block.plsda
res=mixOmics(X=list(data=data,Y=type.id),indY=2,ncomp=c(3,3),study=exp)

#block.splsda
res=mixOmics(X=list(data=data,Y=type.id),indY=2,keepX=list(c(10,5,15)),ncomp=c(3,3),study=exp)



## ======================================================================
## ==      RGCCA/sparse.RGCCA
## ======================================================================

#RGCCA
res=mixOmics(X=A,tau=c(1,1))
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




