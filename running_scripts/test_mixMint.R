## =========================================================================================================
## =========================================================================================================
## ===========================                  test mixMint module               ==============================##
## =========================================================================================================
## =========================================================================================================

opar <- par(no.readonly = TRUE)

data(stemcells) #load gene, celltype and study

data=stemcells$gene
type.id=stemcells$celltype
exp=stemcells$study

Y.mat=unmap(stemcells$celltype)
rownames(Y.mat)=rownames(stemcells$gene)
ind=which(stemcells$study=="1")
study.light=stemcells$study[ind]
type.id.light=stemcells$celltype[ind]
data.light=stemcells$gene[ind,]
Y.mat.light=unmap(type.id.light)
rownames(Y.mat.light)=rownames(data.light)

A=list(X=data,Y=Y.mat)
A.light=list(X=data.light,Y=Y.mat.light)


#add useless column for near.zero.var=TRUE tests

data.near.zero=cbind(stemcells$gene,matrix(0,nrow=nrow(stemcells$gene),ncol=10))
colnames(data.near.zero)=c(colnames(stemcells$gene),paste0("nearzero",1:10))


res=mint.spls(X=data,Y=Y.mat,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)

#source("/Users/florian/Work/git/package-mixomics/mixOmics/R/plotIndiv.mint.R")
#source("/Users/florian/Work/git/package-mixomics/mixOmics/R/check.plotIndiv.R")
#source("/Users/florian/Work/git/package-mixomics/mixOmics/R/internal_graphicModule.R")

plotIndiv(res,study=c(1:2,"all",6:5),group=type.id)
plotIndiv(res,study=c(1:2,"all",6:5),group=type.id,add.legend=T)
plotIndiv(res,study=c("all"),group=type.id,add.legend=T)
plotIndiv(res,study=c("1"),group=type.id,add.legend=T)
plotIndiv(res,study=c(1,2),group=type.id,add.legend=T,main="bla")

plotIndiv(res,study=c(1:2),group=type.id,add.legend=T,legend.position="right",legend.title.size=rel(2),legend.text.size=rel(1.7),size.xlabel=rel(1.2),size.ylabel=rel(1.2),size.axis=rel(1.2),size.subtitle=rel(3),size.title=rel(4),main="bla")

plotIndiv(res,study=c(1:2),group=type.id,add.legend=T,legend.position="bottom",legend.title.size=rel(2),legend.text.size=rel(1.7),size.xlabel=rel(1.2),size.ylabel=rel(1.2),size.axis=rel(1.2),size.subtitle=rel(3),size.title=rel(4),main="bla")

plotIndiv(res,study=c(1,2),group=type.id,add.legend=T,style="lattice")

plotIndiv(res,study=c(1:2,"all",6:5),group=type.id,add.legend=T,pch=14+as.numeric(factor(exp)))


plotIndiv(res,study=c(1:2,"all",6:5),group=type.id,subtitle=paste0("study:",letters[seq( from = 1, to = 5 )]))
plotIndiv(res,study=c(1:2,"all",6:5),group=type.id,subtitle=letters[seq( from = 1, to = 5 )],main="MINT")



#wraper.mint.spls.hybrid, function not available to user, so not a proper check to conduct
#res=wrapper.mint.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX=c(10,5,10),ncomp=3,near.zero.var=FALSE,tol=1e-25)
#res=wrapper.mint.spls.hybrid(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000006576","ENSG00000008226")),keepX=c(10,5),ncomp=3,near.zero.var=FALSE,tol=1e-25)






## ======================================================================
## ==      data for one study
## ======================================================================




#######  wraper.block.pls
res=block.pls(list(data),Y=Y.mat,ncomp=2)
res2=block.pls(list(data),Y=Y.mat,ncomp=2,bias=TRUE)
res3=block.pls(list(data),Y=Y.mat,ncomp=2,scale=FALSE)

selectVar(res,block=2)

#######  wraper.block.spls
res=block.spls(X=A.light,indY=2,keepX=list(block1=c(10,5,15),block2=c(3,2)),ncomp=c(3,3))
res=block.spls(list(data),Y=Y.mat,keepX=list(block1=c(10,5,15),block2=c(3,2)),ncomp=3)

selectVar(res)
selectVar(res,block=2,comp=3)


#######  wraper.block.plsda
res=block.plsda(X=list(X=data,Y=type.id),indY=2,ncomp=c(2,2))
res=block.plsda(list(data),Y=type.id,ncomp=2)

selectVar(res,block=2)


#######  wraper.block.splsda
res=block.splsda(X=list(X=data,Y=type.id),keepX=list(block1=c(10,5)),indY=2,ncomp=c(3,2))
res=block.splsda(X=list(X=data,Y=type.id),keepX=list(block1=c(10,5)),indY=2,ncomp=c(2,2))
res=block.splsda(X=list(X=data),Y=type.id,ncomp=3,keepX=list(c(100)),
keepX.constraint=list(list(comp1=c("ENSG00000187243","ENSG00000175573"),comp2=c("ENSG00000133475"))))

res2=block.splsda(X=list(X=data),Y=type.id,keepX=list(c(10,5,10)),ncomp=3,mode="canonical")


selectVar(res2)
selectVar(res2,block=2,comp=2)




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

selectVar(res)
selectVar(res,block=2,comp=3)
plotVar(res)


#source("mixOmics/R/plotIndiv.pls.R")


#source("mixOmics/R/plotIndiv.pls.R")
#source("mixOmics/R/check.plotIndiv.R")
#source("mixOmics/R/internal_graphicModule.R")
#source("mixOmics/R/plotIndiv.mint.R")

#plotIndiv(res,study=1)
plotIndiv(res,study=c(1:2,"all",5:6))
plotIndiv(res,study=c(1:2,"all",6:5),group=type.id)
plotIndiv(res,study=c(1:2,"all",6:5),group=type.id,subtitle=paste0("study:",letters[seq( from = 1, to = 5 )]))
plotIndiv(res,study=c(1:2,"all",6:5),group=type.id,subtitle=letters[seq( from = 1, to = 5 )],main="MINT")


plotIndiv(res,study=c(1:2,"all",5:6),style="lattice")
plotIndiv(res,study=c(1:2,"all",6:5),group=type.id,style="lattice")

plotIndiv(res,study=c(1:2,"all",5:6),style="graphics")
plotIndiv(res,study=c(1:2,"all",6:5),group=type.id,style="graphics")


plotIndiv(res,ind.names=FALSE)

plotIndiv(res,ind.names=FALSE,group=type.id)

plotIndiv(res,ind.names=FALSE,group=type.id,add.legend=TRUE)



#######  wraper.mint.splsda
res=mint.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)

selectVar(res)
selectVar(res,block=2,comp=3)


#######  wraper.mint.block.pls
res=mint.block.pls(X=list(X=data,Y=Y.mat),indY=2,ncomp=c(2,2))
res=mint.block.pls(list(data),Y=Y.mat,ncomp=2)

#######  wraper.mint.block.spls
res=mint.block.spls(X=list(X=data,Y=Y.mat),indY=2,keepX=list(block1=c(10,5)),ncomp=c(2,2))
res=mint.block.spls(list(data),Y=Y.mat,ncomp=2)

selectVar(res)
selectVar(res,block=2,comp=2)
plotVar(res)


#######  wraper.mint.block.plsda
res=mint.block.plsda(X=list(X=data,Y=type.id),indY=2,ncomp=c(2,2))
res=mint.block.plsda(list(data),Y=type.id,ncomp=2)

#######  wraper.mint.block.splsda
res=mint.block.splsda(list(data),Y=type.id,ncomp=2)
res=mint.block.splsda(X=list(X=data,Y=type.id),indY=2,keepX=list(block1=c(10,5)),ncomp=c(2,2))

selectVar(res)
selectVar(res,block=2,comp=2)
plotVar(res)



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
res1=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c("ENSG00000187243","ENSG00000175573")),keepX=c(100,50),ncomp=3)#with gene names in keepX.constraint
res2=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),keepX=c(100,50),ncomp=3)#with numbers in keepX.constraint
all.equal(res,res1)
res=mixOmics(X=data,Y=Y.mat,study=exp,keepX.constraint=list(c(120,179)),ncomp=3)#mint.spls missing keepX is completed by pls-like
res$keepX

#mint.plsda
res=mixOmics(data,type.id,ncomp=3,study=exp)

#mint.splsda
res=mixOmics(X=data,Y=type.id,study=exp,keepX.constraint=list(c("ENSG00000187243","ENSG00000175573")),keepX=c(10,15),ncomp=3)




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


par(opar)


