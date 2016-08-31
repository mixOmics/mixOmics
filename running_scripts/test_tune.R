# created on 12/03/15
# last modified: 02-03-2016
# Author: F.Rohart
#purpose: test the pls/plsda/spls/splsda function


#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)


# splsda
# ----
data(breast.tumors)
X <- breast.tumors$gene.exp
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(breast.tumors$sample$treatment)

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))

tune= tune.splsda(X,Y,ncomp=1,nrepeat=1,logratio="none",test.keepX = c(5,10),folds=10,dist="max.dist", progressBar = FALSE)

tune= tune.splsda(X,Y,ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5, 10, 15),folds=10,dist="max.dist", progressBar = FALSE, already.tested.X = c(5,10))

tune= tune.splsda(X,Y,ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5, 10, 15),folds=10,dist="max.dist", progressBar = FALSE, already.tested.X = list(comp1=c(5,10)), constraint=TRUE)


tune= tune.splsda(X,Y,ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5, 10, 15),folds=10,dist="max.dist", progressBar = FALSE, light.output=FALSE)


tune= tune.splsda(X,Y,ncomp=2,nrepeat=1,logratio="none",test.keepX = c(5, 15),folds=10,dist="max.dist", progressBar = FALSE)

tune= tune(method="splsda",X,Y,ncomp=2,nrepeat=1,logratio="none",test.keepX = c(5, 15),folds=10,dist="max.dist", progressBar = FALSE)

plot(tune)

#source("mixOmics/R/tune.splsda.R")
#source("mixOmics/R/MCVfold.R")



data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample,
stimul = vac18$stimulation)

# multilevel sPLS-DA model
res.1level = splsda(X, Y = design[,2], ncomp = 3, multilevel = design[,1],
keepX = c(30, 137, 123))

tune= tune.splsda(X,Y=design[,2],ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5,50,100),folds=10,dist="max.dist", progressBar = FALSE, multilevel = design[,1])

tune= tune.splsda(X,Y=design[,2],ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5,10,15),folds=10,dist="max.dist", progressBar = FALSE, multilevel = design[,1])

tune= tune.splsda(X,Y=design[,2],ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5,10,15),folds=10,dist="max.dist", progressBar = FALSE, multilevel = design[,1],constraint=FALSE)
tune$choice.keepX
tune$choice.keepX.constraint

tune= tune.splsda(X,Y=design[,2],ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5,10,15),folds=10,dist="max.dist", progressBar = FALSE, multilevel = design[,1],constraint=TRUE)
tune$choice.keepX
tune$choice.keepX.constraint

tune= tune.splsda(X,Y=design[,2],ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5,10,15),folds=10,dist="max.dist",already.tested.X=c(5,10), progressBar = FALSE, multilevel = design[,1],constraint=FALSE)
tune$choice.keepX
tune$choice.keepX.constraint


# justified error
#tune= tune.splsda(X,Y=design[,2],ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5,10,15),folds=10,dist="max.dist",already.tested.X=c(5,10), progressBar = FALSE, multilevel = design[,1],constraint=TRUE)

tune= tune.splsda(X,Y=design[,2],ncomp=3,nrepeat=5,logratio="none",test.keepX = c(5,10,15),folds=10,dist="max.dist",already.tested.X=list(compA=c(1,3,6)), progressBar = FALSE, multilevel = design[,1],constraint=TRUE)

plot(tune)

# mint.splsda
# ----
library(mixOmics)
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

res=mint.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)
out=perf(res)

plot(out)
plot(out,study="all.partial")
plot(out,study="global", overlay="measure")
plot(out,study="1", overlay="measure")


tt=tune.mint.splsda(X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,study=exp,test.keepX=seq(1,100,10), progressBar = FALSE,constraint=FALSE)


tt=tune(method="mint.splsda",X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,study=exp,test.keepX=seq(1,10,1), progressBar = FALSE)

plot(tt)


if(FALSE)
{
# with the full data, to check that we get the same results: 2+15
load("/Users/florian/Work/git/multi-group/Fibro-ESC-iPSC.15exp.342samples.all.set.Rdata")

experiment.id=factor(experiment.id)
exp.learn=c("Briggs","Bock","Guenther","Takahashi","Ebert","Chung",#"Kim",
"Maherali","Marchetto")
exp.test=c("Kim","Andrade","Hu","Yu","Vitale","Si-Tayeb","Loewer")#levels(experiment.id)[-which(levels(experiment.id)%in%exp.learn)]
ind.learn=which(experiment.id%in%exp.learn)
ind.test=which(experiment.id%in%exp.test)

data.learn=data[ind.learn,]
experiment.id.learn=factor(experiment.id[ind.learn])
type.id.learn=factor(type.id[ind.learn])

data.test=data[ind.test,]
experiment.id.test=factor(experiment.id[ind.test])
type.id.test=factor(type.id[ind.test])

levels(experiment.id.learn);levels(experiment.id.test)


tt=tune.mint.splsda(X=data.learn,Y=type.id.learn,ncomp=3,near.zero.var=FALSE,study=experiment.id.learn,test.keepX=seq(1,20,1))

res=mint.splsda(X=data.learn,Y=type.id.learn,ncomp=2,near.zero.var=FALSE,
keepX.constraint=tt$choice.keepX.constraint[1:2],keepX=NULL,study=experiment.id.learn)

out=perf(res)
}
# rcc
# ----
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.res <- rcc(X, Y)

## Regularized CCA
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res1 <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)


tt=tune.rcc(X, Y)

par(opar)
