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

tune= tune.splsda(X,Y,ncomp=1,nrepeat=1,logratio="none",test.keepX = c(5),folds=10,dist="max.dist",already.tested.X=NULL, progressBar = FALSE)



tune= tune.splsda(X,Y,ncomp=1,nrepeat=10,logratio="none",test.keepX = c(5, 10, 15),folds=10,dist="max.dist",already.tested.X=NULL, progressBar = TRUE, light.output=FALSE)


tune= tune.splsda(X,Y,ncomp=3,nrepeat=10,logratio="none",test.keepX = c(5, 10, 15),folds=10,dist="max.dist",already.tested.X=NULL, progressBar = FALSE)

#source("mixOmics/R/tune.splsda.R")
#source("mixOmics/R/MCVfold.R")



data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample,
stimul = vac18$stimulation)

# multilevel sPLS-DA model
res.1level = splsda(X, ncomp = 3, multilevel = design,
keepX = c(30, 137, 123))

tune= tune.splsda(X,ncomp=3,nrepeat=10,logratio="none",test.keepX = c(5,50,100),folds=10,dist="max.dist",already.tested.X=NULL, progressBar = FALSE, multilevel = design)



# mint.splsda
# ----
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

tt=tune.mint.splsda(X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,study=exp,test.keepX=seq(1,10,1))


tt=tune(method="mint.splsda",X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,study=exp,test.keepX=seq(1,10,1))

# with the full data, to check that we get the same results: 2+15
#tt=tune.mint.splsda(X=data.learn,Y=type.id.learn,ncomp=3,near.zero.var=FALSE,study=experiment.id.learn,test.keepX=seq(1,20,1))


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
