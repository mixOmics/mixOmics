# created on 25/08/16
# Author: F.BARTOLO
#purpose: test the auroc function


#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)


############################AUROC


# plsda/splsda
# ----

#library(mixOmics)
data(breast.tumors)
set.seed(1)
test=sample(1:47,5,replace=FALSE)
X <- breast.tumors$gene.exp
X.test<-breast.tumors$gene.exp[test,]
Y <- breast.tumors$sample$treatment
Y.test<-breast.tumors$sample$treatment[test]

res.plsda <- plsda(X, Y, ncomp = 2)

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)


auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.comp = 2)


res.plsda <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))

res.plsda <- splsda(X, Y, ncomp = 2, keepX = c(2, 2))

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)


auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.comp = 2)




data(liver.toxicity)
set.seed(1)
test=sample(1:64,7,replace=FALSE)
X <- liver.toxicity$gene
X.test<-liver.toxicity$gene[test,]
Y <- liver.toxicity$treatment[, 4]
Y.test <- liver.toxicity$treatment[test, 4]

res.plsda <- plsda(X, Y, ncomp = 2)

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)


auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.comp = 2)



res.plsda <- plsda(X, Y, ncomp = 3)

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)


auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.comp = 2)




# mint.plsda
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

res.plsda=mint.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = TRUE,roc.comp = 3)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 3)

# ###############WARNINGS / RESOLU

auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 3)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 3)





res.plsda=block.plsda(X=list(X=data,Y=type.id),indY=2,ncomp = 2)

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)




res.plsda=block.plsda(list(XX=data),Y=type.id,ncomp=2)

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)



res.plsda=block.splsda(X=list(X=data,Y=type.id),keepX=list(X=c(10,5)),indY=2,ncomp = 2)
auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)


res.plsda=block.splsda(X=list(X=data),Y=type.id,ncomp=3,keepX=list(X=c(100)))#,
#keepX.constraint=list(X=list(comp1=c("ENSG00000164930","ENSG00000044090"),comp2=c("ENSG00000109819"))))

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)



res.plsda=block.splsda(X=list(X=data),Y=type.id,ncomp=3,keepX=list(X=c(100)))#,
#keepX.constraint=list(X=list(comp1=c("ENSG00000164930","ENSG00000044090"),comp2=c("ENSG00000109819"))))

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)


res.plsda=block.splsda(X=list(X2=data),Y=type.id,keepX=list(X2=c(10,5,10)),ncomp=3,mode="canonical")

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(X2=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(X2=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(X2=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(X2=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)



res.plsda=mint.plsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,study=exp)

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = TRUE,roc.comp = 3)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 3)

# ###############WARNINGS / RESOLU

auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 3)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 3)



res.plsda=mint.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = TRUE,roc.comp = 3)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 3)

# ###############WARNINGS/RESOLU

auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 3)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)
auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 3)



res.plsda=mint.block.plsda(X=list(X=data,Y=type.id),indY=2,ncomp = 2)


####PROBLEME ERROR//RESOLU SORTIE res.plsda$ind.mat est NULL  /// RESOLU

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(X=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)


res.plsda=mint.block.plsda(list(XX=data),Y=type.id,ncomp=2)

####PROBLEME ERROR//RESOLU//

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)



res.plsda=mint.block.splsda(list(XX=data),Y=type.id,ncomp=2)

####PROBLEME ERROR//RESOLU

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)



res.plsda=mint.block.splsda(X=list(XX=data,Y=type.id),indY=2,keepX=list(XX=c(10,5)),ncomp = 2)

####PROBLEME ERROR//RESOLU

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)




# sgccda
# ----

set.seed(12)
data(nutrimouse)
train=sample(1:40,40,replace=FALSE)
test=sample(1:40,4,replace=FALSE)
Y = nutrimouse$diet[train]
Y.test=nutrimouse$diet[test]
data = list(gene = nutrimouse$gene[train,], lipid = nutrimouse$lipid[train,])
data.test=list(gene = nutrimouse$gene[test,], lipid = nutrimouse$lipid[test,])

design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutrimouse.sgccda <- wrapper.sgccda(X = data,
                                    Y = Y,
                                    design = design,
                                    keepX = list(gene = c(10,10), lipid = c(15,15)),
                                    ncomp = 3,
                                    scheme = "centroid", tol=1e-30, init="svd")

auroc(nutrimouse.sgccda,plot = TRUE,roc.block = 1,roc.comp = 1)
auroc(nutrimouse.sgccda,plot = TRUE,roc.block = 2,roc.comp = 1)
auroc(nutrimouse.sgccda,plot = TRUE,roc.block = 1,roc.comp = 2)
auroc(nutrimouse.sgccda,plot = TRUE,roc.block = 2,roc.comp = 2)
auroc(nutrimouse.sgccda,plot = TRUE,roc.block = 1,roc.comp = 3)
auroc(nutrimouse.sgccda,plot = TRUE,roc.block = 2,roc.comp = 3)
auroc(nutrimouse.sgccda,plot = FALSE,roc.block = 1,roc.comp = 1)
auroc(nutrimouse.sgccda,plot = FALSE,roc.block = 2,roc.comp = 1)
auroc(nutrimouse.sgccda,plot = FALSE,roc.block = 1,roc.comp = 2)
auroc(nutrimouse.sgccda,plot = FALSE,roc.block = 2,roc.comp = 2)
auroc(nutrimouse.sgccda,plot = FALSE,roc.block = 1,roc.comp = 3)
auroc(nutrimouse.sgccda,plot = FALSE,roc.block = 2,roc.comp = 3)

###PROBLEME FACTOR (levels de facteurs non prÃ©sent): RESOLU en ajoutant factor


auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.block = 1,roc.comp = 1)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.block = 2,roc.comp = 1)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.block = 1,roc.comp = 2)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.block = 2,roc.comp = 2)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.block = 1,roc.comp = 3)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = TRUE,roc.block = 2,roc.comp = 3)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.block = 1,roc.comp = 1)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.block = 2,roc.comp = 1)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.block = 1,roc.comp = 2)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.block = 2,roc.comp = 2)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.block = 1,roc.comp = 3)
auroc(nutrimouse.sgccda,newdata = data.test,outcome.test = as.factor(Y.test),plot = FALSE,roc.block = 2,roc.comp = 3)


############################ e<-perf



# plsda/splsda
# ----

#library(mixOmics)
data(breast.tumors)
test=sample(1:47,5,replace=FALSE)
X <- breast.tumors$gene.exp
X.test<-breast.tumors$gene.exp[test,]
Y <- breast.tumors$sample$treatment
Y.test<-breast.tumors$sample$treatment[test]

res.plsda <- plsda(X, Y, ncomp = 2)

####ERROR A RESOUDRE MFOLD NREPEAT 1

for(validation in c("Mfold","loo"))
{
  for(nrepeat in 1:2)
  {
      #print(paste("validation ",validation, "nrepeat", nrepeat, "\n"))
    e<-perf(res.plsda,dist = c("max.dist"),
            auc = TRUE,validation = validation,nrepeat=nrepeat, progressBar = FALSE)
    e<-perf(res.plsda,dist = c("max.dist"),
            auc = FALSE,validation = validation,nrepeat=nrepeat, progressBar = FALSE)
    
    e<-perf(res.plsda,dist = c("all"),
            auc = TRUE,validation = validation,nrepeat=nrepeat, progressBar = FALSE)
    e<-perf(res.plsda,dist = c("all"),
            auc = FALSE,validation = validation,nrepeat=nrepeat, progressBar = FALSE)}
}


res.plsda <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))


res.plsda <- plsda(X, Y, ncomp = 2)

for(validation in c("Mfold","loo"))
{
  for(nrepeat in 1:2)
  {
      #print(paste("validation ",validation, "nrepeat", nrepeat, "\n"))
    e<-perf(res.plsda,dist = c("max.dist"),
            auc = TRUE,validation = validation,nrepeat=nrepeat, progressBar = FALSE)
    
    e<-perf(res.plsda,dist = c("all"),
            auc = TRUE,validation = validation,nrepeat=nrepeat, progressBar = FALSE)
}
}


data(liver.toxicity)

test=sample(1:64,7,replace=FALSE)
X <- liver.toxicity$gene
X.test<-liver.toxicity$gene[test,]
Y <- liver.toxicity$treatment[, 4]
Y.test <- liver.toxicity$treatment[test, 4]

res.plsda <- plsda(X, Y, ncomp = 2)


res.plsda <- plsda(X, Y, ncomp = 3)


# mint.plsda
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

res.plsda=mint.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)


auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 1, roc.study = 3)


e<-perf(res.plsda,dist = c("max.dist"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("max.dist"),
         auc = TRUE, progressBar = FALSE)

e<-perf(res.plsda,dist = c("all"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("all"),
         auc = TRUE, progressBar = FALSE)

e<-perf(res.plsda,dist = c("centroids.dist"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("centroids.dist"),
         auc = TRUE, progressBar = FALSE)

e<-perf(res.plsda,dist = c("mahalanobis.dist"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("mahalanobis.dist"),
         auc = TRUE, progressBar = FALSE)





res.plsda=mint.plsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,study=exp)

e<-perf(res.plsda,dist = c("max.dist"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("max.dist"),
         auc = TRUE, progressBar = FALSE)

e<-perf(res.plsda,dist = c("all"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("all"),
         auc = TRUE, progressBar = FALSE)

e<-perf(res.plsda,dist = c("centroids.dist"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("centroids.dist"),
         auc = TRUE, progressBar = FALSE)

e<-perf(res.plsda,dist = c("mahalanobis.dist"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("mahalanobis.dist"),
         auc = TRUE, progressBar = FALSE)


res.plsda=mint.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)

e<-perf(res.plsda,dist = c("max.dist"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("max.dist"),
         auc = TRUE, progressBar = FALSE)

e<-perf(res.plsda,dist = c("all"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("all"),
         auc = TRUE, progressBar = FALSE)

e<-perf(res.plsda,dist = c("centroids.dist"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("centroids.dist"),
         auc = TRUE, progressBar = FALSE)

e<-perf(res.plsda,dist = c("mahalanobis.dist"),
         auc = FALSE, progressBar = FALSE)
e<-perf(res.plsda,dist = c("mahalanobis.dist"),
         auc = TRUE, progressBar = FALSE)




############################ simulating data to check AUC <0.5

X = matrix(rnorm(50*100), nrow=50)
Y = factor(c(rep("a",25),rep("b",25)))

mod = plsda(X,Y,ncomp=2)
cv <- perf(mod, validation = "loo", auc = TRUE, progressBar = FALSE)

par(opar)
par(mfrow=c(1,1))
par(new=FALSE)
plot(1:10)
plot(1:10)

cv$auc #should be <0.5
if(cv$auc[[2]][1]>=0.5)
stop("random AUC should be less than 0.5. This is a extra check -see last line of R code")

