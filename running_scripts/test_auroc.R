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

library(mixOmics)
data(breast.tumors)
test=sample(1:47,5,replace=TRUE)
X <- breast.tumors$gene.exp
X.test<-breast.tumors$gene.exp[test,]
Y <- breast.tumors$sample$treatment
Y.test<-breast.tumors$sample$treatment[test]

res.plsda <- plsda(X, Y, ncomp = 2)

for(roc.comp in 1:2){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,plot = plot,roc.comp = roc.comp)
    print(a)
}
}

for(roc.comp in 1:2){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = plot,roc.comp = roc.comp)
    print(a)
  }
}


res.plsda <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))

for(roc.comp in 1:2){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,plot = plot,roc.comp = roc.comp)
    print(a)
  }
}

for(roc.comp in 1:2){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = plot,roc.comp = roc.comp)
    print(a)
  }
}




data(liver.toxicity)

test=sample(1:64,7,replace=TRUE)
X <- liver.toxicity$gene
X.test<-liver.toxicity$gene[test,]
Y <- liver.toxicity$treatment[, 4]
Y.test <- liver.toxicity$treatment[test, 4]

res.plsda <- plsda(X, Y, ncomp = 2)

for(roc.comp in 1:2){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,plot = plot,roc.comp = roc.comp)
    print(a)
  }
}

for(roc.comp in 1:2){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = plot,roc.comp = roc.comp)
    print(a)
  }
}


res.plsda <- plsda(X, Y, ncomp = 3)

for(roc.comp in 1:3){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,plot = plot,roc.comp = roc.comp)
    print(a)
  }
}

for(roc.comp in 1:3){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,newdata = X.test,outcome.test = as.factor(Y.test),plot = plot,roc.comp = roc.comp)
    print(a)
  }
}


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

for(roc.comp in 1:3){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,plot = plot,roc.comp = roc.comp)
    print(a)
  }
}

for(roc.comp in 1:3){
  for(plot in c(TRUE,FALSE)){
    a=auroc(res.plsda,newdata = data.light,outcome.test = as.factor(type.id.light),study.test = study.light,plot = plot,roc.comp = roc.comp)
print(a)
  }
}


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


res.plsda<-mint.plsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,study=exp)

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



res.plsda<-mint.block.splsda(list(XX=data),Y=type.id,ncomp=2)

####PROBLEME ERROR//RESOLU

auroc(res.plsda,plot = TRUE,roc.comp = 1)
auroc(res.plsda,plot = TRUE,roc.comp = 2)
auroc(res.plsda,plot = FALSE,roc.comp = 1)
auroc(res.plsda,plot = FALSE,roc.comp = 2)

auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = TRUE,roc.comp = 2)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 1)
auroc(res.plsda,newdata = list(XX=data.light),outcome.test = as.factor(type.id.light),study.test = study.light,plot = FALSE,roc.comp = 2)



res.plsda<-mint.block.splsda(X=list(XX=data,Y=type.id),indY=2,keepX=list(XX=c(10,5),ncomp = 2))

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

data(nutrimouse)
train=sample(1:40,40,replace=TRUE)
test=sample(1:40,4,replace=TRUE)
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
                                    scheme = "centroid",
                                    verbose = FALSE,
                                    bias = FALSE)

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

###PROBLEME FACTOR (levels de facteurs non présent): RESOLU en ajoutant factor


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


############################perf



# plsda/splsda
# ----

library(mixOmics)
data(breast.tumors)
test=sample(1:47,5,replace=TRUE)
X <- breast.tumors$gene.exp
X.test<-breast.tumors$gene.exp[test,]
Y <- breast.tumors$sample$treatment
Y.test<-breast.tumors$sample$treatment[test]

res.plsda <- plsda(X, Y, ncomp = 2)

####ERROR A RESOUDRE MFOLD NREPEAT 1 // RESOLU PAR FR

####ERROR CONSTRAINT //RESOLU


for(validation in c("Mfold","loo"))
{
  for(nrepeat in 1:2)
{
    for(auc in c(TRUE,FALSE))
    { 
      for(dist in c("all", "max.dist", "centroids.dist", "mahalanobis.dist"))
      {
        for(constraint in c(TRUE,FALSE))
        {
          e<-perf(res.plsda,dist = dist,constraint=constraint,
                  auc = auc,validation = validation,nrepeat=nrepeat)
        }
      }
      
    }
}
}

res.plsda <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))

####ERROR 

for(validation in c("Mfold","loo"))
{
  for(nrepeat in 1:2)
  {
    for(auc in c(TRUE,FALSE))
    { 
      for(dist in c("all", "max.dist", "centroids.dist", "mahalanobis.dist"))
      {
        for(constraint in c(TRUE,FALSE))
        {
          e<-perf(res.plsda,dist = dist,constraint=constraint,
                  auc = auc,validation = validation,nrepeat=nrepeat)
        }
      }
      
    }
  }
}

res.plsda <- plsda(X, Y, ncomp = 2)

####ERROR 

for(validation in c("Mfold","loo"))
{
  for(nrepeat in 1:2)
  {
    for(auc in c(TRUE,FALSE))
    { 
      for(dist in c("all", "max.dist", "centroids.dist", "mahalanobis.dist"))
      {
         for(constraint in c(TRUE,FALSE))
         {
          e<-perf(res.plsda,dist = dist,constraint=constraint,
                  auc = auc,validation = validation,nrepeat=nrepeat)
         }
      }
      
    }
  }
}
                    
                    

data(liver.toxicity)

test=sample(1:64,7,replace=TRUE)
X <- liver.toxicity$gene
X.test<-liver.toxicity$gene[test,]
Y <- liver.toxicity$treatment[, 4]
Y.test <- liver.toxicity$treatment[test, 4]

res.plsda <- plsda(X, Y, ncomp = 2)

for(validation in c("Mfold","loo"))
{
  for(nrepeat in 1:2)
  {
    for(auc in c(TRUE,FALSE))
    { 
      for(dist in c("all", "max.dist", "centroids.dist", "mahalanobis.dist"))
      {
        for(constraint in c(TRUE,FALSE))
        {
          e<-perf(res.plsda,dist = dist,constraint=constraint,
                  auc = auc,validation = validation,nrepeat=nrepeat)
        }
      }
      
    }
  }
}

res.plsda <- plsda(X, Y, ncomp = 3)

for(validation in c("Mfold","loo"))
{
  for(nrepeat in 1:2)
  {
    for(auc in c(TRUE,FALSE))
    { 
      for(dist in c("all", "max.dist", "centroids.dist", "mahalanobis.dist"))
      {
        for(constraint in c(TRUE,FALSE))
        {
          e<-perf(res.plsda,dist = dist,constraint=constraint,
                  auc = auc,validation = validation,nrepeat=nrepeat)
        }
      }
      
    }
  }
}

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

for(auc in c(TRUE,FALSE))
{ 
  for(dist in c("all", "max.dist", "centroids.dist", "mahalanobis.dist"))
  {
    for(constraint in c(TRUE,FALSE))
    {
      e<-perf(res.plsda,dist = dist,constraint=constraint,
              auc = auc)
      print(e)
    }
  }
  
}



res.plsda=mint.plsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,study=exp)

for(auc in c(TRUE,FALSE))
{ 
  for(dist in c("all", "max.dist", "centroids.dist", "mahalanobis.dist"))
  {
    for(constraint in c(TRUE,FALSE))
    {
      e<-perf(res.plsda,dist = dist,constraint=constraint,
              auc = auc)
      print(e)
    }
  }
  
}


res.plsda=mint.splsda(X=data,Y=type.id,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=exp)


    for(auc in c(TRUE,FALSE))
   { 
      for(dist in c("all", "max.dist", "centroids.dist", "mahalanobis.dist"))
      {
        for(constraint in c(TRUE,FALSE))
        {
          e<-perf(res.plsda,dist = dist,constraint=constraint,
            auc = auc)
    print(e)
        }
      }
    
    }