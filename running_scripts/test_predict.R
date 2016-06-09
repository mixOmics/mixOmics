## =========================================================================================================
## =========================================================================================================
## ===========================                  test predict               ==============================##
## =========================================================================================================
## =========================================================================================================


# example with pls
# ------------------
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, ncomp = 2, mode = "classic")

indiv1 <- c(200, 40, 60)
indiv2 <- c(190, 45, 45)
newdata <- rbind(indiv1, indiv2)
colnames(newdata) <- colnames(X)
newdata

pred <- predict(linn.pls, newdata)

plotIndiv(linn.pls, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction PLS")
points(pred$variates[, 1], pred$variates[, 2], pch = 19, cex = 1.2)
text(pred$variates[, 1], pred$variates[, 2],
c("new ind.1", "new ind.2"), pos = 3)

# example with plsda
# ------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- as.factor(liver.toxicity$treatment[, 4])

## if training is perfomed on 4/5th of the original data
samp <- sample(1:5, nrow(X), replace = TRUE)
test <- which(samp == 1)   # testing on the first fold
train <- setdiff(1:nrow(X), test)

plsda.train <- plsda(X[train, ], Y[train], ncomp = 2)
test.predict <- predict(plsda.train, X[test, ], method = "max.dist")

Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
cbind(Y = as.character(Y[test]), Prediction)

plotIndiv(plsda.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction PLS-DA")

points(test.predict$variates[, 1], test.predict$variates[, 2], pch = 19, cex = 1.2,col=color.mixo(as.numeric(factor(Prediction,levels=levels(Y)))))
text(test.predict$variates[, 1], test.predict$variates[, 2],
c(paste0("new ind",1:length(test))), pos = 3)


# example with splsda
# ------------------
splsda.train <- splsda(X[train, ], Y[train], ncomp = 2, keepX = c(30, 30))

test.predict <- predict(splsda.train, X[test, ], method = "max.dist")
Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
cbind(Y = as.character(Y[test]), Prediction)


plotIndiv(splsda.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction sPLS-DA")

points(test.predict$variates[, 1], test.predict$variates[, 2], pch = 19, cex = 1.2,col=color.mixo(as.numeric(factor(Prediction,levels=levels(Y)))))
text(test.predict$variates[, 1], test.predict$variates[, 2],
c(paste0("new ind",1:length(test))), pos = 3)


# ----------------------------
# load data for block analyses
# ----------------------------
data(nutrimouse)
# need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
Y.mat = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y.mat)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
byrow = TRUE, dimnames = list(names(data), names(data)))


# train on 75% of the data
ind.train=NULL
for(i in 1:nlevels(nutrimouse$diet))
ind.train=c(ind.train,which(nutrimouse$diet==levels(nutrimouse$diet)[i])[1:6])


#training set
gene.train=nutrimouse$gene[ind.train,]
lipid.train=nutrimouse$lipid[ind.train,]
Y.mat.train=Y.mat[ind.train,]
Y.train=nutrimouse$diet[ind.train]
data.train=list(gene=gene.train,lipid=lipid.train,Y=Y.mat.train)

#test set
gene.test=nutrimouse$gene[-ind.train,]
lipid.test=nutrimouse$lipid[-ind.train,]
Y.mat.test=Y.mat[-ind.train,]
Y.test=nutrimouse$diet[-ind.train]
data.test=list(gene=gene.test,lipid=lipid.test)



# example with block.pls
# ------------------
res.train=block.pls(X=data.train,indY=3,ncomp = 3)
test.predict <- predict(res.train, newdata=data.test, method = "max.dist")

par(mfrow=c(1,2))
plotIndiv(res.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction block.PLS",blocks="gene")
points(test.predict$variates[["gene"]][, 1], test.predict$variates[["gene"]][, 2], pch = 19, cex = 1.2)
text(test.predict$variates[["gene"]][, 1], test.predict$variates[["gene"]][, 2],
c(paste0("new ind",1:length(test))), pos = 3)

plotIndiv(res.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,blocks="lipid")
points(test.predict$variates[["lipid"]][, 1], test.predict$variates[["lipid"]][, 2], pch = 19, cex = 1.2)
text(test.predict$variates[["lipid"]][, 1], test.predict$variates[["lipid"]][, 2],
c(paste0("new ind",1:length(test))), pos = 3)


# example with block.spls
# ------------------
res.train=block.spls(X=data.train,indY=3,ncomp = 3,keepX=list(gene=c(10,10,10),lipid=c(5,5,5)))
test.predict <- predict(res.train, newdata=data.test, method = "max.dist")

par(mfrow=c(1,2))
plotIndiv(res.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction block.sPLS",blocks="gene")
points(test.predict$variates[["gene"]][, 1], test.predict$variates[["gene"]][, 2], pch = 19, cex = 1.2)
text(test.predict$variates[["gene"]][, 1], test.predict$variates[["gene"]][, 2],
c(paste0("new ind",1:length(test))), pos = 3)

plotIndiv(res.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,blocks="lipid")
points(test.predict$variates[["lipid"]][, 1], test.predict$variates[["lipid"]][, 2], pch = 19, cex = 1.2)
text(test.predict$variates[["lipid"]][, 1], test.predict$variates[["lipid"]][, 2],
c(paste0("new ind",1:length(test))), pos = 3)


# example with block.plsda
# ------------------
res.train=block.plsda(X=list(gene=gene.train,lipid=lipid.train),Y=Y.train,ncomp = 3)
test.predict <- predict(res.train, newdata=data.test, method = "max.dist")

Prediction <- test.predict$vote$max.dist[, 3]
color.test=color.mixo(as.numeric(factor(Prediction,levels=levels(Y.train))))
color.test[is.na(color.test)]=1

par(mfrow=c(1,2))
plotIndiv(res.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction block.sPLS-DA",blocks="gene")
points(test.predict$variates[["gene"]][, 1], test.predict$variates[["gene"]][, 2], pch = 19, cex = 1.2,col=color.test)
text(test.predict$variates[["gene"]][, 1], test.predict$variates[["gene"]][, 2],
c(paste0("new ind",1:length(test))), pos = 3)


plotIndiv(res.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,blocks="lipid")
points(test.predict$variates[["lipid"]][, 1], test.predict$variates[["lipid"]][, 2], pch = 19, cex = 1.2,col=color.test)
text(test.predict$variates[["lipid"]][, 1], test.predict$variates[["lipid"]][, 2],
c(paste0("new ind",1:length(test))), pos = 3)


# need to code the majority vote for block.plsda and block.splsda, to add in predict.mint.block.pls

#Prediction <- levels(Y.train)[test.predict$class$max.dist[, 2]]
#cbind(Y = as.character(Y[test]), Prediction)

# example with block.splsda=diablo=sgccda
# ------------------
res.train=block.splsda(X=list(gene=gene.train,lipid=lipid.train),Y=Y.train,ncomp = 3,keepX=list(gene=c(10,10,10),lipid=c(5,5,5)))
test.predict <- predict(res.train, newdata=data.test, method = "max.dist")

Prediction <- test.predict$vote$max.dist[, 3]
color.test=color.mixo(as.numeric(factor(Prediction,levels=levels(Y.train))))
color.test[is.na(color.test)]=1

par(mfrow=c(1,2))
plotIndiv(res.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction block.sPLS-DA",blocks="gene")
points(test.predict$variates[["gene"]][, 1], test.predict$variates[["gene"]][, 2], pch = 19, cex = 1.2,col=color.test)
text(test.predict$variates[["gene"]][, 1], test.predict$variates[["gene"]][, 2],
c(paste0("new ind",1:length(test))), pos = 3)


plotIndiv(res.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,blocks="lipid")
points(test.predict$variates[["lipid"]][, 1], test.predict$variates[["lipid"]][, 2], pch = 19, cex = 1.2,col=color.test)
text(test.predict$variates[["lipid"]][, 1], test.predict$variates[["lipid"]][, 2],
c(paste0("new ind",1:length(test))), pos = 3)


# example with block.splsda=diablo=sgccda and a missing block
res.train=block.splsda(X=list(gene=gene.train,lipid=lipid.train),Y=Y.train,ncomp = 3,keepX=list(gene=c(10,10,10),lipid=c(5,5,5)))
test.predict <- predict(res.train, newdata=data.test[2], method = "max.dist")

Prediction <- test.predict$vote$max.dist[, 3]
color.test=color.mixo(as.numeric(factor(Prediction,levels=levels(Y.train))))
color.test[is.na(color.test)]=1



# ----------------------------
# load data for MINT analyses
# ----------------------------
data(stemcells) #load gene, celltype and study

gene=stemcells$gene
celltype=stemcells$celltype
study=stemcells$study

# unmap for mint.pls approach
Y.mat=unmap(celltype)
rownames(Y.mat)=rownames(gene)


## if training is perfomed on 4 out of 6studies
ind.train=which(study%in%c(2:6))
gene.train=gene[ind.train,]
study.train=factor(study[ind.train])
celltype.train=celltype[ind.train]
Y.mat.train=Y.mat[ind.train,]

ind.test=which(study%in%c(1))
gene.test=gene[ind.test,]
study.test=factor(study[ind.test])
celltype.test=celltype[ind.test]
Y.mat.test=Y.mat[ind.test,]

par(mfrow=c(1,1))
# example with mint.pls
# ------------------
mint.train=mint.pls(X=gene.train,Y=Y.mat.train,ncomp=3,near.zero.var=FALSE,study=study.train)

test.predict <- predict(mint.train, gene.test, method = "max.dist",study.test=study.test)

plotIndiv(mint.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction mint.PLS",xlim=c(-20,50))
points(test.predict$variates[, 1], test.predict$variates[, 2], pch = 19, cex = 1.2)
text(test.predict$variates[, 1], test.predict$variates[, 2],
c(paste0("new ind",1:length(study.test))), pos = 3)


# example with mint.spls
# ------------------
mint.train=mint.spls(X=gene.train,Y=Y.mat.train,ncomp=3,near.zero.var=FALSE,study=study.train,keepX=c(10,10,10))

test.predict <- predict(mint.train, gene.test, method = "max.dist",study.test=study.test)

plotIndiv(mint.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction mint.sPLS",xlim=c(-3,6))
points(test.predict$variates[, 1], test.predict$variates[, 2], pch = 19, cex = 1.2)
text(test.predict$variates[, 1], test.predict$variates[, 2],
c(paste0("new ind",1:length(study.test))), pos = 3)



# example with mint.plsda
# ------------------
mint.train=mint.plsda(X=gene.train,Y=celltype.train,ncomp=3,near.zero.var=FALSE,study=study.train)

test.predict <- predict(mint.train, gene.test, method = "max.dist",study.test=study.test)
Prediction <- test.predict$class$max.dist[, 2]
cbind(Y = as.character(celltype.test), Prediction)

plotIndiv(mint.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction mint.PLS-DA")
points(test.predict$variates[, 1], test.predict$variates[, 2], pch = 19, cex = 1.2,col=color.mixo(as.numeric(factor(Prediction,levels=levels(celltype)))))
#text(test.predict$variates[, 1], test.predict$variates[, 2],
#c(paste0("new ind",1:length(study.test))), pos = 3)



# example with mint.splsda
# ------------------
mint.train=mint.splsda(X=gene.train,Y=celltype.train,ncomp=3,near.zero.var=FALSE,study=study.train,keepX=c(10,10,10))

test.predict <- predict(mint.train, gene.test, method = "max.dist",study.test=study.test)
Prediction <- test.predict$class$max.dist[, 2]
cbind(Y = as.character(celltype.test), Prediction)

plotIndiv(mint.train, comp = 1:2, rep.space = "X-variate",style="graphics",ind.names=FALSE,title="Prediction mint.sPLS-DA")
points(test.predict$variates[, 1], test.predict$variates[, 2], pch = 19, cex = 1.2,col=color.mixo(as.numeric(factor(Prediction,levels=levels(celltype)))))
#text(test.predict$variates[, 1], test.predict$variates[, 2],
#c(paste0("new ind",1:length(study.test))), pos = 3)

# example with mint.block.pls
# ------------------

# example with mint.block.spls
# ------------------

# example with mint.block.plsda
# ------------------

# example with mint.block.splsda
# ------------------





