#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
if(FALSE)
{
library(rARPACK)
library(matrixStats)
library(mixOmicsDD)

#source("mixOmics/R/plotIndiv.R")
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, outcome=Y, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",bias = FALSE,weighted=FALSE, nrepeat=3)

set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "factorial",verbose = FALSE,bias = FALSE,weighted=FALSE, nrepeat=3)




{set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",dist="max.dist",bias = FALSE,weighted=FALSE, nrepeat=3)
print(tune2$error.rate)

set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",dist="max.dist",bias = FALSE,weighted=TRUE, nrepeat=3)
print("weighted")
print(tune2$error.rate)


set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",dist="max.dist",bias = TRUE,weighted=TRUE, nrepeat=3)
print("bias")
print(tune2$error.rate)


set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "factorial",dist="max.dist",bias = TRUE,weighted=TRUE, nrepeat=3)
print("scheme")
print(tune2$error.rate)



set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "factorial",dist="max.dist",bias = TRUE,weighted=TRUE, nrepeat=3, init = "svd")
print("svd")
print(tune2$error.rate)

set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "factorial",dist="max.dist",bias = TRUE,weighted=TRUE, nrepeat=3, init = "svd.single")
print("svd.single")
print(tune2$error.rate)

set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",dist="max.dist",bias = TRUE,weighted=TRUE, nrepeat=3, init = "svd.single")
print("centroid")
print(tune2$error.rate)


set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",dist="mahalanobis.dist",bias = TRUE,weighted=TRUE, nrepeat=3)
print("mahalanobis.dist")
print(tune2$error.rate)


set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",dist="centroids.dist",bias = TRUE,weighted=TRUE, nrepeat=3)
print("centroids.dist")
print(tune2$error.rate)

design2 = matrix(c(0,0,1,0,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design2,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",dist="centroids.dist",bias = TRUE,weighted=TRUE, nrepeat=3)
print("design2")
print(tune2$error.rate)




}




library(mixOmicsDD)
load("/Users/f.rohart/Downloads/Manuscript/mixomics_org_diablo/casestudy_brca/data/trainTestDatasetsNormalized.RDATA")
Y.train <- droplevels(pam50Train0$Call)
#names(Y.train) <- rownames(pam50Train0)
X.train <- list(mRNA = mrnaTrain0, miRNA = mirnaTrain0, CpGs = methTrain0, Proteins = protTrain0)
all(names(Y.train) == rownames(X.train[[1]]))
all(names(Y.train) == rownames(X.train[[2]]))
all(names(Y.train) == rownames(X.train[[3]]))
all(names(Y.train) == rownames(X.train[[4]]))
dim(X.train[[1]]); dim(X.train[[2]]); dim(X.train[[3]]); dim(X.train[[4]]);

design <- matrix(1/6, nrow = length(X.train), ncol = length(X.train))
rownames(design) <- colnames(design) <- names(X.train)
diag(design) <- 0
design= cbind(design,"Y"=1)
design=rbind(design,"Y"=1)
diag(design) <- 0


test.keepX = list(mRNA = c(10, 15, 20), miRNA = c(10, 15, 20), CpGs = c(10, 15, 20), Proteins = c(10, 15, 20))


set.seed(12)
system.time(tune2 <- tune.block.splsda(X = X.train,Y = Y.train,design = design,ncomp = 2,test.keepX = test.keepX,scheme = "centroid",dist="centroids.dist",bias = FALSE,weighted=TRUE, nrepeat=5, folds=5))



mod = block.splsda(X=data,
Y = Y,
design = design,
keepX = tune2$choice.keepX,
ncomp = 2,#c(2, 2),
scheme = "centroid",
bias = FALSE,
tol=1e-30)

pe=perf(mod)


nutrimouse.sgccda <- block.splsda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,#c(2, 2),
scheme = "centroid",
bias = FALSE,
tol=1e-30)


pred=predict(nutrimouse.sgccda, newdata=data)

pred=predict(nutrimouse.sgccda, newdata.scale=nutrimouse.sgccda$X)


load("temp.Rdata")
x=out.DA$class[[1]]
class.per.comp = lapply(1:min(ncomp), function(y) {matrix(sapply(x, function(z)  z[,y, drop = FALSE]),ncol=J)}) # combine the results per component
names(class.per.comp) = paste0("comp",1:min(ncomp))
class.per.comp = lapply(class.per.comp, function(y){rownames(y) = rownames(out.DA$MajorityVote[[1]]); colnames(y) = names(x); y})

y=class.per.comp[[2]] #second comp

z=y[2,]
temp = aggregate(weights,list(z),sum)



source("mixOmicsDD/R/internal_wrapper.mint.R")
source("mixOmicsDD/R/internal_wrapper.mint.block.R")
source("mixOmicsDD/R/internal_mint.block.R")
source("mixOmicsDD/R/internal_mint.block_helpers.R")
source("mixOmicsDD/R/check_entry.R")

#source("mixOmicsDD/R/tune.splsda.R")
source("mixOmicsDD/R/tune.block.splsda.R")
source("mixOmicsDD/R/block.splsda.R")
source("mixOmicsDD/R/block.spls.R")
source("mixOmicsDD/R/splsda.R")
source("mixOmicsDD/R/spls.R")
source("mixOmicsDD/R/pls.R")
source("mixOmicsDD/R/plsda.R")
source("mixOmicsDD/R/MCV.block.splsda.R")
source("mixOmicsDD/R/MCV.splsda.R")
source("mixOmicsDD/R/perf.R")

nutrimouse.sgccda <- block.splsda(X=data,
Y = Y,
design = design,
#keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,#c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE,
tol=1e-30)


source("mixOmicsDD/R/tune.block.splsda.R")
source("mixOmicsDD/R/MCV.block.splsda.R")
source("mixOmicsDD/R/internal_wrapper.mint.block.R")
source("mixOmicsDD/R/internal_mint.block.R")
source("mixOmicsDD/R/predict.mint.block.pls.R")
source("mixOmicsDD/R/internal_predict.DA.R")
source("mixOmicsDD/R/t.test.process.R")


set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",bias = FALSE,weighted=FALSE, measure="overall", init="svd")


tune2
}
