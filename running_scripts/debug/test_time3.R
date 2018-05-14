#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################

library(mixOmics)
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}
sourceDir("mixOmics/R/",trace=FALSE)
library(rARPACK)
library(matrixStats)

library(profvis)
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, outcome=Y, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)



profvis({
    nutrimouse.sgccda <- block.splsda(X=data,
    Y = Y,
    design = design,
    keepX = list(gene=c(10,10), lipid=c(15,15)),
    ncomp = 2,#c(2, 2),
    scheme = "centroid",
    tol=1e-30)
})

#simu
library(mixOmics)
n=1000;p1=10000;p2=20000
X1 = matrix(rnorm(n*p1),nrow=n)
X2 = matrix(rnorm(n*p2),nrow=n)
rownames(X1) = rownames(X2) = paste0("n",1:n)
colnames(X1) = paste0("p1",1:p1)
colnames(X2) = paste0("p2",1:p2)
Y = rbinom(n,1,0.5)
XNA1=X1
XNA1[sample(1:(n*p1),100)]=NA
XNA2=X2
XNA2[sample(1:(n*p2),100)]=NA

data = list(X1=X1,X2=X2)
dataNA = list(XNA1=XNA1,XNA2=XNA2)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

sourceDir("mixOmics/R/",trace=FALSE)

profvis({
    mod <- block.splsda(X=data,Y = Y,design = design,ncomp = 1,scheme = "centroid",max.iter=2)
})

profvis({
    mod2 <- perf(mod,folds=3)
})


profvis({
    mod <- block.splsda(X=dataNA,
    Y = Y,
    design = design,
    #keepX = list(X1=c(10,10), X2=c(15,15)),
    ncomp = 1,#c(2, 2),
    scheme = "centroid",max.iter=2)
})

profvis({
    tune = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 1,scheme = "centroid",progressBar = TRUE,folds=3,max.iter=2)
})

profvis({
    tune = tune.block.splsda(X = dataNA,Y = Y,design = design,ncomp = 1,scheme = "centroid",progressBar = TRUE,folds=3,max.iter=2)
})



library(profmem)

library(microbenchmark)
library(mixOmics)
library(profvis)
n=100;p=1000
X = matrix(rnorm(n*p),nrow=n)
rownames(X) = paste(1:n)
colnames(X) = paste(1:p)
Y = rbinom(n,1,0.5)
mod=splsda(X,Y)


library(mixOmics)
library(profvis)
n=2000;p=100000
X = matrix(rnorm(n*p),nrow=n)
rownames(X) = paste(1:n)
colnames(X) = paste(1:p)
Y = rbinom(n,1,0.5)

#X[sample(1:1000,10)]=NA

#mod=splsda(X,Y)

profvis({
mod=splsda(X,Y)
})


microbenchmark(splsda(X,Y),times=10)

sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}
sourceDir("mixOmics/R/",trace=FALSE)
library(rARPACK)
library(matrixStats)

profvis({
    mod2=splsda(X,Y)
})
microbenchmark(mixOmics::splsda(X,Y),splsda(X,Y),times=2)



# with NA
library(microbenchmark)
library(mixOmics)
library(profvis)
n=100;p=1000
n=1000;p=50000
X = matrix(rnorm(n*p),nrow=n)
rownames(X) = paste(1:n)
colnames(X) = paste(1:p)
Y = rbinom(n,1,0.5)
XNA=X
XNA[sample(1:(n*p),100)]=NA

microbenchmark(mixOmics::splsda(XNA,Y),times=2)
profvis({
    mod=mixOmics::splsda(XNA,Y)
})

sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}
sourceDir("mixOmics/R/",trace=FALSE)
library(rARPACK)
library(matrixStats)
#mod=splsda(XNA,Y)

n=100;p=1000
X = matrix(rnorm(n*p),nrow=n)
rownames(X) = paste(1:n)
colnames(X) = paste(1:p)

Y = matrix(rnorm(n*2,1,0.5),ncol=2)
colnames(Y) = paste0("Y",1:ncol(Y))
tun = tune.spls(X,Y,ncomp=4, test.keepX = c(5,10,15), nrepeat=5, folds=3, test.keepY=2)



microbenchmark(splsda(XNA,Y),times=2)

profvis({
    mod=splsda(XNA,Y)
})

profvis({
    mod=splsda(X,Y)
})


profvis({
    tune1 <- tune.splsda(X,Y,ncomp=1,nrepeat=1,logratio="none",test.keepX = c(5,10,15),folds=3,dist="max.dist", progressBar = progressBar)
})


profvis({
    tune1 <- mixOmics::tune.splsda(X,Y,ncomp=1,nrepeat=1,logratio="none",test.keepX = c(5,10,15),folds=3,dist="max.dist", progressBar = progressBar)
})

profvis({
    tune1 <- mixOmics::tune.splsda(XNA,Y,ncomp=1,nrepeat=1,logratio="none",test.keepX = c(5,10,15),folds=3,dist="max.dist", progressBar = progressBar)
})

profvis({
    tune2 <- tune.splsda(XNA,Y,ncomp=1,nrepeat=1,logratio="none",test.keepX = c(5,10,15),folds=3,dist="max.dist", progressBar = progressBar)
})


progressBar = TRUE
library(mixOmics)
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- as.factor(breast.tumors$sample$treatment)

profvis({
    tune1 <- mixOmics::tune.splsda(X,Y,ncomp=3,nrepeat=1,logratio="none",test.keepX = seq(5,100,5),folds=3,dist="max.dist", progressBar = progressBar)
})






sourceDir("mixOmics/R/",trace=FALSE)
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
Y=Y[,9]

set.seed(3)
tun = mixOmics::tune.spls(X,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=5, folds=10)

sourceDir("mixOmics/R/",trace=FALSE)
set.seed(3)
tun2=tune.spls(X,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=5, folds=10)

all.equal(tun2,tun)


profvis({
tun = tune.spls(X,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=5, folds=3)
})



library(profvis)
library(mixOmics)
n=100;p=1000
X = matrix(rnorm(n*p),nrow=n)
rownames(X) = paste(1:n)
colnames(X) = paste(1:p)
Y = as.matrix(rnorm(n,sd=2))
XNA=X
XNA[sample(1:(n*p),100)]=NA

tune.spls(XNA,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=1, folds=3)

set.seed(3)
profvis({
    tun1 = mixOmics::tune.spls(X,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=1, folds=3)
})

set.seed(3)
profvis({
    tun1NA = mixOmics::tune.spls(XNA,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=1, folds=3)
})


sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}
sourceDir("mixOmics/R/",trace=FALSE)
library(rARPACK)
library(matrixStats)

set.seed(3)
profvis({
    tun2 = tune.spls(X,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=1, folds=3)
})


set.seed(3)
profvis({
    tun2NA = tune.spls(XNA,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=1, folds=3)
})

all.equal(tun1,tun2)
all.equal(tun1NA,tun2NA)




sourceDir("mixOmics/R/",trace=FALSE)
tun = tune.spls(XNA,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=1, folds=3)







microbenchmark(mixOmics::splsda(XNA,Y),splsda(XNA,Y),times=2)


microbenchmark(splsda(XNA,Y),splsda(X,Y))


save(list=ls(),file="temp.Rdata")

load("temp.Rdata")
Rtemp=R
profvis({
    library(microbenchmark)
    microbenchmark(
    Rtemp = lapply(1:J, function(q){if(misdata[q]) {replace(Rtemp[[q]], is.na.A[[q]], 0)}else{Rtemp[[q]]}}) # if missing data, R is the one replace by 0 where NA are supposed to be
,
if(misdata.all)# replace NA in A[[q]] by 0
for(j in c(1:J)[misdata])
Rtemp[[j]][is.na.A[[j]]]=0,

times=100)
})




source("mixOmics/R/check_entry.R")
source("mixOmics/R/splsda_allinone.R")
source("mixOmics/R/internal_mint.block_helpers.R")
library(rARPACK)
library(matrixStats)

profvis({
    mod=mean_centering_per_study(X,study=rep(1,n),scale=TRUE)
})

source("mixOmics/R/internal_mint.block_helpers.R")
profvis({
    data=X
    study=rep(1,n)
    scale=TRUE
    data.list.study = study_split(data, study)
    
    # center and scale data per group, and concatene the data
    #res = lapply(data.list.study, scale.function2, scale = scale)
    res = scale.function(data.list.study[[1]], scale=scale)
    res2 = scale.function2(data.list.study[[1]], scale=scale)
    
    print(all.equal(res,res2))
})

    temp=data.list.study[[1]]
    meanX = colMeans(temp, na.rm = TRUE)

    sqrt.sdX = colSds(temp, na.rm=TRUE)
    sqrt.sdX2 = colSds(temp,center=meanX, na.rm=TRUE)

})




colSdColMeans <- function(x, na.rm=TRUE) {
    if (na.rm) {
        n <- colSums(!is.na(x)) # thanks @flodel
    } else {
        n <- nrow(x)
    }
    colVar <- colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
    return(sqrt(colVar * n/(n-1)))
}


profvis({
mod2=splsda_allin(X,Y)
})

abc <- function(X){
    X.new=X
    ab(X=X.new,x=2)
    print(ls(),envir=environment())
}

ab <- function(X,x){
    XX=X+x
    print(ls(),envir=environment())
    print(ls(),envir=parent.env(environment()))
    rm(X.new,envir=parent.frame())
    print(ls())
    return(XX)
}
 XX=abc(X)

progressBar=TRUE
p<-profmem(tune1 <- tune.splsda(X,Y,ncomp=2,nrepeat=1,logratio="none",test.keepX = seq(5,100,5),folds=3,dist="max.dist", progressBar = progressBar,cpus=4))


Rprof(tf <- "rprof.log", memory.profiling=TRUE)
tune1 <- tune.splsda(X,Y,ncomp=2,nrepeat=1,logratio="none",test.keepX = seq(5,100,5),folds=3,dist="max.dist", progressBar = progressBar,cpus=4)
Rprof(NULL)
summaryRprof(tf)


Rprofmem("Rprofmem.out", threshold = 1000)
tune1 <- tune.splsda(X,Y,ncomp=2,nrepeat=1,logratio="none",test.keepX = seq(5,100,5),folds=3,dist="max.dist", progressBar = progressBar,cpus=2)
Rprofmem(NULL)
noquote(readLines("Rprofmem.out", n = 5))


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
