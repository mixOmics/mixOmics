#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################

library(mixOmics)
library(mixOmicsD)
library(mixOmicsDD)

set.seed(12)
n=1000
p=10000

system.time(X <- matrix(rnorm(n*p),nrow=n, ncol=p))
beta = rep(0,p)
beta[sample(1:p,10)]=100

#X[10]=NA
NAA=TRUE
if(NAA)
{
    X[sample(1:length(X), 0.01*n*p)] = NA

    X2 = X
    X2[which(is.na(X))] = 0
    temp= X2%*% beta +rnorm(n)
} else {
    temp= X%*% beta +rnorm(n)

}
Y = rep(0,n)
Y[temp>quantile(Y,0.7)]=1

a=crossprod(X,unmap(Y))
b=svd(a, nu=1,nv=1)


system.time(out <- splsda(X,Y, ncomp=2, keepX=c(10,5)))
system.time(out <- splsda(X,Y, ncomp=1, keepX=c(10)))

selectVar(out)


source("mixOmicsD/R/internal_mint.block_helpers.R")
source("mixOmicsD/R/check_entry.R")

system.time(X_scale <- scale(X,scale=TRUE))
system.time(X_scale1 <- scale.function(X,scale=TRUE))

system.time(X_scale2 <- scale.function2(X,scale=TRUE))



source("mixOmicsDD/R/internal_mint.block.R")
source("mixOmicsDD/R/internal_mint.block_helpers.R")
test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3))



library(mixOmicsD)
load("../Diablo/childhoodAsthmaDatasets.RDATA")
X.start <- list(mRNA = mrnaAsthma, CpGs = methAsthma)

lapply(X.start,dim)

#filtering CpGs
#for(i in 1:length(X.start))
#{
#    temp = X.start[[i]]
#    mad_temp = apply(temp, 2, mad)
#    temp <- temp[, names(apply(temp, 2, mad)[order(mad_temp, decreasing = TRUE)][1:30000])]
#    X.start[[i]] = temp
#}
#
lapply(X.start,dim)

colnames(X.start[[1]]) = 1:ncol(X.start[[1]])
colnames(X.start[[2]]) = 1:ncol(X.start[[2]])

Y.start <- factor(groupAsthma, levels = c("Nonasthmatic", "Asthmatic"))



ncomp = 2

#design with 0.1 for blocks
design=matrix(0.1, ncol=3,nrow = 3)
design[3,] = design[,3]=1
diag(design)=0
design

# set parameters
keepXList.diablo = list()
for(indice in 1:length(X.start))
{
    keepXList.diablo[[indice]] = rep(100,ncomp)
}
names(keepXList.diablo) = names(X.start)
progressBar = TRUE
system.time(mod.diablo <- block.splsda(X = X.start, Y = Y.start, ncomp = ncomp, keepX = keepXList.diablo, scheme = "horst", design = design))



library(rARPACK)
library(matrixStats)
A=X.start
A[[3]] = unmap(Y.start)
J=3

source("mixOmicsD/R/internal_mint.block_helpers.R")
system.time(A_scale <- lapply(A, scale.function))

system.time(alpha <-  lapply(1 : J, function(y){initsvd(A_scale[[y]]$data)}))

system.time(alpha <-  lapply(1 : J, function(y){initsvd(A[[y]])}))



ind = 2

system.time(a<-initsvd(A[[ind]]))
system.time(a<-svd(A[[ind]]))
system.time(a<-svd(t(A[[ind]])))

system.time(b<-svds(A[[3]], k=1, nu=0, nv=1)$u)
system.time(b<-svds(t(A[[ind]])))



# WITH NA
library(mixOmicsD)

data(breast.tumors)
X <- breast.tumors$gene.exp
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(breast.tumors$sample$treatment)

sum(is.na(X))

system.time(res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25)))
