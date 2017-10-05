#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################

library(mixOmicsDD)

set.seed(12)
n=1000
p=10000

system.time(X <- matrix(rnorm(n*p),nrow=n, ncol=p))
beta = rep(0,p)
beta[sample(1:p,10)]=100

colnames(X) = paste0("X", 1:ncol(X))
temp= scale(X, scale=TRUE)%*% beta +rnorm(n)
Y = rep(0,n)
Y[temp>quantile(Y,0.7)]=1


#system.time(out <- splsda(X,Y, ncomp=2, keepX=c(10,10)))
library(rARPACK)
library(matrixStats)

dist = "max.dist"
tol=1e-6
max.iter=100
scale=TRUE
near.zero.var=TRUE
study=factor(rep(1,80))
ncomp=3

colnames(X) = paste0("X",1:ncol(X))
keepX = NULL#c(10,15)
test.keepX = 1:30
keepY=2
test.keepY=2

omit = 81:100
X.train = X[-omit, ]
Y.train = Y[-omit]
Y.train.mat = unmap(Y.train)
colnames(Y.train.mat) = levels(Y.train)
X.test = X[omit, , drop = FALSE]#matrix(X[omit, ], nrow = length(omit)) #removed to keep the colnames in X.test
Y.test = Y[omit]

source("mixOmicsDD/R/internal_wrapper.mint.R")
source("mixOmicsDD/R/internal_mint.block.R")
source("mixOmicsDD/R/internal_mint.block_helpers.R")
source("mixOmicsDD/R/check_entry.R")

source("mixOmicsDD/R/tune.splsda.R")
source("mixOmicsDD/R/splsda.R")


source("mixOmicsDD/R/predict.mint.block.pls.R")
source("mixOmicsDD/R/internal_predict.DA.R")
source("mixOmicsDD/R/MCVfold.R")

system.time(tune<-tune.splsda(X,Y,test.keepX=1:50, ncomp=2))
plot(tune)

