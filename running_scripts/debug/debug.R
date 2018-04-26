#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
library(rARPACK)
library(matrixStats)
library(mixOmics)

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
Y = nutrimouse$diet
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
set.seed(12)
# add missing values in the data
data[[1]][sample(1:nrow(data[[1]]),10), sample(1:ncol(data[[1]]),10)] = NA
data[[2]][sample(1:nrow(data[[2]]),10), sample(1:ncol(data[[2]]),10)] = NA

# adding constant variables and almost constant variables
mat = matrix(1:3,nrow=nrow(data[[1]]),ncol=6,byrow=T)
mat[sample(1:nrow(mat),4),4:6] = rnorm(4)
data[[1]] = cbind(data[[1]], mat)
head(data[[1]])


set.seed(45)
# classic tune
tune = tune.block.splsda(
X = data,
Y = Y,
ncomp = 2,
design = design,
folds=5,
test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),
)

nutrimouse.sgccda <- wrapper.sgccda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10), lipid=c(15)),
ncomp = 3,#c(2, 2),
scheme = "centroid",
tol=1e-30)

set.seed(12)
a=perf(nutrimouse.sgccda,validation = "Mfold", folds = 10,
nrepeat = 10, progressBar = TRUE)



# add missing values in the data
data[[1]][sample(1:nrow(data[[1]]),2), sample(1:ncol(data[[1]]),2)] = NA


nutrimouse.sgccda <- wrapper.sgccda(X=data, Y=Y
design = design,
keepX = list(gene=c(10), lipid=c(15)),
ncomp = 3,#c(2, 2),
scheme = "centroid",
tol=1e-30)

source("mixOmics/R/internal_wrapper.mint.R")
source("mixOmics/R/internal_wrapper.mint.block.R")
source("mixOmics/R/internal_mint.block.R")
source("mixOmics/R/internal_mint.block_helpers.R")
source("mixOmics/R/check_entry.R")
library(rARPACK)
library(matrixStats)



source("mixOmics/R/internal_wrapper.mint.R")
source("mixOmics/R/internal_wrapper.mint.block.R")
source("mixOmics/R/internal_mint.block.R")
source("mixOmics/R/internal_mint.block_helpers.R")
source("mixOmics/R/check_entry.R")
source("mixOmics/R/t.test.process.R")

#source("mixOmicsDD/R/tune.splsda.R")
source("mixOmics/R/tune.splsda.R")
source("mixOmics/R/perf.R")
source("mixOmics/R/MCV.splsda.R")
source("mixOmics/R/tune.spls.R")
source("mixOmics/R/MCV.spls.R")
source("mixOmics/R/internal_graphic.perf.R")

source("mixOmics/R/nipals.R")


data(srbct)
X <- srbct$gene
Y <- srbct$class
ncomp=3

srbct.plsda <- plsda(X, Y, ncomp = 10)

set.seed(1111)
perf.srbct.plsda <- perf(srbct.plsda, validation = "Mfold", folds = 5,
                   progressBar = TRUE, nrepeat = 10)
perf.srbct.plsda$choice.ncomp
#Error in ncomp_opt[measure, ijk] <- t.test.process(t(mat.error.rate[[measure_i]][[ijk]])) :
#number of items to replace is not a multiple of replacement length
set.seed(2222)
perf.srbct.plsda <- perf(srbct.plsda, validation = "Mfold", folds = 5,
                   progressBar = TRUE, nrepeat = 10)
perf.srbct.plsda$choice.ncomp
#Error in ncomp_opt[measure, ijk] <- t.test.process(t(mat.error.rate[[measure_i]][[ijk]])) :
#number of items to replace is not a multiple of replacement length

set.seed(2543)
perf.srbct.plsda <- perf(srbct.plsda, validation = "Mfold", folds = 5,
                   progressBar = TRUE, nrepeat = 10)
perf.srbct.plsda$choice.ncomp



data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.pls <- pls(X, Y, ncomp = 3)

Y=Y[,9]
tox=spls(X,Y, keepX=c(1,10))
p=predict(tox,X)
plot(Y,p$predict[,,1])
sum((Y-p$predict[,,1])^2)/length(Y)/var(Y)

tox=spls(X,Y, keepX=c(10,10))
p=predict(tox,X)
plot(Y,p$predict[,,1])
sum((Y-p$predict[,,1])^2)/length(Y)/var(Y)
sum((Y-p$predict[,,2])^2)/length(Y)/var(Y)



tox=spls(X,Y, keepX=c(10,10))

selectVar(tox)





source("mixOmics/R/internal_wrapper.mint.R")
source("mixOmics/R/internal_wrapper.mint.block.R")
source("mixOmics/R/internal_mint.block.R")
source("mixOmics/R/internal_mint.block_helpers.R")
source("mixOmics/R/check_entry.R")
source("mixOmics/R/t.test.process.R")

#source("mixOmicsDD/R/tune.splsda.R")
source("mixOmics/R/tune.spls.R")
source("mixOmics/R/MCV.spls.R")
source("mixOmics/R/internal_graphic.perf.R")


data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic


tun = tune.spls(X,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=5)

par(mfrow=c(2,2))
set.seed(12)
tun = tune.spls(X,Y,ncomp=2, test.keepX = c(5,10,15), nrepeat=5, measure = "MSE")
plot(tun)
set.seed(12)
tun2 = tune.spls(X,Y,ncomp=2, test.keepX = c(5,10,15), nrepeat=5, measure = "MAE")
plot(tun2)
set.seed(12)
tun3 = tune.spls(X,Y,ncomp=2, test.keepX = c(5,10,15), nrepeat=5, measure = "Bias")
plot(tun3)


tun = tune.spls(X,Y,ncomp=5, test.keepX = c(1:10,15,seq(20,200,10),250,seq(300,3100,100),3116), nrepeat=10)

plot(tun)
plot(tun,sd=FALSE)


par(mfrow=c(3,2))
set.seed(12)
tun = tune.spls(X,Y,ncomp=5, test.keepX = c(1:10,15,seq(20,200,10),250,seq(300,3100,100),3116), nrepeat=10, measure = "MSE")
plot(tun)
plot(tun,sd=FALSE)
set.seed(12)
tun2 = tune.spls(X,Y,ncomp=5, test.keepX = c(1:10,15,seq(20,200,10),250,seq(300,3100,100),3116), nrepeat=10, measure = "MAE")
plot(tun2)
plot(tun2,sd=FALSE)
set.seed(12)
tun3 = tune.spls(X,Y,ncomp=5, test.keepX = c(1:10,15,seq(20,200,10),250,seq(300,3100,100),3116), nrepeat=10, measure = "Bias")
plot(tun3)
plot(tun3,sd=FALSE)






# gives same results )where X.test is not scaled
object=spls(X.train,Y.train.mat, ncomp=2, keepX=c(10,10))
p=predict(object, X.test)




load("temp.Rdata")


source("mixOmics/R/block.splsda.R")
source("mixOmics/R/block.spls.R")
source("mixOmics/R/splsda.R")
source("mixOmics/R/spls.R")
source("mixOmics/R/pls.R")
source("mixOmics/R/plsda.R")
source("mixOmics/R/MCV.spls.R")
source("mixOmics/R/MCV.splsda.R")
source("mixOmics/R/perf.R")
source("mixOmics/R/tune.block.splsda.R")
source("mixOmics/R/MCV.block.splsda.R")
source("mixOmics/R/internal_wrapper.mint.block.R")
source("mixOmics/R/internal_mint.block.R")
source("mixOmics/R/predict.mint.block.pls.R")
source("mixOmics/R/internal_predict.DA.R")


nutrimouse.sgccda <- block.splsda(X=data,
Y = Y,
design = design,
#keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,#c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE,
tol=1e-30)




set.seed(12)
tune2 = tune.block.splsda(X = data,Y = Y,design = design,ncomp = 2,test.keepX = list(gene=c(1,5,10,4),lipid=c(1,2,3)),scheme = "centroid",bias = FALSE,weighted=FALSE, measure="overall", init="svd")


tune2
}
