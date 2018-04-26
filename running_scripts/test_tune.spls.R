#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
library(rARPACK)
library(matrixStats)
library(mixOmics)


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




if(FALSE)
{
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
