#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
progressBar = FALSE

#library(rARPACK)
#library(matrixStats)
#library(mixOmics)


data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.matrix(liver.toxicity$clinic)

toxicity.pls <- pls(X, Y, ncomp = 3)

Y=Y[,9,drop=FALSE]
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



data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.matrix(liver.toxicity$clinic)

tun = tune.spls(X,Y,ncomp=4, test.keepX = c(5,10,15), nrepeat=5, test.keepY=c(2,5,8,10), progressBar = progressBar)
plot(tun,keepY=5)
plot(tun,keepY=10)

Y=Y[,9,drop=FALSE]
tun = tune.spls(X,Y,ncomp=1, test.keepX = c(5,10,15), nrepeat=5, progressBar = progressBar)

par(mfrow=c(2,2))
set.seed(12)
tun = tune.spls(X,Y,ncomp=2, test.keepX = c(5,10,15), nrepeat=5, measure = "MSE", progressBar = progressBar)
plot(tun)
set.seed(12)
tun2 = tune.spls(X,Y,ncomp=2, test.keepX = c(5,10,15), nrepeat=5, measure = "MAE", progressBar = progressBar)
plot(tun2)
set.seed(12)
tun3 = tune.spls(X,Y,ncomp=2, test.keepX = c(5,10,15), nrepeat=5, measure = "Bias", progressBar = progressBar)
plot(tun3)
set.seed(12)
tun4 = tune.spls(X,Y,ncomp=2, test.keepX = c(5,10,15), nrepeat=5, measure = "R2", progressBar = progressBar)
plot(tun4)



#test.keepX = c(1:10,15,seq(20,200,10),250,seq(300,3100,100),3116)
test.keepX = c(5,10,15,100,200,300,1000)

par(mfrow=c(2,2))
set.seed(12)
tun = tune.spls(X,Y,ncomp=5, test.keepX = test.keepX, nrepeat=10, measure = "MSE", progressBar = progressBar)
plot(tun,legend.position = "topleft")
plot(tun,sd=FALSE,legend.position = "topleft")
set.seed(12)
tun2 = tune.spls(X,Y,ncomp=5, test.keepX = test.keepX, nrepeat=10, measure = "MAE", progressBar = progressBar)
plot(tun2,legend.position = "topleft")
plot(tun2,sd=FALSE,legend.position = "topleft")
set.seed(12)
tun3 = tune.spls(X,Y,ncomp=5, test.keepX = test.keepX, nrepeat=10, measure = "Bias", progressBar = progressBar)
plot(tun3,legend.position = "bottomright")
plot(tun3,sd=FALSE,legend.position = "bottomright")
set.seed(12)
tun4 = tune.spls(X,Y,ncomp=5, test.keepX = test.keepX, nrepeat=10, measure = "R2", progressBar = progressBar)
plot(tun4,legend.position = "bottomleft")
plot(tun4,sd=FALSE,legend.position = "bottomleft")


### with NA
n=nrow(X)
p=ncol(X)
system.time(tun <- tune.spls(X,Y,ncomp=3, test.keepX = c(5,10,15), nrepeat=1, test.keepY=(5:10), folds=3, progressBar = progressBar))

X[sample(1:(n*p),100)]=NA
system.time(tun2 <- tune.spls(X,Y,ncomp=3, test.keepX = c(5,10,15), nrepeat=1, test.keepY=(5:10), folds=3, progressBar = progressBar))

