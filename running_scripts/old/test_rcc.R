require(corpcor)

data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.res <- rcc(X, Y)

#test ncomp
linn.res <- rcc(X, Y, ncomp = 1)
plotIndiv(linn.res)
linn.res <- rcc(X, Y, ncomp = 2)
plotIndiv(linn.res,comp=c(1,2))
linn.res <- rcc(X, Y, ncomp = 3)
plotIndiv(linn.res,comp=c(1,3))

rcc(X, Y, ncomp = 4)

#test method
linn.res <- rcc(X, Y, method = "ridge")
plotIndiv(linn.res)
linn.res <- rcc(X, Y, method = "shrinkage")
plotIndiv(linn.res)
linn.res <- rcc(X, Y,ncomp=3, method = "shrinkage")
plotIndiv(linn.res,comp=c(1,3))

linn.res <- rcc(X, Y,ncomp=3, method = "alpha")


#test lambda
linn.res <- rcc(X, Y, lambda1=1)
plotIndiv(linn.res)
linn.res <- rcc(X, Y, lambda1=0.1)
plotIndiv(linn.res)
linn.res <- rcc(X, Y, lambda1=0.01)
plotIndiv(linn.res)
linn.res <- rcc(X, Y, lambda2=1)
plotIndiv(linn.res)
linn.res <- rcc(X, Y, lambda2=0.1)
plotIndiv(linn.res)
linn.res <- rcc(X, Y, lambda2=8)
plotIndiv(linn.res)
linn.res <- rcc(X, Y, lambda1=1,lambda2=0.01)
plotIndiv(linn.res)
linn.res <- rcc(X, Y, lambda1=0.1,lambda2=0.1)
plotIndiv(linn.res)
linn.res <- rcc(X, Y, lambda1=0.01,lambda2=1)
plotIndiv(linn.res)

linn.res <- rcc(X, Y, lambda1=-0.01,lambda2=1)


