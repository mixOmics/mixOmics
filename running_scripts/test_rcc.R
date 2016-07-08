#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.res <- rcc(X, Y)

## Regularized CCA
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res1 <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

## using shrinkage parameters
nutri.res2 <- rcc(X, Y, ncomp = 3, method = 'shrinkage')
nutri.res2$lambda # the shrinkage parameters


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{

    data(linnerud)
    X <- linnerud$exercise
    Y <- linnerud$physiological
    linn.res <- rcc(X, Y)

    #test ncomp
    linn.res <- rcc(X, Y, ncomp = 1)
    #plotIndiv(linn.res) # can not plot with only 1 comp
    linn.res <- rcc(X, Y, ncomp = 2)
    plotIndiv(linn.res,comp=c(1,2))
    linn.res <- rcc(X, Y, ncomp = 3)
    plotIndiv(linn.res,comp=c(1,3))

    #rcc(X, Y, ncomp = 4)

    #test method
    linn.res <- rcc(X, Y, method = "ridge")
    plotIndiv(linn.res)
    linn.res <- rcc(X, Y, method = "shrinkage")
    plotIndiv(linn.res)
    linn.res <- rcc(X, Y,ncomp=3, method = "shrinkage")
    plotIndiv(linn.res,comp=c(1,3))

    # linn.res <- rcc(X, Y,ncomp=3, method = "alpha") # error method


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

    #linn.res <- rcc(X, Y, lambda1=-0.01,lambda2=1)

}
par(opar)
