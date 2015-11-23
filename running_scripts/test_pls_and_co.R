# created on 12/03/15
# Author: F.Rohart
#purpose: test the pls/plsda/spls/splsda function


library(mixOmics)

#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################

# pls
# ----
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, mode = "classic")

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.pls <- pls(X, Y, ncomp = 3)

# plsda
# ----
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)
plotIndiv(plsda.breast, ind.names = TRUE, plot.ellipse = TRUE, add.legend = TRUE)


data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$treatment[, 4]

plsda.liver <- plsda(X, Y, ncomp = 2)
plotIndiv(plsda.liver, ind.names = Y, plot.ellipse = TRUE, add.legend =TRUE)


# spls
# ----
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.spls <- spls(X, Y, ncomp = 2, keepX = c(50, 50),
keepY = c(10, 10))

# splsda
# ----
data(breast.tumors)
X <- breast.tumors$gene.exp
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(breast.tumors$sample$treatment)

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))


# individual names appear
plotIndiv(res, ind.names = Y, add.legend = TRUE, plot.ellipse =TRUE)


data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))

# individual name is set to the treatment
plotIndiv(splsda.liver, ind.names = Y, plot.ellipse = TRUE, add.legend = TRUE)

#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################
if(additional.test==TRUE)
{
    
    set.seed(123)
    
    data(liver.toxicity)
    X=liver.toxicity$gene
    Y=liver.toxicity$clinic
    Z=as.factor(liver.toxicity$treatment[, 4])
    
    # result pls
    test1=pls(X, Y, ncomp = 3)
    
    #test2=pls(X=X, ncomp = 3)       # without Y
    #test2=pls(Y=Y, ncomp = 3)       # without X
    #test2=pls(X=X,Y=Z,ncomp = 3)    # with a factor
    
    
    #result plsda
    test3=plsda(X, Y=Z, ncomp = 2)
    
    #test4=plsda(X=X, ncomp = 2)         #without Y
    #test4=plsda(Y=Z, ncomp = 2)         # without X
    #test4=plsda(X=X, Y=Y, ncomp = 2)    # with a non factor
    
    
    # result spls
    test1=spls(X, Y, ncomp = 2, keepX = c(50, 50),keepY = c(10, 10))
    
    #test2=spls(X=X, ncomp = 2, keepX = c(50, 50),keepY = c(10, 10)) #without Y
    #test3=spls(Y=Y, ncomp = 2, keepX = c(50, 50),keepY = c(10, 10)) #without X
    #test4=spls(X, Y=Z, ncomp = 2, keepX = c(50, 50),keepY = c(10, 10)) # with a factor
    #test5=spls(X, Y, ncomp = 2, keepY = c(10, 10)) # without keepX
    test6=spls(X, Y, ncomp = 2, keepX = c(50, 50)) #without keepY
    test7=spls(X, Y, ncomp = 2) #without keepX/Y
    #test8=spls(X, Y, ncomp = 2, keepX = c(50, 50,40),keepY = c(10, 10)) # with no match between keepX/ncomp
    #test9=spls(X, Y, ncomp = 2, keepX = c(50, 50),keepY = c(10,30)) # with keepY too big
    
    
    # result splsda
    test1=splsda(X, Y=Z, ncomp = 2, keepX = c(20, 20))
    
    #test2=splsda(X=X, ncomp = 2, keepX = c(50, 50)) #without Y
    #test3=splsda(Y=Z, ncomp = 2, keepX = c(50, 50)) #without X
    #test4=splsda(X, Y=Y, ncomp = 2, keepX = c(50, 50)) # without a factor
    #test5=splsda(X, Y=Z, ncomp = 2, keepY = c(10, 10)) # without keepX
    test6=splsda(X, Y=Z, ncomp = 2, keepX = c(50, 50)) #without keepY
    test7=splsda(X, Y=Z, ncomp = 2) #without keepX/Y
    #test8=splsda(X, Y=Z, ncomp = 2, keepX = c(50, 50,40)) # with no match between keepX/ncomp
    #test9=splsda(X, Y=Z, ncomp = 2, keepX = c(50, 50000)) # with keepY too big
    #test10=splsda(X=X, ncomp = 2, keepY = c(50, 50)) #with too many arguments
    
    
    
}