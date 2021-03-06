# created on 12/03/15
# last modified: 02-03-2016
# Author: F.Rohart
#purpose: test the pls/plsda/spls/splsda function


#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

# pls
# ----
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, mode = "classic")


#library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.pls <- pls(X, Y, ncomp = 3)

Y1=Y[,1]
tox=spls(X,Y1)

tox=spls(X,Y[,1,drop=FALSE])
plotVar(tox)


plotIndiv(toxicity.pls, group=rep(1:4,16))
plotIndiv(toxicity.pls, group=factor(rep(1:4,16),labels=c("a","b","c","d")),legend=TRUE)
plotIndiv(toxicity.pls, group=rep(1:4,16),ind.names=FALSE)
plotIndiv(toxicity.pls, group=rep(1:4,16),ind.names=FALSE,legend=TRUE)
plotIndiv(toxicity.pls, group=rep(1:4,16),ellipse=TRUE,legend=TRUE)


plotIndiv(toxicity.pls, group=rep(1:4,16),style="lattice")
plotIndiv(toxicity.pls, group=factor(rep(1:4,16),labels=c("a","b","c","d")),legend=TRUE,style="lattice")
plotIndiv(toxicity.pls, group=rep(1:4,16),ind.names=FALSE,style="lattice")
plotIndiv(toxicity.pls, group=rep(1:4,16),ind.names=FALSE,legend=TRUE,style="lattice")
plotIndiv(toxicity.pls, group=rep(1:4,16),ellipse=TRUE,legend=TRUE,style="lattice")
plotIndiv(toxicity.pls, group=rep(1:4,16),star=TRUE,centroid=TRUE,style="lattice")


# plsda
# ----
#library(mixOmics)
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)
plotIndiv(plsda.breast, ind.names = TRUE, ellipse = TRUE, legend = TRUE)

plotLoadings(plsda.breast)


data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$treatment[, 4]

plsda.liver <- plsda(X, Y, ncomp = 2)
plotIndiv(plsda.liver, ind.names = Y, ellipse = TRUE, legend =TRUE)

plotIndiv(plsda.liver, ind.names = FALSE, ellipse = TRUE, legend =TRUE, centroids = TRUE)


# spls
# ----
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.spls <- spls(X, Y, ncomp = 2, keepX = c(50, 50),
keepY = c(10, 10))
plotIndiv(toxicity.spls)

# splsda
# ----
#library(mixOmics)
data(breast.tumors)
X <- breast.tumors$gene.exp
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(breast.tumors$sample$treatment)

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))

plotLoadings(res)
plotLoadings(res, contrib = "min")
plotLoadings(res, contrib = "max")

tune= tune.splsda(X,Y,ncomp=2,nrepeat=5,logratio="none",test.keepX = c(5, 10, 15),folds=10,dist="max.dist")


# individual names appear
plotIndiv(res, ind.names = Y, legend = TRUE, ellipse =TRUE)

#with only 1 sample in a class
Y <- as.factor(breast.tumors$sample$treatment)
ind.remove = which(Y==levels(Y)[1])[-sample(1:sum(Y==levels(Y)[1]),1)]
Y2 = Y[-ind.remove]
X2 = X[-ind.remove, ]

res2 <- splsda(X2, Y2, ncomp = 2, keepX = c(25, 25))
plotIndiv(res2, ind.names = Y2, legend = TRUE, ellipse =TRUE)

#with only 1 sample in a class
Y <- as.factor(breast.tumors$sample$treatment)
ind.remove = which(Y==levels(Y)[2])[-sample(1:sum(Y==levels(Y)[2]),1)]
Y3 = Y[-ind.remove]
X3 = X[-ind.remove, ]

res3 <- splsda(X3, Y3, ncomp = 2, keepX = c(25, 25))
plotIndiv(res3, ind.names = Y3, legend = TRUE, ellipse =TRUE)


data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))

# individual name is set to the treatment
plotIndiv(splsda.liver, ind.names = Y, ellipse = TRUE, legend = TRUE)

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
    
    
    
    #######  wraper.pls
    res=pls(X,Y)
    res=pls(X,Y,ncomp=3)
    res=pls(X,Y,ncomp=3)
    res=pls(X,Y,near.zero.var=TRUE)
    
    
    
    #######  wraper.spls
    res=spls(X,Y)
    res=spls(X,Y,ncomp=3,keepX=c(10,5,15))
    res=spls(X,Y,ncomp=3,keepX=c(10,5)) #complete keepX
    res=spls(X,Y,ncomp=3,keepX=c(10))#,keepX.constraint=list(comp1=c(100,1,3),comp2=c(10)))
    res=spls(X,Y,ncomp=3,keepX=c(1))#,keepX.constraint=list(comp1=c(100),comp2=c(10)))
    res=spls(X,Y,ncomp=3)#,keepX.constraint=list(comp1=c(100),comp2=c(10)))
    res=spls(X,Y,ncomp=3,keepY=c(3))#,keepY.constraint=list(comp1=c(1),comp2=c(2)))
    
    #######  wraper.plsda
    res=plsda(X,Z,ncomp=3)
    
    
    #######  wraper.splsda
    res=splsda(X,Z)
    res=splsda(X,Z,ncomp=3,keepX=c(10,5,15))
    res=splsda(X,Z,ncomp=3,keepX=c(10,5)) #complete keepX
    res=splsda(X,Z,ncomp=3,keepX=c(10))#,keepX.constraint=list(comp1=c(100,1,3),comp2=c(10)))
    res=splsda(X,Z,ncomp=3,keepX=c(1))#,keepX.constraint=list(comp1=c(100),comp2=c(10)))
    res=splsda(X,Z,ncomp=3)#,keepX.constraint=list(comp1=c(100),comp2=c(10)))
    

    
    
}
par(opar)
