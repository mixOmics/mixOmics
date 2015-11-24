# created on 12/03/15
# Author: F.Rohart
#purpose: test the pls/plsda/spls/splsda function


set.seed(123)
source("package-mixomics/mixOmics/R/pls.R")
source("package-mixomics/mixOmics/R/plsda.R")
source("package-mixomics/mixOmics/R/spls.R")
source("package-mixomics/mixOmics/R/splsda.R")


data(liver.toxicity)
X=liver.toxicity$gene
Y=liver.toxicity$clinic
Z=as.factor(liver.toxicity$treatment[, 4])

# result pls
test1=pls(X, Y, ncomp = 3)


test2=pls(X=X, ncomp = 3)       # without Y
test2=pls(Y=Y, ncomp = 3)       # without X
test2=pls(X=X,Y=Z,ncomp = 3)    # with a factor


#result plsda
test3=plsda(X, Y=Z, ncomp = 2)

test4=plsda(X=X, ncomp = 2)         #without Y
test4=plsda(Y=Z, ncomp = 2)         # without X
test4=plsda(X=X, Y=Y, ncomp = 2)    # with a non factor


# result spls
test1=spls(X, Y, ncomp = 2, keepX = c(50, 50),keepY = c(10, 10))

test2=spls(X=X, ncomp = 2, keepX = c(50, 50),keepY = c(10, 10)) #without Y
test3=spls(Y=Y, ncomp = 2, keepX = c(50, 50),keepY = c(10, 10)) #without X
test4=spls(X, Y=Z, ncomp = 2, keepX = c(50, 50),keepY = c(10, 10)) # with a factor
test5=spls(X, Y, ncomp = 2, keepY = c(10, 10)) # without keepX
test6=spls(X, Y, ncomp = 2, keepX = c(50, 50)) #without keepY
test7=spls(X, Y, ncomp = 2) #without keepX/Y
test8=spls(X, Y, ncomp = 2, keepX = c(50, 50,40),keepY = c(10, 10)) # with no match between keepX/ncomp
test9=spls(X, Y, ncomp = 2, keepX = c(50, 50),keepY = c(10,30)) # with keepY too big


# result splsda
test1=splsda(X, Y=Z, ncomp = 2, keepX = c(20, 20))

test2=splsda(X=X, ncomp = 2, keepX = c(50, 50)) #without Y
test3=splsda(Y=Z, ncomp = 2, keepX = c(50, 50)) #without X
test4=splsda(X, Y=Y, ncomp = 2, keepX = c(50, 50)) # without a factor
test5=splsda(X, Y=Z, ncomp = 2, keepY = c(10, 10)) # without keepX
test6=splsda(X, Y=Z, ncomp = 2, keepX = c(50, 50)) #without keepY
test7=splsda(X, Y=Z, ncomp = 2) #without keepX/Y
test8=splsda(X, Y=Z, ncomp = 2, keepX = c(50, 50,40)) # with no match between keepX/ncomp
test9=splsda(X, Y=Z, ncomp = 2, keepX = c(50, 50000)) # with keepY too big
test10=splsda(X=X, ncomp = 2, keepY = c(50, 50)) #with too many arguments



