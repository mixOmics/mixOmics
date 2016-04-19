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


# splsda
# ----
data(breast.tumors)
X <- breast.tumors$gene.exp
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(breast.tumors$sample$treatment)

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))


source("mixOmics/R/tune.splsda.R")
source("mixOmics/R/MCVfold.R")

tune= tune.splsda(X,Y,ncomp=1,nrepeat=1,logratio="none",test.keepX = c(5),folds=10,dist="max.dist",already.tested.X=NULL, progressBar = TRUE)



tune= tune.splsda(X,Y,ncomp=1,nrepeat=10,logratio="none",test.keepX = c(5, 10, 15),folds=10,dist="max.dist",already.tested.X=NULL, progressBar = TRUE)


tune= tune.splsda(X,Y,ncomp=3,nrepeat=10,logratio="none",test.keepX = c(5, 10, 15),folds=10,dist="max.dist",already.tested.X=NULL, progressBar = TRUE)

#source("mixOmics/R/tune.splsda.R")
#source("mixOmics/R/MCVfold.R")
tune= tune.splsda(X,Y,ncomp=3,nrepeat=10,logratio="none",test.keepX = c(5,50,100),folds=10,dist="max.dist",already.tested.X=NULL, progressBar = TRUE)


par(opar)
