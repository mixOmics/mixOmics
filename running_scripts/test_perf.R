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


data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic


# try tune the number of component to choose
# ---------------------
# first learn the full model
liver.pls <- pls(X, Y, ncomp = 10)

# with 5-fold cross validation: we use the same parameters as in model above
# but we perform cross validation to compute the MSEP, Q2 and R2 criteria
# ---------------------------
liver.val <- perf(liver.pls, validation = "Mfold", folds = 5)

# Q2 total should decrease until it reaches a threshold
liver.val$Q2.total

# ncomp = 2 is enough
plot(liver.val$Q2.total, type = 'l', col = 'red', ylim = c(-0.5, 0.5),
xlab = 'PLS components', ylab = 'Q2 total')
abline(h = 0.0975, col = 'darkgreen')
legend('topright', col = c('red', 'darkgreen'),
legend = c('Q2 total', 'threshold 0.0975'), lty = 1)
title('Liver toxicity PLS 5-fold, Q2 total values')

#have a look at the other criteria
# ----------------------
# R2
liver.val$R2
matplot(t(liver.val$R2), type = 'l', xlab = 'PLS components', ylab = 'R2 for each variable')
title('Liver toxicity PLS 5-fold, R2 values')

# MSEP
liver.val$MSEP
matplot(t(liver.val$MSEP), type = 'l', xlab = 'PLS components', ylab = 'MSEP for each variable')
title('Liver toxicity PLS 5-fold, MSEP values')


## validation for objects of class 'spls' (regression)
# ----------------------------------------
ncomp = 7
# first, learn the model on the whole data set
model.spls = spls(X, Y, ncomp = ncomp, mode = 'regression',
keepX = c(rep(10, ncomp)), keepY = c(rep(4,ncomp)))


# with leave-one-out cross validation
##set.seed(45)
model.spls.val <- perf(model.spls, validation = "Mfold", folds = 5 )#validation = "loo")

#Q2 total
model.spls.val$Q2.total

# R2:we can see how the performance degrades when ncomp increases
model.spls.val$R2
plot(model.spls.val, criterion="R2", type = 'l')
plot(model.spls.val, criterion="Q2", type = 'l')


## validation for objects of class 'plsda' (classification)
# ----------------------------------------
data(srbct)
X <- srbct$gene
Y <- srbct$class

ncomp = 5

srbct.plsda <- plsda(X, Y, ncomp = ncomp)

# with Mfold
# ---------
set.seed(45)
error <- perf(srbct.plsda, validation = "Mfold", folds = 8, dist = "all", progressBar = FALSE)

plot(error)

error <- perf(srbct.plsda, validation = "Mfold", folds = 8, dist = "all", progressBar = FALSE, constraint = TRUE)
#warning


#source("mixOmics/R/perf.R")
#source("mixOmics/R/MCVfold.R")
error.overall <- perf(srbct.plsda, validation = "Mfold", folds = 8, dist = "all", progressBar = FALSE)
error.BER <- perf(srbct.plsda, validation = "Mfold", folds = 8, dist = "all", progressBar = FALSE)
plot(error)

#source("mixOmics/R/perf.R")
#source("mixOmics/R/MCVfold.R")
error.both <- perf(srbct.plsda, validation = "Mfold", folds = 8, dist = "all", progressBar = FALSE)

plot(error.both)
plot(error.both, measure = "BER")



## validation for objects of class 'splsda' (classification)
# ----------------------------------------
data(srbct)
X <- srbct$gene
Y <- srbct$class

ncomp = 5

srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))

# with Mfold
# ---------
set.seed(45)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8, dist = "all", progressBar = FALSE)

plot(error )
plot(error.both , measure = "BER")

error <- perf(srbct.splsda, validation = "Mfold", folds = 8, dist = "all", progressBar = FALSE, constraint = TRUE)



par(opar)
