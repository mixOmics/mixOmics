library(mixOmics)
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
legend('topright', col = c('red', 'darkgreen'), legend = c('Q2 total', 'threshold 0.0975') , lty = 1)
title('Liver toxicity PLS 5-fold, Q2 values')
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
keepX = c(rep(100, ncomp)), keepY = c(rep(4,ncomp)))
# with leave-one-out cross validation
set.seed(45)
model.spls.loo.val <- perf(model.spls, validation = "Mfold",folds=5)
#Q2 total
model.spls.loo.val$Q2.total

plot(model.spls.loo.val, criterion="R2")
# R2:we can see how the performance degrades when ncomp increases
# results are similar to 5-fold
model.spls.loo.val$R2
## validation for objects of class 'splsda' (classification)
# ----------------------------------------
library(mixOmics)
data(srbct)
X <- srbct$gene
Y <- srbct$class
ncomp = 5
srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))
# with Mfold
# ---------
set.seed(45)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8,
method.predict = "all")
plot(error, criterion="R2",type = "l",pred.method="all")
## End(Not run)

