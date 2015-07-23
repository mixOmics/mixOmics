# ----------------------------------------
# testing orkshop material
# date: 23/03/2015
# latest update: 23/07/2015 for V5.1 update
# ----------------------------------------

# setwd("~/Documents/k.lecao/Packages/mixOmics/GIT/package-mixomics")
# 
# # ------ notes for me to compile the package (if need be) on a terminal
# R CMD build --resave-data mixOmics
# R CMD INSTALL -l MyR/ mixOmics_5.0-4.tar.gz
# R CMD check mixOmics --as-cran --timings
# # ---------------------------------
# 
# 
# # now in R, load the package
# detach("package:mixOmics", unload=TRUE)
# 
# library(mixOmics, lib.loc = 'MyR/')
# # check that version 5.0-4 is loaded
# sessionInfo() # ok mixOmics 5.0-4


library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# KA
sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside mixOmics/R


# ---------------------------------------
# help files first
# ---------------------------------------
## validation for objects of class 'pls' (regression)
# ----------------------------------------
\dontrun{
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
  legend('topright', col = c('red', 'darkgreen'), legend = c('Q2 total', 'threshold 0.0975'), lty = 1)
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
  error <- perf(srbct.splsda, validation = "Mfold", folds = 8, 
                method.predict = "all")
  
  plot(error, type = "l")
}



###################################################
### code chunk number 19: meth_PLS.Rnw:102-108
###################################################
data(nutrimouse)
X <- nutrimouse$gene
Y <- nutrimouse$lipid
# checking IDs are matching
#cbind(rownames(X), rownames(Y))
res.pls1 <- pls(X, Y, ncomp = 10, mode = "regression")


###################################################
### code chunk number 20: meth_PLS.Rnw:115-118
###################################################

tune.pls <- perf(res.pls1, validation = 'loo', #validation = 'Mfold', folds = 10,
                 criterion = 'all', progressBar = FALSE)

###################################################
### code chunk number 21: Q2.total
###################################################
plot(tune.pls$Q2.total)
abline(h = 0.0975)


###################################################
### code chunk number 22: meth_PLS.Rnw:133-136
###################################################
res.pls2 <- pls(X, Y, ncomp = 5, mode = 'regression')
perf.pls <- perf(res.pls2, validation = 'loo', #'Mfold', folds = 10,
                 criterion = 'all', progressBar = FALSE)


###################################################
### code chunk number 23: meth_PLS.Rnw:141-142
###################################################
head(perf.pls$MSEP)


###################################################
### code chunk number 24: meth_PLS.Rnw:147-148
###################################################
head(perf.pls$R2)


###################################################
### code chunk number 25: meth_PLS.Rnw:154-162
###################################################
ncomp = 10
res.spls <- spls(X, Y, ncomp = ncomp, keepX = c(rep(10, ncomp)), keepY = c(rep(5, ncomp)),
                 mode = 'regression')
perf.spls <- perf(res.spls, validation = 'loo',#, folds = 10,
                  criterion = 'all', progressBar = FALSE)

# !! changed
plot(perf.spls$Q2.total)
abline(h = 0.0975)

head(perf.spls$MSEP)
head(perf.spls$R2)


###################################################
### code chunk number 26: meth_PLS.Rnw:166-167 (eval = FALSE)
###################################################
#plot(perf.spls, criterion = c('RMSEP'), type = 'l')





###################################################
### code chunk number 29: meth_PLSDA.Rnw:99-105
###################################################
##library(mixOmics)
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
# here, Y is a factor
Y <- as.factor(liver.toxicity$treatment[, "Dose.Group"]) 
summary(Y)


###################################################
### code chunk number 30: meth_PLSDA.Rnw:112-121
###################################################
res.plsda <- plsda(X, Y, ncomp = 10)
set.seed(45)
tune.plsda <- perf(res.plsda, method.predict = 'all', 
                   validation = 'Mfold', folds = 5, 
                   progressBar = FALSE) 

# !note that the function perf should be rerun several times 
# and the error rates obtained below should then be averaged
tune.plsda$error.rate


###################################################
### code chunk number 31: meth_PLSDA.Rnw:127-132
###################################################

# !! changed
matplot(tune.plsda$error.rate, type = 'l', lty = 1, lwd = 2, 
        xlab = 'ncomp', ylab = 'Classification error rate with PLSDA')
legend('topright', col = c(1:3),  
       legend = c('max.dist', 'centroid.dist', 'Mahalanobis.dist'), 
       lty = 1, lwd = 1)

# might need to save plot?
plot(tune.plsda, type = 'l')






###################################################
### code chunk number 96: cs_liver_toxicity.Rnw:22-29
###################################################
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

# ! always check that the subjects are matched in the two data sets
# just a head
head(cbind(rownames(X), rownames(Y)))


###################################################
### code chunk number 97: cs_liver_toxicity.Rnw:35-40
###################################################
#!! mfold changed
res.pls1 <- pls(X, Y, ncomp = 10, mode = 'regression')
tune.pls <- perf(res.pls1, validation = 'Mfold', folds = 10, 
                 criterion = 'all', progressBar = FALSE) 
plot(tune.pls$Q2.total)
abline(h = 0.0975)


###################################################
### code chunk number 98: cs_liver_toxicity.Rnw:47-50
###################################################
liver.pls <- pls(X, Y, ncomp = 3, mode = "regression")
liver.spls <- spls(X, Y, ncomp = 3, keepX = c(10, 5, 10), 
                   mode = "regression")


###################################################
### code chunk number 99: cs_liver_toxicity.Rnw:59-81
###################################################
M <- 10 # number of folds
perf.pls <- perf(liver.pls, validation = "Mfold", folds = M, 
                 progressBar = FALSE)
perf.spls <- perf(liver.spls, validation = "Mfold", folds = M, 
                  progressBar = FALSE)

plot(perf.spls$Q2.total)
abline(h = 0.0975)

plot(perf.pls$Q2.total)
abline(h = 0.0975)

perf.pls$Q2.total
perf.spls$Q2.total

##plot(perf.pls, criterion = 'MSEP', type = 'l')

## ! changed
par(mfrow = c(3, 3))
for(i in 1:9) {
  spls.rmsep <- sqrt(perf.spls$MSEP[i, ])
  pls.rmsep <- sqrt(perf.pls$MSEP[i, ])
  matplot(cbind(spls.rmsep, pls.rmsep), lwd = 2, xlab = "dim", 
          ylab = "MSEP", type = "l", lty = c(1, 2), 
          axes = FALSE)
  axis(1, 1:3, labels = 1:3)
  axis(2)
  title(main = paste(rownames(perf.spls$MSEP)[i]))
}
par(mfrow = c(1, 1))


## this is to checl V5.1, making sure this output does not output 'others'
perf.spls$features







