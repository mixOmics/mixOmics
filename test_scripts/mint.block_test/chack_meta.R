###############################
############# PLS #############
###############################
#### Difference because of not the same starting point between the 2 algo "pls" and "spls" in mixOmics package when tol = 1e-06
rm(list=ls())

source("/Users/bgautier/Dropbox/mixOmics/Updates_latest/Francois_parallel/paral.options.R")
source("pls.sgcca.R"); source("predict.sgcca.R"); source("predict.sgcca.R"); source("mix.options.R")
source("perf.sgcca.R"); source("helpers.R"); 
library(mixOmics);

### Example 1 ###
data(linnerud)                            
X <- linnerud$exercise
Y <- linnerud$physiological

### Classic
classic.linn.pls <- pls(X, Y, mode = "regression", tol = 1e-25, ncomp = 3)
classic.linn.pls$loadings$X

classic.linn.pls.sgcca <- pls.sgcca(X, Y, mode = "regression", tol = 1e-25, ncomp = 3)
classic.linn.pls.sgcca$loadings$X

checklist = names(classic.linn.pls)
checklist = checklist[!checklist %in% c("iter")]

classic.linn.pls.abs = abslist(classic.linn.pls, c("variates", "loadings", "mat.c", "mat.d"))
classic.linn.pls.sgcca.abs = abslist(classic.linn.pls.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(classic.linn.pls.abs[names(classic.linn.pls.abs) == x], 
                                        classic.linn.pls.sgcca.abs[names(classic.linn.pls.sgcca.abs) == x])})

head(classic.linn.pls.abs$loadings$X); head(classic.linn.pls.sgcca.abs$loadings$X);

indiv1 <- c(200, 40, 60); indiv2 <- c(190, 45, 45)
newdata <- rbind(indiv1, indiv2)
colnames(newdata) <- colnames(X)
newdata

pred.classic.linn.pls <- predict(classic.linn.pls, newdata)
pred.classic.linn.pls.sgcca <- unlist(predict.sgcca(classic.linn.pls.sgcca, newdata), recursive = FALSE)

names(pred.classic.linn.pls)

max(abs(pred.classic.linn.pls$predict) - abs(pred.classic.linn.pls.sgcca$predict))
max(abs(pred.classic.linn.pls$variates) - abs(pred.classic.linn.pls.sgcca$variates))
max(abs(pred.classic.linn.pls$B.hat) - abs(pred.classic.linn.pls.sgcca$B.hat))

system.time(perf(classic.linn.pls, method.predict = "all", validation = "loo", progressBar = FALSE))
system.time(perf.sgcca(classic.linn.pls.sgcca, validation = "loo"))

classic.linn.pls.perf = perf(classic.linn.pls, method.predict = "all", validation = "loo")
classic.linn.pls.perf.sgcca = perf.sgcca(classic.linn.pls.sgcca, validation = "loo")

names(classic.linn.pls.perf)
all.equal(classic.linn.pls.perf$MSEP, classic.linn.pls.perf.sgcca$MSEP$X.1)
all.equal(classic.linn.pls.perf$R2, classic.linn.pls.perf.sgcca$R2$X.1)
all.equal(classic.linn.pls.perf$Q2, classic.linn.pls.perf.sgcca$Q2$X.1) # Difference
all.equal(classic.linn.pls.perf$Q2.total, classic.linn.pls.perf.sgcca$Q2.total$X.1) # Difference
max(abs(classic.linn.pls.perf$press.mat[, , 1] - classic.linn.pls.perf.sgcca$press.mat$X.1[[1]]))
max(abs(classic.linn.pls.perf$press.mat[, , 2] - classic.linn.pls.perf.sgcca$press.mat$X.1[[2]]))
max(abs(classic.linn.pls.perf$RSS.indiv[, , 2] - classic.linn.pls.perf.sgcca$press.mat$X.1[[1]]))
max(abs(classic.linn.pls.perf$RSS.indiv[, , 3] - classic.linn.pls.perf.sgcca$press.mat$X.1[[2]]))
max(abs(classic.linn.pls.perf$PRESS.inside - classic.linn.pls.perf.sgcca$PRESS.inside$X.1))
max(abs(classic.linn.pls.perf$RSS - classic.linn.pls.perf.sgcca$RSS$X.1)) # Difference

classic.linn.pls.perf$Q2

classic.linn.pls.perf.sgcca$Q2

### Invariant
invariant.linn.pls <- pls(X, Y, mode = "invariant", tol = 1e-25)
invariant.linn.pls.sgcca <- pls.sgcca(X, Y, mode = "invariant", tol = 1e-25)

invariant.linn.pls.abs = abslist(invariant.linn.pls, c("variates", "loadings", "mat.c"))
invariant.linn.pls.sgcca.abs = abslist(invariant.linn.pls.sgcca, c("variates", "loadings", "mat.c"))

lapply(checklist, function(x){all.equal(invariant.linn.pls.abs[names(invariant.linn.pls.abs) == x], invariant.linn.pls.sgcca.abs[names(invariant.linn.pls.sgcca.abs) == x])})

head(invariant.linn.pls$loadings$X); head(invariant.linn.pls.sgcca$loadings$X);

pred.invariant.linn.pls <- predict(invariant.linn.pls, newdata)
pred.invariant.linn.pls.sgcca <- unlist(predict.sgcca(invariant.linn.pls.sgcca, newdata), recursive = FALSE)

names(pred.invariant.linn.pls)

max(abs(pred.invariant.linn.pls$predict) - abs(pred.invariant.linn.pls.sgcca$predict))
max(abs(pred.invariant.linn.pls$variates) - abs(pred.invariant.linn.pls.sgcca$variates))
max(abs(pred.invariant.linn.pls$B.hat) - abs(pred.invariant.linn.pls.sgcca$B.hat))

### Canonical
canonical.linn.pls <- pls(X, Y, mode = "canonical", tol = 1e-25)
canonical.linn.pls.sgcca <- pls.sgcca(X, Y, mode = "canonical", tol = 1e-25)

canonical.linn.pls.abs = abslist(canonical.linn.pls, c("variates", "loadings", "mat.c", "mat.e"))
canonical.linn.pls.sgcca.abs = abslist(canonical.linn.pls.sgcca, c("variates", "loadings", "mat.c", "mat.e"))

lapply(checklist, function(x){all.equal(canonical.linn.pls.abs[names(canonical.linn.pls.abs) == x], canonical.linn.pls.sgcca.abs[names(canonical.linn.pls.sgcca.abs) == x])})

head(canonical.linn.pls$loadings$X); head(canonical.linn.pls.sgcca$loadings$X);

pred.canonical.linn.pls <- predict(canonical.linn.pls, newdata)
pred.canonical.linn.pls.sgcca <- unlist(predict.sgcca(canonical.linn.pls.sgcca, newdata), recursive = FALSE)

names(pred.canonical.linn.pls)

max(abs(pred.canonical.linn.pls$predict) - abs(pred.canonical.linn.pls.sgcca$predict))
max(abs(pred.canonical.linn.pls$variates) - abs(pred.canonical.linn.pls.sgcca$variates))
max(abs(pred.canonical.linn.pls$B.hat) - abs(pred.canonical.linn.pls.sgcca$B.hat))

### Regression
regression.linn.pls <- pls(X, Y, mode = "regression", tol = 1e-25, ncomp = 3)
regression.linn.pls.sgcca <- pls.sgcca(X, Y, mode = "regression", tol = 1e-25, ncomp = 3)

regression.linn.pls.abs = abslist(regression.linn.pls, c("variates", "loadings", "mat.c", "mat.d"))
regression.linn.pls.sgcca.abs = abslist(regression.linn.pls.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(regression.linn.pls.abs[names(regression.linn.pls.abs) == x],
                                        regression.linn.pls.sgcca.abs[names(regression.linn.pls.sgcca.abs) == x])})

pred.regression.linn.pls <- predict(regression.linn.pls, newdata)
pred.regression.linn.pls.sgcca <- unlist(predict.sgcca(regression.linn.pls.sgcca, newdata), recursive = FALSE)

names(pred.regression.linn.pls)

max(abs(pred.regression.linn.pls$predict) - abs(pred.regression.linn.pls.sgcca$predict))
max(abs(pred.regression.linn.pls$variates) - abs(pred.regression.linn.pls.sgcca$variates))
max(abs(pred.regression.linn.pls$B.hat) - abs(pred.regression.linn.pls.sgcca$B.hat))

system.time(perf(regression.linn.pls, method.predict = "all", validation = "loo", progressBar = FALSE))
system.time(perf.sgcca(regression.linn.pls.sgcca, validation = "loo"))

### Example 2 ###
rm(list=ls())
source("pls.sgcca.R")

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
row.names(Y) <- paste("X", row.names(Y), sep = "")

### Classic
classic.tox.pls <- pls(X, Y, ncomp = 3, mode = "classic", tol = 1e-25)
classic.tox.pls.sgcca <- pls.sgcca(X, Y, ncomp = 3, mode = "classic", tol = 1e-25)

checklist = names(classic.tox.pls)
checklist = checklist[!checklist %in% c("iter")]

classic.tox.pls = abslist(classic.tox.pls, c("variates", "loadings", "mat.c"))
classic.tox.pls.sgcca = abslist(classic.tox.pls.sgcca, c("variates", "loadings", "mat.c"))

lapply(checklist, function(x){all.equal(classic.tox.pls[names(classic.tox.pls) == x], classic.tox.pls.sgcca[names(classic.tox.pls.sgcca) == x])})

### Invariant
invariant.tox.pls <- pls(X, Y, ncomp = 3, mode = "invariant", tol = 1e-25)
invariant.tox.pls.sgcca <- pls.sgcca(X, Y, ncomp = 3, mode = "invariant", tol = 1e-25)

invariant.tox.pls = abslist(invariant.tox.pls, c("variates", "loadings", "mat.c"))
invariant.tox.pls.sgcca = abslist(invariant.tox.pls.sgcca, c("variates", "loadings", "mat.c"))

lapply(checklist, function(x){all.equal(invariant.tox.pls[names(invariant.tox.pls) == x], invariant.tox.pls.sgcca[names(invariant.tox.pls.sgcca) == x])})

### Canonical
canonical.tox.pls <- pls(X, Y, ncomp = 3, mode = "canonical", tol = 1e-25)
canonical.tox.pls.sgcca <- pls.sgcca(X, Y, ncomp = 3, mode = "canonical", tol = 1e-25)

canonical.tox.pls = abslist(canonical.tox.pls, c("variates", "loadings", "mat.c", "mat.e"))
canonical.tox.pls.sgcca = abslist(canonical.tox.pls.sgcca, c("variates", "loadings", "mat.c", "mat.e"))

lapply(checklist, function(x){all.equal(canonical.tox.pls[names(canonical.tox.pls) == x], canonical.tox.pls.sgcca[names(canonical.tox.pls.sgcca) == x])})

### Regression
regression.tox.pls <- pls(X, Y, ncomp = 3, mode = "regression", tol = 1e-25)
regression.tox.pls.sgcca <- pls.sgcca(X, Y, ncomp = 3, mode = "regression", tol = 1e-25)

regression.tox.pls = abslist(regression.tox.pls, c("variates", "loadings", "mat.c", "mat.d"))
regression.tox.pls.sgcca = abslist(regression.tox.pls.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(regression.tox.pls[names(regression.tox.pls) == x], regression.tox.pls.sgcca[names(regression.tox.pls.sgcca) == x])})


### Example 3 ### Book M. Tenenhaus "La regression PLS - Theorie et pratique" (page 78)
rm(list=ls())
source("pls.sgcca.R"); source("predict.sgcca.R"); source("predict.sgcca.R");
source("perf.sgcca.R"); source("mix.options.R")

X <- matrix(c(0.00, 0.23, 0.00, 0.00, 0.00, 0.74, 0.03,
              0.00, 0.10, 0.00, 0.00, 0.12, 0.74, 0.04,
              0.00, 0.00, 0.00, 0.10, 0.12, 0.74, 0.04,
              0.00, 0.49, 0.00, 0.00, 0.12, 0.37, 0.02,
              0.00, 0.00, 0.00, 0.62, 0.12, 0.18, 0.08,
              0.00, 0.62, 0.00, 0.00, 0.00, 0.37, 0.01,
              0.17, 0.27, 0.10, 0.38, 0.00, 0.00, 0.08,
              0.17, 0.19, 0.10, 0.38, 0.02, 0.06, 0.08,
              0.17, 0.21, 0.10, 0.38, 0.00, 0.06, 0.08,
              0.17, 0.15, 0.10, 0.38, 0.02, 0.10, 0.08, 
              0.21, 0.36, 0.12, 0.25, 0.00, 0.00, 0.06,
              0.00, 0.00, 0.00, 0.55, 0.00, 0.37, 0.08), nrow = 12, ncol = 7, byrow = TRUE,
            dimnames = list(paste("Ind", 1:12), paste("Col", 1 : 7)))

Y <- c(98.7, 97.8, 96.6, 92.0, 86.6, 91.2, 81.9, 83.1, 82.4, 83.2, 81.4, 88.1)

example3.pls.reg <- spls(X, Y, mode = "regression", ncomp = 4)
example3.pls.reg.perf = perf(example3.pls.reg, validation = "loo")
example3.pls.reg.perf$Q2
example3.pls.reg.perf$Q2.total

example3.pls.sgcca.reg <- pls.sgcca(X, Y, mode = "regression", tol = 1e-25, ncomp = 4)
example3.pred.sgcca.reg <- predict.sgcca(example3.pls.sgcca.reg, X)
example3.pls.sgcca.reg.perf = perf.sgcca(example3.pls.sgcca.reg, validation = "loo")
example3.pls.sgcca.reg.perf$Q2
example3.pls.sgcca.reg.perf$Q2.total

example3.pls.sgcca.reg$loadings
example3.pls.sgcca.reg$variates

################################
############# SPLS #############
################################
rm(list=ls())
source("spls.sgcca.R"); source("helpers.R")
library(mixOmics);

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
row.names(Y) <- paste("X", row.names(Y), sep = "")

### Regression
regression.tox.spls <- spls(X, Y, ncomp = 2, mode = "regression", keepX = c(50, 50), keepY = c(10, 10))
regression.tox.spls.sgcca <- spls.sgcca(X, Y, ncomp = 2, mode = "regression", keepX = c(50, 50), keepY = c(10, 10))

checklist = names(regression.tox.spls)
checklist = checklist[!checklist %in% c("iter")]

regression.tox.spls = abslist(regression.tox.spls, c("variates", "loadings", "mat.c", "mat.d"))
regression.tox.spls.sgcca = abslist(regression.tox.spls.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(regression.tox.spls[names(regression.tox.spls) == x], regression.tox.spls.sgcca[names(regression.tox.spls.sgcca) == x])})

### Canonical
canonical.tox.spls <- spls(X, Y, ncomp = 2, mode = "canonical", keepX = c(50, 50), keepY = c(10, 10))
canonical.tox.spls.sgcca <- spls.sgcca(X, Y, ncomp = 2, mode = "canonical", keepX = c(50, 50), keepY = c(10, 10))

canonical.tox.spls = abslist(canonical.tox.spls, c("variates", "loadings", "mat.c", "mat.e"))
canonical.tox.spls.sgcca = abslist(canonical.tox.spls.sgcca, c("variates", "loadings", "mat.c", "mat.e"))

lapply(checklist, function(x){all.equal(canonical.tox.spls[names(canonical.tox.spls) == x], canonical.tox.spls.sgcca[names(canonical.tox.spls.sgcca) == x])})

##################################
############# PLS-DA #############
##################################
rm(list=ls())
source("splsda.sgcca.R"); source("sgccda.sgcca.R"); source("predict.sgcca.R");
source("perf.sgcca.R"); source("mix.options.R")
library(mixOmics);

### Example 1 ###
data(breast.tumors)
X <- breast.tumors$gene.exp
# X <- t(na.omit(t(X)))  ### discrepancy with missing data
Y <- breast.tumors$sample$treatment

breast.plsda <- plsda(X, Y, ncomp = 1)
breast.plsda.sgcca <- plsda.sgcca(X, Y, ncomp = 1)

checklist = names(breast.plsda)
checklist = checklist[!checklist %in% c("iter")]

breast.plsda.abs = abslist(breast.plsda, c("variates", "loadings", "mat.c", "mat.d"))
breast.plsda.sgcca.abs = abslist(breast.plsda.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(breast.plsda.abs[names(breast.plsda.abs) == x], 
                                        breast.plsda.sgcca.abs[names(breast.plsda.sgcca.abs) == x])})

breast.plsda.predict = predict(breast.plsda, X)
breast.plsda.predict.sgcca <- predict.sgcca(breast.plsda.sgcca, newdata = X, method = "all")

names(breast.plsda.predict); names(breast.plsda.predict.sgcca)
max(abs(breast.plsda.predict$predict) - abs(breast.plsda.predict.sgcca$predict[[1]]))
max(abs(breast.plsda.predict$variates) - abs(breast.plsda.predict.sgcca$variates[[1]]))
max(abs(breast.plsda.predict$B.hat) - abs(breast.plsda.predict.sgcca$B.hat[[1]]))
max(abs(breast.plsda.predict$centroids) - abs(breast.plsda.predict.sgcca$centroids[[1]]))
breast.plsda.predict$method == breast.plsda.predict.sgcca$method
max(abs(breast.plsda.predict$class[[1]]) - abs(unlist(breast.plsda.predict.sgcca$class, recursive = FALSE)[[1]]))

breast.plsda.perf = perf(breast.plsda, validation = "loo")
breast.plsda.perf.sgcca = perf.sgcca(breast.plsda.sgcca, validation = "loo")
max(breast.plsda.perf$error.rate - breast.plsda.perf.sgcca$error.rate.X$X)

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "breast.plsda.predict.sgcca", "map", "plsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "miscrossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

for (i in 1 : 10){
  folds = breast.plsda.perf = breast.plsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); breast.plsda.perf = perf(breast.plsda, validation = "Mfold", folds = folds)
  set.seed(i); breast.plsda.perf.sgcca = perf.sgcca(breast.plsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(breast.plsda.perf$error.rate - breast.plsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(breast.plsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(breast.plsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

### Example 2 ###
rm(list = ls())
source("plsda.sgcca.R"); source("sgccda.sgcca.R"); source("predict.sgcca.R")
source("perf.sgcca.R"); source("mix.options.R")
require(mixOmics)

data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.factor(liver.toxicity$treatment[, 4])

#### Difference because not the same starting point in the algo between function splsda and plsda in mixOmics package (indeed between pls and spls, see above)
#### To play down this discrepancy, set up tol = 1e-25 instead of 1e-06 
liver.plsda <- plsda(X, Y, ncomp = 2, tol = 1e-25)               
liver.plsda.sgcca <- plsda.sgcca(X, Y, ncomp = 2, tol = 1e-25)

checklist = names(liver.plsda)
checklist = checklist[!checklist %in% c("iter")]

liver.plsda.abs = abslist(liver.plsda, c("variates", "loadings", "mat.d", "mat.c"))
liver.plsda.sgcca.abs = abslist(liver.plsda.sgcca, c("variates", "loadings", "mat.d", "mat.c"))

lapply(checklist, function(x){all.equal(liver.plsda.abs[names(liver.plsda.abs) == x], 
                                        liver.plsda.sgcca.abs[names(liver.plsda.sgcca.abs) == x])})

# Comparison prediction
samp <- sample(1:5, nrow(X), replace = TRUE)  
test <- which(samp == 1)   # testing on the first fold
train <- setdiff(1:nrow(X), test)

plsda.train <- plsda(X[train, ], Y[train], ncomp = 2, tol = 1e-25)
test.predict <- predict(plsda.train, X[test, ], method = "max.dist")

plsda.train.sgcca <- plsda.sgcca(X[train, ], Y[train], ncomp = 2, tol = 1e-25)
test.predict.sgcca <- predict.sgcca(plsda.train.sgcca, X[test, ], method = "max.dist")

names(test.predict); names(test.predict.sgcca)
max(abs(test.predict$predict) - abs(test.predict.sgcca$predict[[1]]))
max(abs(test.predict$variates) - abs(test.predict.sgcca$variates[[1]]))
max(abs(test.predict$B.hat) - abs(test.predict.sgcca$B.hat[[1]]))
max(abs(test.predict$centroids) - abs(test.predict.sgcca$centroids[[1]]))
test.predict$method == test.predict.sgcca$method
max(abs(test.predict$class[[1]]) - abs(unlist(test.predict.sgcca$class, recursive = FALSE)[[1]]))

Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
cbind(Y = as.character(Y[test]), Prediction)

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "liver.plsda.sgcca", "map", "plsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "miscrossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

liver.plsda.perf = perf(liver.plsda, validation = "loo")
liver.plsda.perf.sgcca = perf.sgcca(liver.plsda.sgcca, validation = "loo")
max(liver.plsda.perf$error.rate - liver.plsda.perf.sgcca$error.rate.X$X)

for (i in 1 : 5){
  folds = liver.plsda.perf = liver.plsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); liver.plsda.perf = perf(liver.plsda, validation = "Mfold", folds = folds)
  set.seed(i); liver.plsda.perf.sgcca = perf.sgcca(liver.plsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(liver.plsda.perf$error.rate - liver.plsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(liver.plsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(liver.plsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

###################################
############# SPLS-DA #############
###################################
rm(list=ls())
source("splsda.sgcca.R"); source("sgccda.sgcca.R"); source("predict.sgcca.R")
source("perf.sgcca.R"); source("mix.options.R")
library(mixOmics);

### Example 1 ###
data(breast.tumors)
X <- breast.tumors$gene.exp
#X <- t(na.omit(t(X)))  ### issue with missing data
Y <- breast.tumors$sample$treatment

breast.splsda <- splsda(X, Y, ncomp = 3)#, keepX = c(25, 25))
breast.splsda.sgcca <- splsda.sgcca(X, Y, ncomp =2,  keepX = c(20, 20))
breast.plsda.sgcca <- plsda.sgcca(X, Y, ncomp =2)

checklist = names(breast.splsda)
checklist = checklist[!checklist %in% c("iter")]

breast.splsda.abs = abslist(breast.splsda, c("variates", "loadings", "mat.c", "mat.d"))
breast.splsda.sgcca.abs = abslist(breast.splsda.sgcca, c("variates", "loadings", "mat.c", "mat.d"))

lapply(checklist, function(x){all.equal(breast.splsda.abs[names(breast.splsda.abs) == x], 
                                        breast.splsda.sgcca.abs[names(breast.splsda.sgcca.abs) == x])})

# Comparison prediction
breast.splsda.predict = predict(breast.splsda, X)
breast.splsda.predict.sgcca <- predict.sgcca(breast.splsda.sgcca, newdata = X, method = "all")

names(breast.splsda.predict); names(breast.splsda.predict.sgcca)
max(abs(breast.splsda.predict$predict) - abs(breast.splsda.predict.sgcca$predict[[1]]))
max(abs(breast.splsda.predict$variates) - abs(breast.splsda.predict.sgcca$variates[[1]]))
max(abs(breast.splsda.predict$B.hat) - abs(breast.splsda.predict.sgcca$B.hat[[1]]))
max(abs(breast.splsda.predict$centroids) - abs(breast.splsda.predict.sgcca$centroids[[1]]))
breast.splsda.predict$method == breast.splsda.predict.sgcca$method
max(abs(breast.splsda.predict$class[[1]]) - abs(unlist(breast.splsda.predict.sgcca$class, recursive = FALSE)[[1]]))

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "breast.splsda.predict.sgcca", "map", "splsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "miscrossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

breast.splsda.perf = perf(breast.splsda, validation = "loo")
breast.splsda.perf.sgcca = perf.sgcca(breast.splsda.sgcca, validation = "loo")
max(breast.splsda.perf$error.rate - breast.splsda.perf.sgcca$error.rate.X$X)

for (i in 1 : 5){
  folds = breast.splsda.perf = breast.splsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); breast.splsda.perf = perf(breast.splsda, validation = "Mfold", folds = folds)
  set.seed(i); breast.splsda.perf.sgcca = perf.sgcca(breast.splsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(breast.splsda.perf$error.rate - breast.splsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(breast.splsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(breast.splsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

### Example 2 ###
rm(list=ls())
source("splsda.sgcca.R"); source("sgccda.sgcca.R"); source("predict.sgcca.R")
source("perf.sgcca.R"); source("mix.options.R")
library(mixOmics);

data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
Y <- as.factor(liver.toxicity$treatment[, 4])

liver.splsda = splsda(X, Y, ncomp = 2, keepX = c(5, 5))
liver.splsda.sgcca = splsda.sgcca(X, Y, ncomp = 2, keepX = c(5, 5))

checklist = names(liver.splsda)
checklist = checklist[!checklist %in% c("iter")]

liver.splsda.abs = abslist(liver.splsda, c("variates", "loadings", "mat.d", "mat.c"))
liver.splsda.sgcca.abs = abslist(liver.splsda.sgcca, c("variates", "loadings", "mat.d", "mat.c"))

lapply(checklist, function(x){all.equal(liver.splsda.abs[names(liver.splsda.abs) == x],
                                        liver.splsda.sgcca.abs[names(liver.splsda.sgcca.abs) == x])})

# Comparison prediction
samp <- sample(1:5, nrow(X), replace = TRUE)  
test <- which(samp == 1)   # testing on the first fold
train <- setdiff(1:nrow(X), test)

splsda.train <- splsda(X[train, ], Y[train], ncomp = 2, keepX = c(30, 30))
test.predict <- predict(splsda.train, X[test, ], method = "max.dist")

splsda.train.sgcca <- splsda.sgcca(X[train, ], Y[train], ncomp = 2, keepX = c(30, 30))
test.predict.sgcca <- predict.sgcca(splsda.train.sgcca, X[test, ], method = "max.dist")

names(test.predict); names(test.predict.sgcca)
max(abs(test.predict$predict) - abs(test.predict.sgcca$predict[[1]]))
max(abs(test.predict$variates) - abs(test.predict.sgcca$variates[[1]]))
max(abs(test.predict$B.hat) - abs(test.predict.sgcca$B.hat[[1]]))
max(abs(test.predict$centroid) - abs(test.predict.sgcca$centroid[[1]]))
test.predict$method == test.predict.sgcca$method
max(abs(test.predict$class[[1]]) - abs(unlist(test.predict.sgcca$class, recursive = FALSE)[[1]]))

Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
cbind(Y = as.character(Y[test]), Prediction)

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "liver.splsda.sgcca", "map", "splsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "miscrossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

liver.splsda.perf = perf(liver.splsda, validation = "loo")
liver.splsda.perf.sgcca = perf.sgcca(liver.splsda.sgcca, validation = "loo")
max(liver.splsda.perf$error.rate - liver.splsda.perf.sgcca$error.rate.X$X)

for (i in 1 : 5){
  folds = liver.splsda.perf = liver.splsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); liver.splsda.perf = perf(liver.splsda, validation = "Mfold", folds = folds)
  set.seed(i); liver.splsda.perf.sgcca = perf.sgcca(liver.splsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(liver.splsda.perf$error.rate - liver.splsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(liver.splsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(liver.splsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

### Test perf.splsda ###
rm(list=ls())
source("splsda.sgcca.R"); source("perf.sgcca.R"); source("predict.sgcca.R")
source("perf.sgcca.R"); source("mix.options.R")
library(mixOmics);

data(srbct)
X <- srbct$gene
Y <- srbct$class  

ncomp = 5

srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))  
srbct.splsda.sgcca <- splsda.sgcca(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))  

# with Mfold
# ---------
set.seed(2)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8, method.predict = "all")
error$error.rate

set.seed(2)
error.sgcca <- perf.sgcca(srbct.splsda.sgcca, validation = "Mfold", folds = 8, method.predict = "all")
error.sgcca$error.rate.X

max(abs(error$error.rate - error.sgcca$error.rate.X$X))

cluster = makeCluster(4, type = "SOCK")
mix.options(parallel = cluster)
clusterExport(cluster, c("perf.sgcca", "srbct.splsda.sgcca", "map", "splsda.sgcca", "unmap",
                         "sgcca", "scale2", "cov2", "sgccak", "norm2", "miscrossprod", "BinarySearch",
                         "predict.sgcca", "defl.select", "deflation", "soft.threshold", "soft"))

srbct.splsda.perf = perf(srbct.splsda, validation = "loo")
srbct.splsda.perf.sgcca = perf.sgcca(srbct.splsda.sgcca, validation = "loo")
max(srbct.splsda.perf$error.rate - srbct.splsda.perf.sgcca$error.rate.X$X)

all.equal(srbct.splsda.perf$features$stable, srbct.splsda.perf.sgcca$features$stable$X)
all.equal(srbct.splsda.perf$features$final, srbct.splsda.perf.sgcca$features$final$X)

for (i in 1 : 5){
  folds = srbct.splsda.perf = srbct.splsda.perf.sgcca = NULL
  folds = round(runif(n = 1, min = 1, max = nrow(X)/2))
  set.seed(i); srbct.splsda.perf = perf(srbct.splsda, validation = "Mfold", folds = folds)
  set.seed(i); srbct.splsda.perf.sgcca = perf.sgcca(srbct.splsda.sgcca, validation = "Mfold", folds = folds)
  print(max(abs(srbct.splsda.perf$error.rate - srbct.splsda.perf.sgcca$error.rate.X$X)))
  print(system.time(perf(srbct.splsda, validation = "Mfold", folds = folds)))
  print(system.time(perf.sgcca(srbct.splsda.sgcca, validation = "Mfold", folds = folds)))
}
stopCluster(mix.options()$parallel)

#########################################
############# wrapper.rgcca #############
#########################################
rm(list=ls())
source("/Users/bgautier/Desktop/rGCCA-DA/rGCCA codes/sgccdaFunctions_BG3.R")
library(mixOmics); #library(RGCCA)

data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid, Y)
# with this design, gene expression and lipids are connected to the diet factor
design = matrix(c(0,0,1,
                  0,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

### centroid ###
centroid.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0.5), ncomp = c(2, 2, 1), scheme = "centroid", verbose = FALSE, tol = 1e-25)
centroid.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0.5), verbose = FALSE, scheme = "centroid", tol = 1e-25) 

checklist = names(centroid.wrap.result.rgcca)
checklist = checklist[!checklist %in% c("class", "data", "names")]

lapply(checklist, function(x){all.equal(centroid.wrap.result.rgcca[names(centroid.wrap.result.rgcca) == x], centroid.sgcca.nutrimouse.rgccak[names(centroid.sgcca.nutrimouse.rgccak) == x])})

### horst ###
horst.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0.5), ncomp = c(2, 2, 1), scheme = "horst", verbose = FALSE, tol = 1e-25)
horst.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0.5), verbose = FALSE, scheme = "horst", tol = 1e-25) 

lapply(checklist, function(x){all.equal(horst.wrap.result.rgcca[names(horst.wrap.result.rgcca) == x], horst.sgcca.nutrimouse.rgccak[names(horst.sgcca.nutrimouse.rgccak) == x])})

### factorial ###
factorial.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0.5), ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE, tol = 1e-25)
factorial.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0.5), verbose = FALSE, scheme = "factorial", tol = 1e-25) 

lapply(checklist, function(x){all.equal(factorial.wrap.result.rgcca[names(factorial.wrap.result.rgcca) == x], factorial.sgcca.nutrimouse.rgccak[names(factorial.sgcca.nutrimouse.rgccak) == x])})

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#note: the tau parameter is the regularization parameter

### centroid ###
centroid.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0), ncomp = c(2, 2, 1), scheme = "centroid", verbose = FALSE, tol = 1e-25)
centroid.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0), verbose = FALSE, scheme = "centroid", tol = 1e-25) 

lapply(checklist, function(x){all.equal(centroid.wrap.result.rgcca[names(centroid.wrap.result.rgcca) == x], centroid.sgcca.nutrimouse.rgccak[names(centroid.sgcca.nutrimouse.rgccak) == x])})

### horst ###
horst.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0), ncomp = c(2, 2, 1), scheme = "horst", verbose = FALSE, tol = 1e-25)
horst.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0), verbose = FALSE, scheme = "horst", tol = 1e-25) 

lapply(checklist, function(x){all.equal(horst.wrap.result.rgcca[names(horst.wrap.result.rgcca) == x], horst.sgcca.nutrimouse.rgccak[names(horst.sgcca.nutrimouse.rgccak) == x])})

### factorial ###
factorial.wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0), ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE, tol = 1e-25)
factorial.sgcca.nutrimouse.rgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(1, 1, 0), verbose = FALSE, scheme = "factorial", tol = 1e-25) 

lapply(checklist, function(x){all.equal(factorial.wrap.result.rgcca[names(factorial.wrap.result.rgcca) == x], factorial.sgcca.nutrimouse.rgccak[names(factorial.sgcca.nutrimouse.rgccak) == x])})

#########################################
############# wrapper.sgcca #############
#########################################
rm(list=ls())
source("/Users/bgautier/Desktop/rGCCA-DA/rGCCA codes/sgccdaFunctions_BG.R")
library(mixOmics); #library(RGCCA)

data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
# with this design, gene expression and lipids are connected to the diet factor
design = matrix(c(0,0,1,
                  0,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

### centroid ###
centroid.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "centroid", verbose = FALSE, tol = 1e-25)
centroid.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(.3, .5, 1), verbose = FALSE, scheme = "centroid", tol = 1e-25, sparse = TRUE) 

checklist = names(centroid.wrap.result.sgcca)
checklist = checklist[!checklist %in% c("class", "data", "names")]

lapply(checklist, function(x){all.equal(centroid.wrap.result.sgcca[names(centroid.wrap.result.sgcca) == x], centroid.sgcca.nutrimouse.sgccak[names(centroid.sgcca.nutrimouse.sgccak) == x])})

### horst ###
horst.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "horst", verbose = FALSE, tol = 1e-25)
horst.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(.3, .5, 1), verbose = FALSE, scheme = "horst", tol = 1e-25, sparse = TRUE) 

lapply(checklist, function(x){all.equal(horst.wrap.result.sgcca[names(horst.wrap.result.sgcca) == x], horst.sgcca.nutrimouse.sgccak[names(horst.sgcca.nutrimouse.sgccak) == x])})

### factorial ###
factorial.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE, tol = 1e-25)
factorial.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(.3, .5, 1), verbose = FALSE, scheme = "factorial", tol = 1e-25, sparse = TRUE) 

lapply(checklist, function(x){all.equal(factorial.wrap.result.sgcca[names(factorial.wrap.result.sgcca) == x], factorial.sgcca.nutrimouse.sgccak[names(factorial.sgcca.nutrimouse.sgccak) == x])})

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#note: the penalty parameters will need to be tuned

### centroid ###
centroid.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "centroid", verbose = FALSE, tol = 1e-25)
centroid.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(0.3, 0.5, 1), verbose = FALSE, scheme = "centroid", tol = 1e-25, sparse = TRUE) 

lapply(checklist, function(x){all.equal(centroid.wrap.result.sgcca[names(centroid.wrap.result.sgcca) == x], centroid.sgcca.nutrimouse.sgccak[names(centroid.sgcca.nutrimouse.sgccak) == x])})

### horst ###
horst.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "horst", verbose = FALSE, tol = 1e-25)
horst.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(0.3, 0.5, 1), verbose = FALSE, scheme = "horst", tol = 1e-25, sparse = TRUE) 

lapply(checklist, function(x){all.equal(horst.wrap.result.sgcca[names(horst.wrap.result.sgcca) == x], horst.sgcca.nutrimouse.sgccak[names(horst.sgcca.nutrimouse.sgccak) == x])})

### factorial ###
factorial.wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE, tol = 1e-25)
factorial.sgcca.nutrimouse.sgccak <- sgcca(A = data, ncomp = c(2, 2, 1), C = design, c1 = c(0.3, 0.5, 1), verbose = FALSE, scheme = "factorial", tol = 1e-25, sparse = TRUE, max.iter = 1000) 

lapply(checklist, function(x){all.equal(factorial.wrap.result.sgcca[names(factorial.wrap.result.sgcca) == x], factorial.sgcca.nutrimouse.sgccak[names(factorial.sgcca.nutrimouse.sgccak) == x])})
