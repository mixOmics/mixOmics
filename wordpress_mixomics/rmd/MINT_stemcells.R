## ----global_options, include=FALSE---------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, #dev = 'jpeg',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 8, fig.width=9)

## ----message = TRUE------------------------------------------------------
library(mixOmics)

## ------------------------------------------------------------------------
data(stemcells)

#the combined data set X
X = stemcells$gene
dim(X) 

# the outcome vector Y:  
Y = stemcells$celltype 
length(Y) 
summary(Y)

# the vector indicating each independent study
study = stemcells$study
# number of samples per study:
summary(study)

# experimental design
table(Y,study)

## ------------------------------------------------------------------------
mint.plsda.res.perf = mint.plsda(X = X, Y = Y, study = study, ncomp = 5)

set.seed(2543)  # for reproducible result in this example
perf.mint.plsda.cell <- perf(mint.plsda.res.perf, validation = "Mfold", folds = 5, 
                  progressBar = FALSE, auc = TRUE) 

## ------------------------------------------------------------------------
plot(perf.mint.plsda.cell, col = color.mixo(5:7))

## ------------------------------------------------------------------------
perf.mint.plsda.cell$global.error

## ------------------------------------------------------------------------
perf.mint.plsda.cell$choice.ncomp

## ------------------------------------------------------------------------
mint.plsda.res = mint.plsda(X = X, Y = Y, study = study, ncomp = 2)
#mint.plsda.res # lists the different functions
plotIndiv(mint.plsda.res, legend = TRUE, title = 'MINT PLS-DA', 
          subtitle = 'stem cell study', ellipse = T)

## ---- eval = TRUE, include = TRUE----------------------------------------
tune.mint = tune(X = X, Y = Y, study = study, ncomp = 2, test.keepX = seq(1, 100, 1), 
method = 'mint.splsda', dist = "max.dist", progressBar = FALSE)

# tune.mint   # lists the different types of outputs

# mean error rate per component and per tested keepX value
# tune.mint$error.rate

## ------------------------------------------------------------------------
# optimal number of components
tune.mint$choice.ncomp #tune.mint$choice.ncomp # tell us again than ncomp=1 is sufficient

# optimal keepX
tune.mint$choice.keepX

plot(tune.mint, col = color.jet(2))

## ------------------------------------------------------------------------
mint.splsda.res = mint.splsda(X = X, Y = Y, study = study, ncomp = 2,  
                              keepX = tune.mint$choice.keepX)

#mint.splsda.res   # lists useful functions that can be used with a MINT object

## ------------------------------------------------------------------------
selectVar(mint.splsda.res, comp = 1)

## ------------------------------------------------------------------------
plotIndiv(mint.splsda.res, study = 'global', legend = TRUE, title = 'MINT sPLS-DA', 
          subtitle = 'Global', ellipse=T)

## ------------------------------------------------------------------------
plotIndiv(mint.splsda.res, study = 'all.partial',  title = 'MINT sPLS-DA', 
          subtitle = paste("Study",1:4))

## ------------------------------------------------------------------------
plotArrow(mint.splsda.res)

## ------------------------------------------------------------------------
cim(mint.splsda.res, comp = 1, margins=c(10,5), 
    row.sideColors = color.mixo(as.numeric(Y)), row.names = FALSE,
    title = "MINT sPLS-DA, component 1")

## ------------------------------------------------------------------------
plotLoadings(mint.splsda.res, contrib="max", method = 'mean', comp=1, 
             study="all.partial", legend=FALSE, title="Contribution on comp 1", 
             subtitle = paste("Study",1:4))

## ------------------------------------------------------------------------
set.seed(123)  # for reproducibility of the results
perf.mint = perf(mint.splsda.res, progressBar = FALSE, dist = 'max.dist')

perf.mint$global.error

## ------------------------------------------------------------------------
plot(perf.mint, col = color.mixo(5))

## ------------------------------------------------------------------------
# we predict on study 3
ind.test = which(study == "3")
test.predict <- predict(mint.splsda.res, newdata = X[ind.test, ], dist = "max.dist",
                        study.test = factor(study[ind.test]))
Prediction <- test.predict$class$max.dist[, 2]

# the confusion table compares the real subtypes with the predicted subtypes
get.confusion_matrix(truth = Y[ind.test],
                     predicted = Prediction)

## ------------------------------------------------------------------------
auc.mint.splsda = auroc(mint.splsda.res, roc.comp = 2)

## ------------------------------------------------------------------------
auc.mint.splsda = auroc(mint.splsda.res, roc.comp = 2, roc.study = '2')

## ------------------------------------------------------------------------
sessionInfo()

## ---- include = FALSE----------------------------------------------------
# extract R code
#purl("MINT_stemcells.Rmd")

