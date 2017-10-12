## ----global_options, include=FALSE---------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center', 
                      fig.show=TRUE, fig.keep = 'all', out.width = '50%') 

## ----message = FALSE-----------------------------------------------------
library(mixOmics)

## ------------------------------------------------------------------------
data(srbct)
X = srbct$gene  #the gene expression data
dim(X)

summary(srbct$class)

## ------------------------------------------------------------------------
pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE)
#pca.srbct #outputs the explained variance per component
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)

## ------------------------------------------------------------------------
plotIndiv(pca.srbct, group = srbct$class, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

## ------------------------------------------------------------------------
Y = srbct$class 
summary(Y)        #outcome categories

## ------------------------------------------------------------------------
srbct.plsda <- plsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = srbct$class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA on SRBCT')

## ------------------------------------------------------------------------
# with background
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
#optional: xlim = c(-40,40), ylim = c(-30,30))

plotIndiv(srbct.plsda, comp = 1:2,
          group = srbct$class, ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)

## ------------------------------------------------------------------------
# takes a couple of minutes to run
set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                  progressBar = FALSE, auc = TRUE, nrepeat = 10) 

## ------------------------------------------------------------------------
# perf.plsda.srbct$error.rate  # error rates
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

## ------------------------------------------------------------------------
auc.plsda = auroc(srbct.plsda, roc.comp = 6)

## ----include = FALSE, eval = FALSE---------------------------------------
## # run internally and saved
## 
## #set.seed(1234) # for reproducibility, only when the `cpus' argument is not used
## # grid of possible keepX values that will be tested for each component
## list.keepX <- c(1:10,  seq(20, 300, 10))
## 
## t1 = proc.time()
## tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 6, validation = 'Mfold', folds = 5,
##                            progressBar = FALSE, dist = 'max.dist', measure = "BER",
##                           test.keepX = list.keepX, nrepeat = 10, cpus = 2)
## t2 = proc.time()
## running_time = t2 - t1; running_time # running time
## 
## error <- tune.splsda.srbct$error.rate
## ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
## select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]
## 
## save(tune.splsda.srbct, ncomp, list.keepX, error, select.keepX, running_time,  file = 'RData/result-SRBCT-sPLSDA.RData')

## ----include = FALSE-----------------------------------------------------
load('RData/result-SRBCT-sPLSDA.Rdata')

## ----eval = FALSE--------------------------------------------------------
## # grid of possible keepX values that will be tested for each component
## list.keepX <- c(1:10,  seq(20, 300, 10))
## 
## tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 6, validation = 'Mfold', folds = 5,
##                            progressBar = TRUE, dist = 'max.dist', measure = "BER",
##                           test.keepX = list.keepX, nrepeat = 10, cpus = 2)
## 

## ------------------------------------------------------------------------
error <- tune.splsda.srbct$error.rate  # error rate per component for the keepX grid

## ------------------------------------------------------------------------
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests
ncomp

## ------------------------------------------------------------------------
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

## ------------------------------------------------------------------------
plot(tune.splsda.srbct, col = color.jet(6))

## ------------------------------------------------------------------------
splsda.srbct <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX) 

## ------------------------------------------------------------------------
plotIndiv(splsda.srbct, comp = c(1,2),
          group = srbct$class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 1 & 2')

## ------------------------------------------------------------------------
plotIndiv(splsda.srbct, comp = c(1,3),
          group = srbct$class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 1 & 3')

## ------------------------------------------------------------------------
auc.splsda = auroc(splsda.srbct, roc.comp = 2)

## ------------------------------------------------------------------------
auc.splsda = auroc(splsda.srbct, roc.comp = ncomp)

## ------------------------------------------------------------------------
set.seed(40) # for reproducibility, only when the `cpus' argument is not used
# takes about 1 min to run
perf.srbct <- perf(splsda.srbct, validation = "Mfold", folds = 5,
                   dist = 'max.dist', nrepeat = 10,
                   progressBar = FALSE) 

## ------------------------------------------------------------------------
# perf.srbct  # lists the different outputs
perf.srbct$error.rate
plot(perf.srbct, col = color.mixo(5))

## ------------------------------------------------------------------------
par(mfrow=c(1,3))
plot(perf.srbct$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)
plot(perf.srbct$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)
plot(perf.srbct$features$stable[[3]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 3', las =2)
par(mfrow=c(1,1))

## ------------------------------------------------------------------------
# here we match the selected variables to the stable features
ind.match = match(selectVar(splsda.srbct, comp = 1)$name, 
                  names(perf.srbct$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.srbct$features$stable[[1]][ind.match])

data.frame(selectVar(splsda.srbct, comp = 1)$value, Freq)

## ------------------------------------------------------------------------
plotLoadings(splsda.srbct, comp = 1, title = 'Loadings on comp 1', 
             contrib = 'max', method = 'mean')

## ------------------------------------------------------------------------
plotLoadings(splsda.srbct, comp = 2, title = 'Loadings on comp 2', 
             contrib = 'max', method = 'mean')

## ------------------------------------------------------------------------
plotLoadings(splsda.srbct, comp = 3, title = 'Loadings on comp 3', 
             contrib = 'max', method = 'mean')

## ------------------------------------------------------------------------
plotLoadings(splsda.srbct, comp = 1, title = 'Loadings on comp 1', 
             contrib = 'min', method = 'mean')

## ------------------------------------------------------------------------
cim(splsda.srbct)

## ------------------------------------------------------------------------
cim(splsda.srbct, comp=1, title ="Component 1")

## ------------------------------------------------------------------------
plotArrow(splsda.srbct, legend=T)

## ------------------------------------------------------------------------
sessionInfo()

## ---- include = FALSE----------------------------------------------------
# extract R code
#purl("PLSDA_SRBCT.Rmd")

