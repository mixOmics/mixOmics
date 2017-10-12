## ----global_options, include=FALSE---------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center', 
                      fig.show=TRUE, fig.keep = 'all', out.width = '50%') 

## ----message = TRUE------------------------------------------------------
library(mixOmics)

## ------------------------------------------------------------------------
data('breast.TCGA')
# extract training data
data = list(mRNA = breast.TCGA$data.train$mrna, 
            miRNA = breast.TCGA$data.train$mirna, 
            proteomics = breast.TCGA$data.train$protein)

# check dimension
lapply(data, dim)

# outcome
Y = breast.TCGA$data.train$subtype
summary(Y)

## ------------------------------------------------------------------------
design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0

design 

## ------------------------------------------------------------------------
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, 
                           design = design)

set.seed(123) # for reproducibility, only when the `cpus' argument is not used
# this code takes a couple of min to run
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)

#perf.diablo  # lists the different outputs
plot(perf.diablo) 

## ------------------------------------------------------------------------
perf.diablo$choice.ncomp$WeightedVote
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

## ---- eval=FALSE, include = FALSE----------------------------------------
## #set.seed(123) # for reproducibility, only when the `cpus' argument is not used
## test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
## 
## t1 = proc.time()
## tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
##                               test.keepX = test.keepX, design = design,
##                               validation = 'Mfold', folds = 10, nrepeat = 1,
##                               cpus = 2, dist = "centroids.dist")
## t2 = proc.time()
## running_time = t2 - t1; running_time
## 
## list.keepX = tune.TCGA$choice.keepX
## list.keepX
## 
## save(tune.TCGA,list.keepX, running_time, file = 'RData/result-TCGA-diablo_design0.1.RData')

## ---- include =TRUE, eval=FALSE------------------------------------------
## #set.seed(123) # for reproducibility, only when the `cpus' argument is not used
## test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
## 
## tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
##                               test.keepX = test.keepX, design = design,
##                               validation = 'Mfold', folds = 10, nrepeat = 1,
##                               cpus = 2, dist = "centroids.dist")
## 
## 
## list.keepX = tune.TCGA$choice.keepX
## list.keepX

## ---- include =FALSE, eval=TRUE------------------------------------------
load('RData/result-TCGA-diablo_design0.1.RData')
list.keepX = tune.TCGA$choice.keepX
tune.TCGA$choice.keepX

## ------------------------------------------------------------------------
list.keepX = list(mRNA = c(6,14), miRNA = c(5,18), proteomics = c(6,7)) # from tuning step

## ------------------------------------------------------------------------
sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)
#sgccda.res   # list the different functions of interest related to that object

## ------------------------------------------------------------------------
sgccda.res$design

## ------------------------------------------------------------------------
# mrna variables selected on component 1
selectVar(sgccda.res, block = 'mRNA', comp = 1)$mRNA$name 

## ------------------------------------------------------------------------
plotDiablo(sgccda.res, ncomp = 1)

## ------------------------------------------------------------------------
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

## ------------------------------------------------------------------------
plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

## ------------------------------------------------------------------------
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(2,2,2), col = c('darkorchid', 'brown1', 'lightgreen'))

## ------------------------------------------------------------------------
circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

## ---- eval = TRUE--------------------------------------------------------
network(sgccda.res, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)

## ----eval = FALSE--------------------------------------------------------
## # not run
## library(igraph)
## my.network = network(sgccda.res, blocks = c(1,2,3),
##         color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)
## write.graph(my.network$gR, file = "myNetwork.gml", format = "gml")

## ------------------------------------------------------------------------
plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')

## ---- eval = TRUE--------------------------------------------------------
cimDiablo(sgccda.res)

## ------------------------------------------------------------------------
set.seed(123)# for reproducibility, only when the `cpus' argument is not used
perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, 
                   dist = 'centroids.dist')
#perf.diablo  # lists the different outputs

# Performance with Majority vote
perf.diablo$MajorityVote.error.rate

# Performance with Weighted prediction
perf.diablo$WeightedVote.error.rate

## ------------------------------------------------------------------------
auc.splsda = auroc(sgccda.res, roc.block = "miRNA", roc.comp = 2)

## ------------------------------------------------------------------------
# prepare test set data: here one block (proteins) is missing
data.test.TCGA = list(mRNA = breast.TCGA$data.test$mrna, 
                      miRNA = breast.TCGA$data.test$mirna)

predict.diablo = predict(sgccda.res, newdata = data.test.TCGA)
# the warning message will inform us that one block is missing
#predict.diablo # list the different outputs

## ------------------------------------------------------------------------
confusion.mat = get.confusion_matrix(truth = breast.TCGA$data.test$subtype, 
                     predicted = predict.diablo$WeightedVote$centroids.dist[,2])
confusion.mat
get.BER(confusion.mat)

## ------------------------------------------------------------------------
sessionInfo()

## ---- include = FALSE----------------------------------------------------
# extract R code
#purl("DIABLO_TCGA.Rmd")

