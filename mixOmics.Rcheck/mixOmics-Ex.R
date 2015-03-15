pkgname <- "mixOmics"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "mixOmics-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('mixOmics')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("cim")
### * cim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cim
### Title: Clustered Image Maps (CIMs) ("heat maps")
### Aliases: cim cim.default cim.rcc cim.spls cim.pls
### Keywords: multivariate iplot hplot graphs cluster

### ** Examples

## default method
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

cim(cor(X, Y), dendrogram = "none")

## CIM representation for objects of class 'rcc'
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

dends <- cim(nutri.res, comp = 1:3, xlab = "genes", 
             ylab = "lipids", margins = c(5, 6))

op <- par(mar = c(5, 4, 4, 4), cex = 0.8)			 
plot(dends$ddr, axes = FALSE, horiz = TRUE)
par(op)

## interactive 'zoom' 
## Not run: 
##D cim(nutri.res, comp = 1:3, zoom = TRUE)
##D ## select the region and "see" the zoom-out region
## End(Not run)

## CIM representation for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

cim(toxicity.spls, comp = 1:3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("image.tune.rcc")
### * image.tune.rcc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: image
### Title: Plot the cross-validation score.
### Aliases: image.tune.rcc
### Keywords: dplot hplot

### ** Examples

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

## this can take some seconds
cv.score <- tune.rcc(X, Y, validation = "Mfold", plt = FALSE)
image(cv.score)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("image.tune.rcc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("imgCor")
### * imgCor

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: imgCor
### Title: Image Maps of Correlation Matrices between two Data Sets
### Aliases: imgCor
### Keywords: multivariate dplot

### ** Examples

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

## 'combine' type plot (default)
imgCor(X, Y)

## 'separate' type plot
## Not run: 
##D imgCor(X, Y, type = "separate")
##D 
##D ## 'separate' type plot without the name of datas
##D imgCor(X, Y, X.names = FALSE, Y.names = FALSE, type = "separate")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("imgCor", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ipca")
### * ipca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ipca
### Title: Independent Principal Component Analysis
### Aliases: ipca
### Keywords: algebra

### ** Examples

data(liver.toxicity)

# implement IPCA on a microarray dataset
ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode="deflation")
ipca.res

# samples representation
plotIndiv(ipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 0.5, 
          col = as.numeric(as.factor(liver.toxicity$treatment[, 4])))
## Not run: 
##D plot3dIndiv(ipca.res, cex = 0.01,
##D             col = as.numeric(as.factor(liver.toxicity$treatment[, 4])))
## End(Not run)

# variables representation
plotVar(ipca.res, var.label = TRUE, cex = 0.5)

## Not run: 
##D plot3dVar(ipca.res, rad.in = 0.5, cex = 0.5, 
##D           col = as.numeric(as.factor(liver.toxicity$treatment[, 4])))
##D           
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ipca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jet.colors")
### * jet.colors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jet.colors
### Title: Jet Colors Palette
### Aliases: jet.colors
### Keywords: color

### ** Examples

## Jet Color Scales Strips
def.par <- par(no.readonly = TRUE)
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
    image(matrix(z, ncol = 1), col = jet.colors(n), 
          xaxt = "n", yaxt = "n", main = paste("n = ", n))
    box()
    par(usr = c(-1, 1, -1, 1))		
    axis(1, at = c(-1, 0, 1))
}
par(def.par)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jet.colors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("mat.rank")
### * mat.rank

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mat.rank
### Title: Matrix Rank
### Aliases: mat.rank
### Keywords: algebra

### ** Examples

## Hilbert matrix
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
mat <- hilbert(16)
mat.rank(mat)

## Hilbert matrix with missing data
idx.na <- matrix(sample(c(0, 1, 1, 1, 1), 36, replace = TRUE), ncol = 6)
m.na <- m <- hilbert(9)[, 1:6]
m.na[idx.na == 0] <- NA
mat.rank(m)
mat.rank(m.na)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mat.rank", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("multilevel")
### * multilevel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: multilevel
### Title: Multilevel analysis for repeated measurements (cross-over
###   design)
### Aliases: multilevel multilevel.spls multilevel.splsda
### Keywords: regression multivariate

### ** Examples

## First example: one-factor analysis with sPLS-DA, selecting a subset of variables
# as in the paper Liquet et al.
#--------------------------------------------------------------
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)

# multilevel sPLS-DA model
res.1level <- multilevel(X, ncomp = 3, design = design,
                         method = "splsda", keepX = c(30, 137, 123))

# set up colors for plotIndiv
col.stim <- c("darkblue", "purple", "green4","red3")
col.stim <- col.stim[as.numeric(Y)]
plotIndiv(res.1level, ind.names = Y, col = col.stim)


## Second example: two-factor analysis with sPLS-DA, selecting a subset of variables
# as in the paper Liquet et al.
#--------------------------------------------------------------
## Not run: 
##D data(vac18.simulated) # simulated data
##D 
##D X <- vac18.simulated$genes
##D design <- data.frame(sample = vac18.simulated$sample,
##D                      stimu = vac18.simulated$stimulation,
##D                      time = vac18.simulated$time)
##D 
##D res.2level <- multilevel(X, ncomp = 2, design = design,
##D                          keepX = c(200, 200), method = 'splsda')
##D 
##D # set up colors and pch for plotIndiv
##D col.stimu <- as.numeric(design$stimu)
##D pch.time <- c(20, 4)[as.numeric(design$time)]
##D 
##D plotIndiv(res.2level, col = col.stimu, ind.names = FALSE,
##D           pch = pch.time)
##D legend('bottomright', legend = levels(design$stimu),
##D        col = unique(col.stimu), pch = 20, cex = 0.8, 
##D        title = "Stimulation")
##D legend('topright', col = 'black', legend = levels(design$time),  
##D        pch = unique(pch.time), cex = 0.8, title = "Time")
## End(Not run)       

## Third example: one-factor analysis with sPLS, selecting a subset of variables
#--------------------------------------------------------------
## Not run: 
##D data(liver.toxicity)
##D # note: we made up those data, pretending they are repeated measurements
##D repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
##D                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
##D                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
##D                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
##D summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
##D 
##D # this is a spls (unsupervised analysis) so no need to mention any factor in design
##D # we only perform a one level variation split
##D design <- data.frame(sample = repeat.indiv) 
##D res.spls.1level <- multilevel(X = liver.toxicity$gene,
##D                                        Y=liver.toxicity$clinic,
##D                                        design = design,
##D                                        ncomp = 3,
##D                                        keepX = c(50, 50, 50), keepY = c(5, 5, 5),
##D                                        method = 'spls', mode = 'canonical')
##D 
##D # set up colors and pch for plotIndiv
##D col.stimu <- as.numeric(as.factor(design$stimu))
##D 
##D plotIndiv(res.spls.1level, rep.space = 'X-variate', ind.names = FALSE, 
##D           col = col.stimu, pch = 20)
##D title(main = 'Gene expression space')
##D plotIndiv(res.spls.1level, rep.space = 'Y-variate', ind.names = FALSE,
##D           col = col.stimu, pch = 20)
##D title(main = 'Clnical measurements space')
##D legend('bottomright', legend = levels(as.factor(design$stimu)),
##D        col = unique(col.stimu), pch = 20, cex = 0.8, 
##D        title = "Dose")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("multilevel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("network")
### * network

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: network
### Title: Relevance Network for (r)CCA and (s)PLS regression
### Aliases: network.default network.rcc network.pls network.spls network
### Keywords: multivariate graphs dplot hplot iplot

### ** Examples


## network representation for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

## Not run: 
##D    # may not work on the Linux version, use Windows instead
##D network(nutri.res, comp = 1:3, threshold = 0.6)
## End(Not run)

## Changing the attributes of the network
## Not run: 
##D network(nutri.res, comp = 1:3, threshold = 0.45,
##D         color.node = c("mistyrose", "lightcyan"),,
##D         shape.node = c("circle", "rectangle"), 
##D         color.edge = jet.colors(8),
##D         lty.edge = c("solid", "solid"), lwd.edge = c(2, 2), 
##D         show.edge.labels = FALSE)
## End(Not run)

## interactive 'threshold' 
## Not run: 
##D network(nutri.res, comp = 1:3, threshold = 0.55, interactive = TRUE)
##D ## select the 'threshold' and "see" the new network
## End(Not run)

## network representation for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
## Not run: 
##D network(toxicity.spls, comp = 1:3, threshold = 0.8, 
##D         X.names = NULL, Y.names = NULL, keep.var = TRUE,
##D         color.node = c("mistyrose", "lightcyan"),
##D         shape.node = c("rectangle", "circle"),
##D         color.edge = c("red", "blue"),
##D         lty.edge = c("solid", "solid"), lwd.edge = c(1, 1), 
##D         show.edge.labels = FALSE, interactive = FALSE)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("network", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("nipals")
### * nipals

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: nipals
### Title: Non-linear Iterative Partial Least Squares (NIPALS) algorithm
### Aliases: nipals
### Keywords: algebra multivariate

### ** Examples

## Hilbert matrix
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
X.na <- X <- hilbert(9)[, 1:6]

## Hilbert matrix with missing data
idx.na <- matrix(sample(c(0, 1, 1, 1, 1), 36, replace = TRUE), ncol = 6)
X.na[idx.na == 0] <- NA
X.rec <- nipals(X.na, reconst = TRUE)$rec
round(X, 2)
round(X.rec, 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("nipals", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pca")
### * pca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pca
### Title: Principal Components Analysis
### Aliases: pca
### Keywords: algebra

### ** Examples

data(multidrug)

## this data set contains missing values, therefore 
## the 'prcomp' function cannot be applied
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
plot(pca.res)
print(pca.res)
biplot(pca.res, xlabs = multidrug$cell.line$Class, cex = 0.7)

# samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 0.5, 
          col = as.numeric(as.factor(multidrug$cell.line$Class)))
## Not run: 
##D plot3dIndiv(pca.res, cex = 0.2,
##D             col = as.numeric(as.factor(multidrug$cell.line$Class)))
## End(Not run)
# variables representation
plotVar(pca.res, var.label = TRUE)
## Not run: 
##D plot3dVar(pca.res, rad.in = 0.5, var.label = TRUE, cex = 0.5)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("perf")
### * perf

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: perf
### Title: Compute evaluation criteria for PLS, sPLS, PLS-DA and sPLS-DA
### Aliases: perf perf.pls perf.spls perf.plsda perf.splsda
### Keywords: regression multivariate

### ** Examples

## validation for objects of class 'pls' (regression)
# ----------------------------------------
## Not run: 
##D data(liver.toxicity)
##D X <- liver.toxicity$gene
##D Y <- liver.toxicity$clinic
##D 
##D 
##D # try tune the number of component to choose
##D # ---------------------
##D # first learn the full model
##D liver.pls <- pls(X, Y, ncomp = 10)
##D 
##D # with 5-fold cross validation: we use the same parameters as in model above
##D # but we perform cross validation to compute the MSEP, Q2 and R2 criteria
##D # ---------------------------
##D liver.val <- perf(liver.pls, validation = "Mfold", folds = 5)
##D 
##D # Q2 total should decrease until it reaches a threshold
##D liver.val$Q2.total
##D 
##D # ncomp = 3 is enough
##D plot(liver.val$Q2.total, type = 'l', col = 'red', ylim = c(-0.1, 0.5), 
##D 	xlab = 'PLS components', ylab = 'Q2 total')
##D abline(h = 0.0975, col = 'darkgreen')
##D legend('topright', col = c('red', 'darkgreen'), legend = c('Q2 total', 'threshold 0.0975')
##D 	, lty = 1)
##D title('Liver toxicity PLS 5-fold, Q2 values')
##D 
##D #have a look at the other criteria
##D # ----------------------
##D # R2
##D liver.val$R2
##D matplot(t(liver.val$R2), type = 'l', xlab = 'PLS components', ylab = 'R2 for each variable')
##D title('Liver toxicity PLS 5-fold, R2 values')
##D 
##D # MSEP
##D liver.val$MSEP
##D matplot(t(liver.val$MSEP), type = 'l', xlab = 'PLS components', ylab = 'MSEP for each variable')
##D title('Liver toxicity PLS 5-fold, MSEP values')
##D 
##D 
##D ## validation for objects of class 'spls' (regression)
##D # ----------------------------------------
##D ncomp = 7
##D # first, learn the model on the whole data set
##D model.spls = spls(X, Y, ncomp = ncomp, mode = 'regression',
##D 	 keepX = c(rep(5, ncomp)), keepY = c(rep(2,ncomp)))
##D 
##D 
##D # with leave-one-out cross validation
##D set.seed(45)
##D model.spls.loo.val <- perf(model.spls, validation = "loo")
##D 
##D #Q2 total
##D model.spls.loo.val$Q2.total
##D 
##D # R2:we can see how the performance degrades when ncomp increases
##D # results are similar to 5-fold
##D model.spls.loo.val$R2
##D 
##D 
##D ## validation for objects of class 'splsda' (classification)
##D # ----------------------------------------
##D data(srbct)
##D X <- srbct$gene
##D Y <- srbct$class  
##D 
##D ncomp = 5
##D 
##D srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))  
##D 
##D # with Mfold
##D # ---------
##D set.seed(45)
##D error <- perf(srbct.splsda, validation = "Mfold", folds = 8, 
##D                method.predict = "all")
##D 
##D plot(error, type = "l")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("perf", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pheatmap.multilevel")
### * pheatmap.multilevel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pheatmap.multilevel
### Title: Clustered heatmap
### Aliases: pheatmap.multilevel pheatmap.multilevel.splsda1fact
###   pheatmap.multilevel.splsda2fact
### Keywords: regression multivariate

### ** Examples

## First example: one-factor analysis with sPLS-DA
# -------------------
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation

design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)
vac18.splsda.multilevel <- multilevel(X, ncomp = 3, design = design,
                                         method = "splsda", keepX = c(30, 137, 123))


# set up colors for pheatmap
col.samp <- c("lightgreen", "red", "lightblue", "darkorange",
              "purple", "maroon", "blue", "chocolate", "turquoise",
              "tomato1", "pink2", "aquamarine")
col.stimu <- c("darkblue", "purple", "green4","red3")
col.stimu <- col.stimu[as.numeric(Y)]
col.stimu <- unique(col.stimu)

pheatmap.multilevel(vac18.splsda.multilevel, 
                    # colors:
                    col_sample = col.samp, 
                    col_stimulation = col.stimu,
                    #labels:
                    label_annotation = c("Subject", "Stimulus"),
                    # scaling:
                    scale = 'row',
                    # distances and clutering
                    clustering_distance_rows = "euclidean", 
                    clustering_distance_cols = "euclidean", 
                    clustering_method = "complete",
                    #  show col/row names and font
                    show_colnames = FALSE,
                    show_rownames = FALSE, 
                    fontsize = 8, 
                    fontsize_row = 3,
                    fontsize_col = 2,
                    border = FALSE, 
                    width = 10)

## Second example: two-factor analysis with sPLS-DA
# --------------------
## Not run: 
##D data(vac18.simulated) # on the simulated data this time
##D 
##D X <- vac18.simulated$genes
##D design <- data.frame(sample = vac18.simulated$sample,
##D                      stimul = vac18.simulated$stimulation,
##D                      time = vac18.simulated$time)
##D 
##D vac18.splsda2.multilevel <- multilevel(X, ncomp = 2, design = design,
##D                             keepX = c(200, 200), method = 'splsda')
##D 
##D # set up colors for each level of pheatmap 
##D col.sample <- c("lightgreen", "red","lightblue","darkorange","purple","maroon") # 6 samples
##D col.time <- c("pink","lightblue1") # two time points
##D col.stimu <- c('green', 'black', 'red', 'blue') # 4 stimulations
##D # set up labels for the 2 levels in design matrix
##D label.stimu <- unique(design[, 2])
##D label.time <- unique(design$time)
##D 
##D pheatmap.multilevel(vac18.splsda2.multilevel,
##D                                 # colors:
##D                                 col_sample=col.sample, 
##D                                 col_stimulation=col.stimu, 
##D                                 col_time=col.time,
##D                                 #labels for each level
##D                                 label_color_stimulation=label.stimu,
##D                                 label_color_time=label.time, 
##D                                 #clustering method
##D                                 clustering_method="ward",
##D                                 #show col/row names and font size
##D                                 show_colnames = FALSE,
##D                                 show_rownames = TRUE,
##D                                 fontsize_row=2)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pheatmap.multilevel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.perf")
### * plot.perf

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.perf
### Title: Plot for model performance
### Aliases: plot.perf
### Keywords: regression multivariate hplot

### ** Examples

require(lattice)

## validation for objects of class 'pls' or 'spls'
## Not run: 
##D data(liver.toxicity)
##D X <- liver.toxicity$gene
##D Y <- liver.toxicity$clinic
##D 
##D liver.pls <- pls(X, Y, ncomp = 3)
##D liver.perf <- perf(liver.pls, validation = "Mfold")
##D 				   
##D plot(liver.perf, criterion = "R2", type = "l", layout = c(2, 2))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.perf", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.rcc")
### * plot.rcc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.rcc
### Title: Canonical Correlations Plot
### Aliases: plot.rcc
### Keywords: multivariate hplot

### ** Examples

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, lambda1 = 0.064, lambda2 = 0.008)

## 'pointplot' type scree 
plot(nutri.res) #(default)

plot(nutri.res, pch = 19, cex = 1.2, 
     col = c(rep("red", 3), rep("darkblue", 18)))

## 'barplot' type scree
plot(nutri.res, scree.type = "barplot")

plot(nutri.res, scree.type = "barplot", density = 20, col = "black")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.rcc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot3dIndiv")
### * plot3dIndiv

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot3dIndiv
### Title: Plot of Individuals (Experimental Units) in three dimensions
### Aliases: plot3dIndiv plot3dIndiv.rcc plot3dIndiv.pls plot3dIndiv.spls
###   plot3dIndiv.plsda plot3dIndiv.splsda plot3dIndiv.pca plot3dIndiv.spca
### Keywords: multivariate hplot dplot

### ** Examples

require(rgl)

## Not run: 
##D ## plot3d of individuals for objects of class 'rcc' 
##D data(nutrimouse)
##D X <- nutrimouse$lipid
##D Y <- nutrimouse$gene
##D nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
##D 
##D col = nutrimouse$diet
##D font = c(rep(1, 20), rep(3, 20))
##D plot3dIndiv(nutri.res, ind.names = nutrimouse$diet, 
##D                 axes.box = "box", font = font, col = col)
##D 				
##D pch = c(rep("s", 20), rep("t", 20))
##D plot3dIndiv(nutri.res, ind.names = FALSE, axes.box = "both", 
##D                 col = col, cex = 1.5, pch = pch)
##D 
##D ## plot3d of individuals for objects of class 'pls' or 'spls'	   
##D data(liver.toxicity)
##D X <- liver.toxicity$gene
##D Y <- liver.toxicity$clinic
##D toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
##D                       keepY = c(10, 10, 10))
##D 		  	
##D Time.Group = liver.toxicity$treatment[, "Time.Group"]				  
##D col <- rep(c("blue", "red", "darkgreen", "darkviolet"), rep(16, 4))
##D plot3dIndiv(toxicity.spls, ind.names = Time.Group, 
##D                  col = col, cex = 0.8)		  
##D 		  				  
##D col <- rainbow(48)[Time.Group]
##D plot3dIndiv(toxicity.spls, ind.names = FALSE, 
##D                  col = col, cex = 0.3, axes.box = "both")	
##D 
##D ## plot3d of individuals for objects of class 'pca'	   
##D data(multidrug)
##D pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
##D 
##D palette(rainbow(9))
##D col = as.numeric(as.factor(multidrug$cell.line$Class))
##D plot3dIndiv(pca.res, cex = 0.25, col = col)
##D palette("default")				 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot3dIndiv", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot3dVar")
### * plot3dVar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot3dVar
### Title: Plot of Variables in three dimensions
### Aliases: plot3dVar plot3dVar.rcc plot3dVar.pls plot3dVar.spls
###   plot3dVar.plsda plot3dVar.splsda plot3dVar.pca plot3dVar.spca
### Keywords: multivariate hplot dplot

### ** Examples

require(rgl)

## Not run: 
##D ## 3D variable representation for objects of class 'rcc'
##D data(nutrimouse)
##D X <- nutrimouse$lipid
##D Y <- nutrimouse$gene
##D nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
##D 
##D \dontrun{
##D # default
##D plot3dVar(nutri.res)
##D 
##D # cutoff active, labeling the variables
##D plot3dVar(nutri.res, cutoff = 0.7, X.label = TRUE, cex = c(0.8, 0.8))
##D }
##D ## 3D variable representation for objects of class 'pls' or 'spls'
##D data(liver.toxicity)
##D X <- liver.toxicity$gene
##D Y <- liver.toxicity$clinic
##D toxi.spls.1 <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
##D                     keepY = c(10, 10, 10))
##D \dontrun{        
##D plot3dVar(toxi.spls.1, rad.in = 0.5, keep.var = TRUE, cex = c(1, 0.8), 
##D           main = "Variables 3D representation") 
##D }
##D toxi.spls.2 <- spls(X, Y, ncomp = 3, keepX = c(10, 10, 10), 
##D                     keepY = c(10, 10, 10))
##D 					  
##D plot3dVar(toxi.spls.2, rad.in = 0.5, Y.label = TRUE, 
##D           main = "Variables 3D representation", 
##D           label.axes.box = "axes")
##D 	
##D ## 3D variable representation for objects of class 'pca'	   
##D data(multidrug)
##D pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
##D 
##D plot3dVar(pca.res)
##D 
##D ## variable representation for objects of class 'splsda'
##D data(liver.toxicity)
##D X <- liver.toxicity$gene
##D Y <- as.factor(liver.toxicity$treatment[, 4])
##D 
##D ncomp <- 3
##D keepX <- rep(20, ncomp)
##D 
##D splsda.liver <- splsda(X, Y, ncomp = ncomp, keepX = keepX)
##D \dontrun{
##D plot3dVar(splsda.liver, var.label = FALSE, Y.label = TRUE, keep.var = TRUE)
##D }
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot3dVar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotIndiv")
### * plotIndiv

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotIndiv
### Title: Plot of Individuals (Experimental Units)
### Aliases: plotIndiv plotIndiv.rcc plotIndiv.pls plotIndiv.spls
###   plotIndiv.plsda plotIndiv.splsda plotIndiv.pca plotIndiv.spca
###   plotIndiv.rgcca plotIndiv.sgcca
### Keywords: multivariate hplot dplot

### ** Examples

## plot of individuals for objects of class 'rcc' 
# ----------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

plotIndiv(nutri.res) #(default) 

col <- rep(c("blue", "red"), c(20, 20))
plotIndiv(nutri.res, ind.names = nutrimouse$diet, col = col)
legend(-2.2, -1.1, c("WT", "PPAR"), pch = c(16, 16), 
       col = c("blue", "red"), text.col = c("blue", "red"),
       cex = 1, pt.cex = c(1.2, 1.2))

## plot of individuals for objects of class 'pls' or 'spls'	
# ----------------------------------------------------   
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
		  
col <- rep(c("blue", "red", "darkgreen", "darkviolet"), rep(16, 4))
cex <- rep(c(1, 1.2, 1, 1.4), c(16, 16, 16, 16))
pch <- rep(c(15, 16, 17, 18), c(16, 16, 16, 16))
plotIndiv(toxicity.spls, comp = 1:2, ind.names = FALSE,
          rep.space = "X-variate", col = col, cex = cex, pch = pch)
legend("topright", c("50 mg/kg", "150 mg/kg", "1500 mg/kg", "2000 mg/kg"), 
       col = c("blue", "red", "darkgreen", "darkviolet"), 
       pch = c(15, 16, 17, 18), pt.cex = c(1, 1.2, 1, 1.4), 
       title = "Treatment")
       
## variable representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------
data(nutrimouse)

# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

#note: the penalty parameters will need to be tuned
wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)


# on the first data set
plotIndiv(wrap.result.sgcca, rep.space = 1, ind.names = TRUE, 
          col = as.numeric(nutrimouse$diet), cex = .6)
# on the second data set
plotIndiv(wrap.result.sgcca, rep.space = 2, ind.names = TRUE, 
          col = as.numeric(nutrimouse$diet), cex = .6)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotIndiv", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotVar")
### * plotVar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotVar
### Title: Plot of Variables
### Aliases: plotVar plotVar.rcc plotVar.pls plotVar.spls plotVar.plsda
###   plotVar.splsda plotVar.pca plotVar.spca plotVar.sgcca plotVar.rgcca
### Keywords: multivariate hplot dplot

### ** Examples

## variable representation for objects of class 'rcc'
# ----------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

plotVar(nutri.res) #(default)

plotVar(nutri.res, comp = 1:2, cutoff = 0.5, 
        X.label = TRUE, Y.label = TRUE)

## variable representation for objects of class 'pls' or 'spls'
# ----------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
	
plotVar(toxicity.spls, keep.var = TRUE, Y.label = TRUE, cex = c(1,0.8))	

## variable representation for objects of class 'splsda'
# ----------------------------------------------------
## Not run: 
##D data(liver.toxicity)
##D X <- liver.toxicity$gene
##D Y <- as.factor(liver.toxicity$treatment[, 4])
##D 
##D ncomp <- 2
##D keepX <- rep(20, ncomp)
##D 
##D splsda.liver <- splsda(X, Y, ncomp = ncomp, keepX = keepX)
##D plotVar(splsda.liver, var.label = FALSE)
## End(Not run)

## variable representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------
data(nutrimouse)

# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

#note: the penalty parameters will need to be tuned
wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
par(mfrow=c(2,2))
plotVar(wrap.result.sgcca, comp = c(1,2), block = c(1,2), ncomp.select = c(1,2))
plotVar(wrap.result.sgcca, comp = c(1,2), block = c(1,2), ncomp.select = c(1,2), labels = TRUE)
plotVar(wrap.result.sgcca, comp = c(1,2), block = c(1,2), ncomp.select = 1, labels = TRUE)
plotVar(wrap.result.sgcca, comp = c(1,2), block = 1, ncomp.select = 1, labels = TRUE)
par(mfrow=c(1,1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotVar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("pls")
### * pls

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pls
### Title: Partial Least Squares (PLS) Regression
### Aliases: pls
### Keywords: regression multivariate

### ** Examples

data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, mode = "classic")

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.pls <- pls(X, Y, ncomp = 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pls", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plsda")
### * plsda

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plsda
### Title: Partial Least Squares Discriminate Analysis (PLS-DA).
### Aliases: plsda
### Keywords: regression multivariate

### ** Examples

## First example
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)
palette(c("red", "blue"))
col.breast <- as.numeric(as.factor(Y))
plotIndiv(plsda.breast, ind.names = TRUE, col = col.breast)
legend('bottomleft', c("After", "Before"), pch = c(16, 16), 
       col = unique(col.breast), cex = 1, pt.cex = c(1.2, 1.2), 
       title = "Treatment")
palette("default")

## Not run: 
##D ## Second example
##D data(liver.toxicity)
##D X <- liver.toxicity$gene
##D Y <- liver.toxicity$treatment[, 4]
##D 
##D plsda.liver <- plsda(X, Y, ncomp = 2)
##D col.rat <- as.numeric(as.factor(Y))
##D plotIndiv(plsda.liver, col = col.rat, ind.names = Y)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plsda", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict")
### * predict

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict
### Title: Predict Method for PLS, sPLS, PLS-DA or sPLS-DA
### Aliases: predict.splsda predict.plsda predict.pls predict.spls
### Keywords: regression multivariate

### ** Examples

data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y, ncomp = 2, mode = "classic")

indiv1 <- c(200, 40, 60)
indiv2 <- c(190, 45, 45)
newdata <- rbind(indiv1, indiv2)
colnames(newdata) <- colnames(X)
newdata

pred <- predict(linn.pls, newdata)

plotIndiv(linn.pls, comp = 1:2, rep.space = "X-variate")
points(pred$variates[, 1], pred$variates[, 2], pch = 19, cex = 1.2)
text(pred$variates[, 1], pred$variates[, 2], 
     c("new ind.1", "new ind.2"), pos = 3)
	 
## First example with plsda
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- as.factor(liver.toxicity$treatment[, 4])

## if training is perfomed on 4/5th of the original data
samp <- sample(1:5, nrow(X), replace = TRUE)  
test <- which(samp == 1)   # testing on the first fold
train <- setdiff(1:nrow(X), test)

plsda.train <- plsda(X[train, ], Y[train], ncomp = 2)
test.predict <- predict(plsda.train, X[test, ], method = "max.dist")
Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
cbind(Y = as.character(Y[test]), Prediction)

## Not run: 
##D ## Second example with splsda
##D splsda.train <- splsda(X[train, ], Y[train], ncomp = 2, keepX = c(30, 30))
##D test.predict <- predict(splsda.train, X[test, ], method = "max.dist")
##D Prediction <- levels(Y)[test.predict$class$max.dist[, 2]]
##D cbind(Y = as.character(Y[test]), Prediction)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.methods")
### * print.methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print
### Title: Print Methods for CCA, (s)PLS, PCA and Summary objects
### Aliases: print print.rcc print.pls print.spls print.summary print.pca
###   print.spca print.rgcca print.sgcca
### Keywords: regression multivariate

### ** Examples

## print for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
print(nutri.res)

## print for objects of class 'summary'
more <- summary(nutri.res, cutoff = 0.65)
print(more)

## print for objects of class 'pls'
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y)
print(linn.pls)

## print for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
print(toxicity.spls)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rcc")
### * rcc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rcc
### Title: Regularized Canonical Correlation Analysis
### Aliases: rcc rcc.default
### Keywords: multivariate

### ** Examples

## Classic CCA
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.res <- rcc(X, Y)

## Regularized CCA
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rcc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("s.match")
### * s.match

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: s.match
### Title: Plot of Paired Coordinates
### Aliases: s.match
### Keywords: multivariate hplot

### ** Examples


## relevant only for canonical mode
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, mode = "canonical", ncomp = 3, 
                      keepX = c(50, 50, 50), keepY = c(10, 10, 10))

color.toxicity <- as.numeric(liver.toxicity$treatment[, 2])
label.toxicity <- liver.toxicity$treatment[, 1]
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)], 
        clabel = 0.5, label = label.toxicity, col = color.toxicity)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("s.match", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("select.var")
### * select.var

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: select.var
### Title: Output of selected variables
### Aliases: select.var select.var.spls select.var.splsda select.var.spca
###   select.var.sipca select.var.sgcca

### ** Examples

data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

# example with sPCA
# ------------------
liver.spca <- spca(X, ncomp = 1, keepX = 10)
select.var(liver.spca, comp = 1)$name
select.var(liver.spca, comp = 1)$value

#example with sIPCA
# -----------------
## Not run: 
##D liver.sipca <- sipca(X, ncomp = 3, keepX = rep(10, 3))
##D select.var(liver.sipca, comp = 1)
## End(Not run)

# example with sPLS
# -----------------
## Not run: 
##D liver.spls = spls(X, Y, ncomp = 2, keepX = c(20, 40),keepY = c(5, 5))
##D select.var(liver.spls, comp = 2)
##D 
##D # example with sPLS-DA
##D data(srbct)   # an example with no gene name in the data
##D X = srbct$gene
##D Y = srbct$class
##D 
##D srbct.splsda = splsda(X, Y, ncomp = 2, keepX = c(5, 10))
##D select = select.var(srbct.splsda, comp = 2)
##D select
##D # this is a very specific case where a data set has no rownames. 
##D srbct$gene.name[substr(select$select, 2,5),]  
## End(Not run)

# example with sGCCA
# -----------------
## Not run: 
##D data(nutrimouse)
##D 
##D # ! need to unmap the Y factor
##D Y = unmap(nutrimouse$diet)
##D data = list(nutrimouse$gene, nutrimouse$lipid,Y)
##D # in this design, gene expression and lipids are connected to the diet factor
##D # and gene expression and lipids are also connected
##D design = matrix(c(0,1,1,
##D                   1,0,1,
##D                   1,1,0), ncol = 3, nrow = 3, byrow = T)
##D #note: the penalty parameters need to be tuned
##D wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), 
##D                                   ncomp = c(2, 2, 1),
##D                                   scheme = "centroid", verbose = FALSE)
##D 
##D select = select.var(wrap.result.sgcca, comp = 1)
##D # loading value of the variables selected on the first block
##D select$value.var[[1]]
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("select.var", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sipca")
### * sipca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sipca
### Title: Independent Principal Component Analysis
### Aliases: sipca
### Keywords: algebra

### ** Examples

data(liver.toxicity)

# implement IPCA on a microarray dataset
sipca.res <- sipca(liver.toxicity$gene, ncomp = 3, mode="deflation", keepX=c(50,50,50))
sipca.res

# samples representation
plotIndiv(sipca.res, ind.names = liver.toxicity$treatment[, 4], cex = 0.5, 
          col = as.numeric(as.factor(liver.toxicity$treatment[, 4])))
## Not run: 
##D plot3dIndiv(sipca.res, cex = 0.01,
##D             col = as.numeric(as.factor(liver.toxicity$treatment[, 4])))
## End(Not run)
# variables representation
plotVar(sipca.res, var.label = TRUE, cex = 0.5)
## Not run: 
##D plot3dVar(sipca.res, rad.in = 0.5, var.label = TRUE, cex = 0.5)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sipca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spca")
### * spca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spca
### Title: Sparse Principal Components Analysis
### Aliases: spca
### Keywords: algebra

### ** Examples

data(liver.toxicity)
spca.rat <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))
spca.rat

## variable representation
plotVar(spca.rat, X.label = TRUE, cex = 0.5)
## Not run: plot3dVar(spca.rat)

## samples representation
plotIndiv(spca.rat, ind.names = liver.toxicity$treatment[, 3], cex = 0.5, 
          col = as.numeric(liver.toxicity$treatment[, 3]))
## Not run: 
##D plot3dIndiv(spca.rat, cex = 0.01, 
##D             col = as.numeric(liver.toxicity$treatment[, 3]))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spls")
### * spls

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spls
### Title: Sparse Partial Least Squares (sPLS)
### Aliases: spls
### Keywords: regression multivariate

### ** Examples

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.spls <- spls(X, Y, ncomp = 2, keepX = c(50, 50), 
                      keepY = c(10, 10))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spls", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("splsda")
### * splsda

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: splsda
### Title: Sparse Partial Least Squares Discriminate Analysis (sPLS-DA)
### Aliases: splsda
### Keywords: regression multivariate

### ** Examples

## First example
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))
palette(c("red", "blue"))
col.breast <- as.numeric(as.factor(Y))
plotIndiv(res, ind.names = TRUE, col = col.breast)
legend('bottomleft', c("After", "Before"), pch = c(16, 16), 
       col = unique(col.breast), cex = 1, pt.cex = c(1.2, 1.2), 
       title = "Treatment")
palette("default")

## Second example
## Not run: 
##D data(liver.toxicity)
##D X <- as.matrix(liver.toxicity$gene)
##D Y <- liver.toxicity$treatment[, 4]
##D 
##D splsda.liver = splsda(X, Y, ncomp = 2, keepX = c(20, 20))
##D col.rat <- as.numeric(as.factor(Y))
##D plotIndiv(splsda.liver, col = col.rat, ind.names = Y)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("splsda", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary")
### * summary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary
### Title: Summary Methods for CCA and PLS objects
### Aliases: summary summary.rcc summary.pls summary.spls
### Keywords: regression multivariate

### ** Examples

## summary for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
more <- summary(nutri.res, cutoff = 0.65)

## summary for objects of class 'pls'
data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y)
more <- summary(linn.pls)

## summary for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))
more <- summary(toxicity.spls, what = "redundancy", keep.var = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tune.multilevel")
### * tune.multilevel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tune.multilevel
### Title: Tuning functions for multilevel analyses
### Aliases: tune.multilevel tune.splsdalevel1 tune.splsdalevel2
###   tune.splslevel
### Keywords: regression multivariate

### ** Examples

## First example: one-factor analysis with sPLS-DA
## Not run: 
##D   data(vac18.simulated) # simulated data
##D   design <- data.frame(sample = vac18.simulated$sample,
##D                        stimu = vac18.simulated$stimulation)
##D   
##D     result.ex1 = tune.multilevel(vac18.simulated$genes,
##D                                design = design,
##D                                ncomp=2,
##D                                test.keepX=c(5, 10, 15), 
##D                                already.tested.X = c(50),
##D                                method = 'splsda',
##D                                dist = 'mahalanobis.dist',
##D                                validation = 'loo') 
##D   
##D   # error rate for the tested parameters est.keepX=c(5, 10, 15)
##D   result.ex1$error
##D   # prediction for ncomp = 2 and keepX = c(50, 15) (15 is the last tested parameter)
##D   result.ex1$prediction.all
##D   table(vac18.simulated$stimulation, result.ex1$prediction.all)
## End(Not run)



## Second example: two-factor analysis with sPLS-DA
## Not run: 
##D   data(liver.toxicity)
##D   dose <- as.factor(liver.toxicity$treatment$Dose.Group)
##D   time <- as.factor(liver.toxicity$treatment$Time.Group)
##D   # note: we made up those data, pretending they are repeated measurements
##D   repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
##D                     6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
##D                     10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
##D                     13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
##D   summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
##D   
##D   design <- data.frame(sample = repeat.indiv,
##D                        dose = dose,
##D                        time = time)
##D   
##D   result.ex2 = tune.multilevel(liver.toxicity$gene,
##D                                 design = design, 
##D                                 ncomp=2,
##D                                 test.keepX=c(5, 10, 15), 
##D                                 already.tested.X = c(50),
##D                                 method = 'splsda',
##D                                 dist = 'mahalanobis.dist') 
##D   result.ex2
## End(Not run)

## Third example: one-factor integrative analysis with sPLS
## Not run: 
##D   data(liver.toxicity)
##D   # note: we made up those data, pretending they are repeated measurements
##D   repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
##D                     6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
##D                     10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
##D                     13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
##D   summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
##D   
##D   # here we are only interested in a one level variation split since spls is an unsupervised method
##D   design <- data.frame(sample = repeat.indiv)
##D   
##D   result.ex3 = tune.multilevel(X = liver.toxicity$gene, Y = liver.toxicity$clinic, 
##D                                 design = design,
##D                                 mode = 'canonical',
##D                                 ncomp=2,
##D                                 test.keepX=c(5, 10, 15), 
##D                                 test.keepY=c(2,3), 
##D                                 already.tested.X = c(50), already.tested.Y = c(5),
##D                                 method = 'spls') 
##D   
##D   result.ex3
##D 
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tune.multilevel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tune.pca")
### * tune.pca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tune.pca
### Title: Tune the number of principal components in PCA
### Aliases: tune.pca
### Keywords: algebra

### ** Examples

data(liver.toxicity)
tune <- tune.pca(liver.toxicity$gene, center = TRUE, scale = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tune.pca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tune.rcc")
### * tune.rcc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tune.rcc
### Title: Estimate the parameters of regularization for Regularized CCA
### Aliases: tune.rcc tune.rcc.default
### Keywords: multivariate dplot

### ** Examples

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

## this can take some seconds
## Not run: 
##D tune.rcc(X, Y, validation = "Mfold")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tune.rcc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("vip")
### * vip

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vip
### Title: Variable Importance in the Projection (VIP)
### Aliases: vip
### Keywords: regression multivariate

### ** Examples

data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological
linn.pls <- pls(X, Y)

linn.vip <- vip(linn.pls)

barplot(linn.vip,
        beside = TRUE, col = c("lightblue", "mistyrose", "lightcyan"),
        ylim = c(0, 1.7), legend = rownames(linn.vip),
        main = "Variable Importance in the Projection", font.main = 4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("vip", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("withinVariation")
### * withinVariation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: withinVariation
### Title: Within matrix decomposition for repeated measurements
###   (cross-over design)
### Aliases: withinVariation
### Keywords: regression multivariate

### ** Examples

## Example: one-factor analysis matrix decomposition
#--------------------------------------------------------------
data(vac18)
X <- vac18$genes
# in design we only need to mention the repeated measurements to split the one level variation
design <- data.frame(sample = vac18$sample)

Xw <- withinVariation(X = X, design = design)
# multilevel PCA
res.pca.1level <- pca(Xw, ncomp = 3)

# compare a normal PCA with a multilevel PCA for repeated measurements.
# note: PCA makes the assumptions that all samples are independent, 
# so this analysis is flawed and you should use a multilevel PCA instead
res.pca <- pca(X, ncomp = 3)

# set up colors for plotIndiv
col.stim <- c("darkblue", "purple", "green4","red3")
col.stim <- col.stim[as.numeric(vac18$stimulation)]

# plotIndiv comparing both PCA and PCA multilevel
plotIndiv(res.pca, ind.names = vac18$stimulation, col = col.stim)
title(main = 'PCA ')
plotIndiv(res.pca.1level, ind.names = vac18$stimulation, col = col.stim)
title(main = 'PCA multilevel')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("withinVariation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("wrapper.rgcca")
### * wrapper.rgcca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: wrapper.rgcca
### Title: mixOmics wrapper for Regularised Generalised Canonical
###   Correlation Analysis (rgcca)
### Aliases: wrapper.rgcca
### Keywords: multivariate

### ** Examples

data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
# with this design, gene expression and lipids are connected to the diet factor
# design = matrix(c(0,0,1,
#                   0,0,1,
#                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#note: the tau parameter is the regularization parameter
wrap.result.rgcca = wrapper.rgcca(data = data, design = design, tau = c(1, 1, 0), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
#wrap.result.rgcca



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("wrapper.rgcca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("wrapper.sgcca")
### * wrapper.sgcca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: wrapper.sgcca
### Title: mixOmics wrapper for Sparse Generalised Canonical Correlation
###   Analysis (sgcca)
### Aliases: wrapper.sgcca
### Keywords: multivariate

### ** Examples

data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid,Y)
# with this design, gene expression and lipids are connected to the diet factor
# design = matrix(c(0,0,1,
#                   0,0,1,
#                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

#note: the penalty parameters will need to be tuned
wrap.result.sgcca = wrapper.sgcca(data = data, design = design, penalty = c(.3,.5, 1), 
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)
wrap.result.sgcca
#did the algo converge?
wrap.result.sgcca$crit  # yes



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("wrapper.sgcca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
