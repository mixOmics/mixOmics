# ----------------------------------------
# testing after the NEW IG update from Jan 2015
# date: 15/02/2015
# latest update: 12/03/2015: testing CIM
# latest update: 13/03/2015 : reformatting for branch multilevel in bitbucket
# ----------------------------------------


# ------ notes for me to compile the package (if need be)
R CMD build --resave-data mixOmics
R CMD INSTALL -l MyR/ mixOmics_5.0-4.tar.gz 
R CMD check mixOmics --as-cran --timings
# --------------------------------

# list of authors in DESCRIPTION file
Author: Sebastien Dejean, Ignacio Gonzalez, Kim-Anh Le Cao with
contributions from Pierre Monget, Jeff Coquery, FangZou Yao,
Benoit Liquet, Florian Rohart, Benoit Gautier

# woud need to put as a Author@R but did not work
# see http://r-pkgs.had.co.nz/description.html help!)
Authors@R: c(
  person("Kim-Anh", "Le Cao", email = "k.lecao@uq.edu.au", role = "cre"),
  person("Ignacio", "Gonzalez", email = "igonzalez@toulouse.inra.fr", role = "aut"),
  person("Sébastien", "Déjean", email = "sebastien.dejean@math.univ-toulouse.fr", role = "aut"),
  person("Florian", "Rohart", email = "k.lecao@uq.edu.au", role = "aut"),
  person("Benoit", "Gautier", email = "k.lecao@uq.edu.au", role = "aut"),
  person("Benoit", "Liquet", email = "b.liquet@uq.edu.au", role = "ctb", comment = 'contribution to the multilevel module'),
  person("FangZou", "Yao", role = "ctb", comment = 'contribution to the IPCA module')
  person("Pierre", "Monget", role = "ctb", comment = 'contribution to the sPLS-DA module')
  person("Jeff", "Coquery", role = "ctb", comment = 'contribution to the IPCA module')
)


sessionInfo()
# R version 3.1.0 (2014-04-10)
# Platform: x86_64-apple-darwin10.8.0 (64-bit)
# 
# locale:
#   [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] tools_3.1.0


# ------------------------------------------------------------------------
# -------summary of old changes from BEFORE Jan 2015 ---------------------
# ------------------------------------------------------------------------
# changes that were made by Ignacio (1/12/2014) and need to be checked:
# - new code for multilevel() function to sPLS and sPLS-DA methods. 

# - new code for Split.variation.one.level() and Split.variation.two.level() 
  # functions.

# - internal functions multilevel.spls(), multilevel.splsda() and check.one.level() 
  # have been removed.

# - pheatmap clustering S3 method have been removed and also the dependence to the
  # pheatmap package.

# - the vac18 data has been renamed vac18.study.

# - the data.simu data has been renamed vac18.simulated. Components has been 
  # renamed like in vac18 data and the 'time' factor has been added.

# -------------
# changes made so far
# ------------
# 1 - revert back to using the split.variation functions in the multilevel function
# 2 - created a generic split.variation function that includes both Split.variation.one.level and Split.variation.two level
# 3 - to do: reinclude heatmaps
# 4 - to do: test tune.multilevel


# ------------------------------------------------------------------------
# -------summary of my new testing/changes changes AFTER Jan 2015 --------
# ------------------------------------------------------------------------
# - see list of changes in ../../Summary_Changed_Multilevel.rtf

# ---------------
# load data
# ---------------
load('mixomics/data/vac18.simulated.rda')
attributes(vac18.simulated)
# $names
# [1] "genes"       "stimulation" "sample"      "time"       

load('mixomics/data/vac18.rda')

source('mixOmics/R/multilevel.R'); source('mixOmics/R/withinVariation.R')
# this is after update with spls/pls with ... removed
source('mixOmics/R/pls.R'); source('mixOmics/R/plsda.R')
source('mixOmics/R/spls.R'); source('mixOmics/R/splsda.R')


# ------------------------------------------
# testing examples for help file multilevel.Rd
# ------------------------------------------

## First example: one-factor analysis with sPLS-DA, selecting a subset of variables
# as in the paper Liquet et al.
#--------------------------------------------------------------
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample, 
                     stimu = vac18$stimulation)

# multilevel sPLS-DA model
res.1level <- multilevel(X, ncomp = 3, design = design,
                         method = "splsda", keepX = c(30, 137, 123))

# set up colors for plotIndiv
col.stimu <- c("darkblue", "purple", "green4","red3")
col.stimu <- col.stimu[as.numeric(Y)]
plotIndiv(res.1level, ind.names = Y, col = col.stimu)

## Second example: two-factor analysis with sPLS-DA, selecting a subset of variables
# as in the paper Liquet et al.
#--------------------------------------------------------------
# set to dontrun{} for elapse time
\dontrun{
data(vac18.simulated) # simulated data

X <- vac18.simulated$genes
design <- data.frame(sample = vac18.simulated$sample,
                     stimu = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

res.2level <- multilevel(X, ncomp = 2, design = design,
                         keepX = c(200, 200), method = 'splsda')

# set up colors and pch for plotIndiv
col.stimu <- as.numeric(design$stimu)
pch.time <- c(20, 4)[as.numeric(design$time)]

plotIndiv(res.2level, col = col.stimu, ind.names = FALSE,
          pch = pch.time)
legend('bottomright', legend = levels(as.factor(design$stimu)),
       col = unique(col.stimu), pch = 20, cex = 0.8, 
       title = "Stimulation")
legend('topright', col = 'black', legend = levels(design$time),  
       pch = unique(pch.time), cex = 0.8, title = "Time")
}

## Third example: one-factor analysis with sPLS, selecting a subset of variables
#--------------------------------------------------------------
# set to dontrun{} for elapse time
\dontrun{
data(liver.toxicity)
# note: we made up those data, pretending they are repeated measurements
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                 6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                 10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                 13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each

design <- data.frame(sample = repeat.indiv, 
                     stimu = liver.toxicity$treatment$Dose.Group)
res.spls.1level <- multilevel(X = liver.toxicity$gene,
                                       Y=liver.toxicity$clinic,
                                       design = design,
                                       ncomp = 3,
                                       keepX = c(50, 50, 50), keepY = c(5, 5, 5),
                                       method = 'spls', mode = 'canonical')

# set up colors and pch for plotIndiv
col.stimu <- as.numeric(as.factor(design$stimu))

plotIndiv(res.spls.1level, rep.space = 'X-variate', ind.names = FALSE, 
          col = col.stimu, pch = 20)
title(main = 'Gene expression variates')
plotIndiv(res.spls.1level, rep.space = 'Y-variate', ind.names = FALSE,
          col = col.stimu, pch = 20)
title(main = 'Clinical measurements variates')
legend('bottomright', legend = levels(as.factor(design$stimu)),
       col = unique(col.stimu), pch = 20, cex = 0.8, 
       title = "Dose")
}

# -------------------------------------------------
# end testing multilevel.Rd examples
# -------------------------------------------------

# ------------------------------------------
# testing examples for help file withinVariation.Rd
# ------------------------------------------

## Example: one-factor analysis matrix decomposition
#--------------------------------------------------------------
data(vac18)
X <- vac18$genes
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)
Xw <- withinVariation(X = X, design = design)
# multilevel PCA
res.pca.1level <- pca(Xw, ncomp = 3)

# compare a normal PCA with a multilevel PCA for repeated measurements.
# note: PCA makes the assumptions that all samples are independent, so this analysis is flawed and you should use a multilevel PCA instead
res.pca <- pca(X, ncomp = 3)

# set up colors for plotIndiv
col.stim <- c("darkblue", "purple", "green4","red3")
col.stim <- col.stim[as.numeric(vac18$stimulation)]

# plotIndiv comparing both PCA and PCA multilevel
plotIndiv(res.pca, ind.names = vac18$stimulation, col = col.stim)
title(main = 'PCA ')
plotIndiv(res.pca.1level, ind.names = vac18$stimulation, col = col.stim)
title(main = 'PCA multilevel')

# -------------------------------------------------
# end testing multilevel.Rd examples
# -------------------------------------------------


# ------------------------------------------
# testing examples for help file tune.multilevel.Rd
# ------------------------------------------
load('mixomics/data/vac18.rda')

source('mixOmics/R/multilevel.R'); source('mixOmics/R/withinVariation.R')
source('mixOmics/R/tune.multilevel.R'); source('mixOmics/R/tune.splsdalevel1.R'); source('mixOmics/R/tune.splsdalevel2.R')
source('mixOmics/R/tune.splslevel.R')

# this is after update with spls/pls with ... removed
source('mixOmics/R/pls.R'); source('mixOmics/R/plsda.R')
source('mixOmics/R/spls.R'); source('mixOmics/R/splsda.R')
load('mixomics/data/vac18.simulated.rda')
attributes(vac18.simulated)

## First example: one-factor analysis with sPLS-DA
\dontrun{
  data(vac18.simulated) # simulated data
  design <- data.frame(sample = vac18.simulated$sample,
                       stimu = vac18.simulated$stimulation)
  
    result.ex1 = tune.multilevel(vac18.simulated$genes,
                               design = design,
                               ncomp=2,
                               test.keepX=c(5, 10, 15), 
                               already.tested.X = c(50),
                               method = 'splsda',
                               dist = 'mahalanobis.dist',
                               validation = 'loo') 
  
  # error rate for the tested parameters est.keepX=c(5, 10, 15)
  result.ex1$error
  # prediction for ncomp = 2 and keepX = c(50, 15) (15 is the last tested parameter)
  result.ex1$prediction.all
  table(vac18.simulated$stimulation, result.ex1$prediction.all)
}



## Second example: two-factor analysis with sPLS-DA
\dontrun{
  data(liver.toxicity)
  dose <- as.factor(liver.toxicity$treatment$Dose.Group)
  time <- as.factor(liver.toxicity$treatment$Time.Group)
  # note: we made up those data, pretending they are repeated measurements
  repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                    6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                    10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                    13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
  summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
  
  design <- data.frame(sample = repeat.indiv,
                       dose = dose,
                       time = time)
  
  result.ex2 = tune.multilevel(liver.toxicity$gene,
                                design = design, 
                                ncomp=2,
                                test.keepX=c(5, 10, 15), 
                                already.tested.X = c(50),
                                method = 'splsda',
                                dist = 'mahalanobis.dist') 
  result.ex2
}

## Third example: one-factor integrative analysis with sPLS
\dontrun{
  data(liver.toxicity)
  dose <- as.factor(liver.toxicity$treatment$Dose.Group)
  # note: we made up those data, pretending they are repeated measurements
  repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                    6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                    10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                    13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
  summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
  
  design <- data.frame(sample = repeat.indiv,
                       dose = dose)
  
  result.ex3 = tune.multilevel(X = liver.toxicity$gene, Y = liver.toxicity$clinic, 
                                design = design,
                                mode = 'canonical',
                                ncomp=2,
                                test.keepX=c(5, 10, 15), 
                                test.keepY=c(2,3), 
                                already.tested.X = c(50), already.tested.Y = c(5),
                                method = 'spls') 
  
  result.ex3
}

# ------------------------------------------
# end testing examples for help file tune.multilevel.Rd
# ------------------------------------------

# ------------------------------------------
# testing examples for help file pheatmap.multileve.Rd
# ------------------------------------------
load('mixomics/data/vac18.rda')

source('mixOmics/R/multilevel.R'); source('mixOmics/R/withinVariation.R')
source('mixOmics/R/pheatmap.multilevel.R')
source('mixOmics/R/pheatmap.multilevel.splsda1fact.R')
source('mixOmics/R/pheatmap.multilevel.splsda2fact.R')
require(pheatmap)

load('mixomics/data/vac18.simulated.rda')
attributes(vac18.simulated)


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
col.stimu = unique(col.stimu)
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
\dontrun{
data(vac18.simulated) 

X <- vac18.simulated$genes
design <- data.frame(sample = vac18.simulated$sample,
                     stimul = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

vac18.splsda2.multilevel <- multilevel(X, ncomp = 2, design = design,
                            keepX = c(200, 200), method = 'splsda')

# set up colors for each level of pheatmap 
col.sample <- c("lightgreen", "red","lightblue","darkorange","purple","maroon") # 6 samples
col.time <- c("pink","lightblue1") # two time points
col.stimu <- c('green', 'black', 'red', 'blue') # 4 stimulations
# set up labels for the 2 levels in design matrix
label.stimu <- unique(design[, 2])
label.time <- unique(design$time)

pheatmap.multilevel(vac18.splsda2.multilevel,
                                # colors:
                                col_sample=col.sample, 
                                col_stimulation=col.stimu, 
                                col_time=col.time,
                                #labels for each level
                                label_color_stimulation=label.stimu,
                                label_color_time=label.time, 
                                #clustering method
                                clustering_method="ward",
                                #show col/row names and font size
                                show_colnames = FALSE,
                                show_rownames = TRUE,
                                fontsize_row=2)
}

# ======== that's it for now, KA 14/03/2015 ======================











# 11.2 multilevel PCA
# ----------------------
rm(list=ls())
source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')
library(mixOmics)
data(vac18)

# testing example from the workshop material (with modifications)
X <- vac18$genes
Y <- vac18$stimulation

design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)

# 11.2.1
#KA# I think I disagree with the fact that withinVariation does not have a 'design' argument.
# ?! to be changed?!#
Xw_IG <- withinVariation(X, design = design)
dim(Xw_IG)  #[1]   42 2500

Xw_mixomics <- Split.variation.one.level(X, Y, sample = vac18$sample)
dim(Xw_mixomics$Xw)  #[1]   42 2500

names(Xw_IG)
names(Xw_mixomics)

all.equal(Xw_IG, Xw_mixomics$Xw)
### Dicrepancy found between IG update and mixomics
# missing outputs (Xb, Xm) in IG update ??? should be added ???

pca.multilevel_IG <- pca(Xw_IG, ncomp = 3, scale = TRUE, center = TRUE)
pca.multilevel_IG
pca.multilevel_mixomics <- pca(Xw_mixomics$Xw, ncomp = 3, scale = TRUE, center = TRUE)
pca.multilevel_mixomics

names(pca.multilevel_IG)
names(pca.multilevel_mixomics)

all.equal(pca.multilevel_IG, pca.multilevel_mixomics)
for (i in names(pca.multilevel_mixomics)[!names(pca.multilevel_mixomics) %in% c("call")]){
  if (!all.equal(pca.multilevel_IG[[i]], pca.multilevel_mixomics[[i]]) == TRUE) {
    stop(i)
  }
}
### Dicrepancy found between update and mixomics
# call ok

# -> ok results same as in tutorial

# 11.2.2 Sample plots (with current existing plot functions)
col.stimu = as.numeric(Y)
plotIndiv(pca.multilevel_IG, col = col.stimu, ind.names = vac18$sample)
title('VAC18, multilevel PCA, comp 1 - 2')
legend("topleft", col = unique(col.stimu), pch = 16,
 legend = c("LIPO5", "GAG+", "GAG-", "NS"))

plotIndiv(pca.multilevel_mixomics, col = col.stimu, ind.names = vac18$sample)
title('VAC18, multilevel PCA, comp 1 - 2')
legend("topleft", col = unique(col.stimu), pch = 16,
       legend = c("LIPO5", "GAG+", "GAG-", "NS"))
# -> ok results same as in tutorial

# 11.3.1 multilevel PLS-DA
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)

vac18.plsda.multilevel_IG <- multilevel(X, ncomp = 3, design = design,
                                        method = "splsda", keepX = c(ncol(X), ncol(X), ncol(X)))
rm(multilevel)

vac18.plsda.multilevel_mixomics <- multilevel(X, cond = Y, ncomp = 3, tab.prob.gene = vac18$tab.prob.gene,
                                              sample = vac18$sample, method = "splsda",
                                              keepX = c(ncol(X), ncol(X), ncol(X)))

names(vac18.plsda.multilevel_IG)
names(vac18.plsda.multilevel_mixomics)

for (i in names(vac18.plsda.multilevel_mixomics)[!names(vac18.plsda.multilevel_mixomics) %in% c("call", "sample", "name.condition", "tab.prob.gene")]){
  if (!all.equal(vac18.plsda.multilevel_IG[[i]], vac18.plsda.multilevel_mixomics[[i]]) == TRUE) {
    stop(i)
  }
}
### Dicrepancy found between IG update and mixomics
# call: ok 
# sample: ok input / output present only in multilevel (mixOmics)
# name.condition: output missing in multilevel (update_IG) ??? should be added ???
# tab.prob.gene: ok input / output present only in multilevel (mixOmics)

# Sample plots (with current existing plot functions)
plotIndiv(vac18.plsda.multilevel_IG, col = col.stimu, ind.names = vac18$sample)
title('VAC18, multilevel PLS-DA, comp 1 - 2')

plotIndiv(vac18.plsda.multilevel_mixomics, col = col.stimu, ind.names = vac18$sample)
title('VAC18, multilevel PLS-DA, comp 1 - 2')
# -> plot slightly different from tutorial as the components were inverted, but ok # *BG* ??? not anymore fix pb by keeping the order of the levels in the design matrix 

# 11.3.2 multilevel sPLS-DA: note for tutorial, need to change name spls.multilevel so splsda.multilevel
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)

source('multilevel_new/multilevel.R');
vac18.splsda.multilevel_IG <- multilevel(X, ncomp = 3, design = design,
                                     method = "splsda", keepX = c(30, 137, 123))
rm(multilevel)

vac18.splsda.multilevel_mixomics <- multilevel(X, cond = Y, ncomp = 3, tab.prob.gene = vac18$tab.prob.gene,
                                             sample = vac18$sample, method = "splsda",
                                             keepX = c(30, 137, 123))

names(vac18.splsda.multilevel_IG)
names(vac18.splsda.multilevel_mixomics)

for (i in names(vac18.splsda.multilevel_mixomics)[!names(vac18.splsda.multilevel_mixomics) %in% c("call", "sample", "name.condition", "tab.prob.gene")]){
  if (!all.equal(vac18.splsda.multilevel_IG[[i]], vac18.splsda.multilevel_mixomics[[i]]) == TRUE) {
    stop(i)
  }
}
### Dicrepancy found between IG update and mixomics
# call: ok 
# sample: ok input / output present only in multilevel (mixOmics)
# name.condition: output missing in multilevel (update_IG) ??? should be added ???
# tab.prob.gene: ok input / output present only in multilevel (mixOmics)

# Sample plots (with current existing plot functions)
plot3dIndiv(vac18.splsda.multilevel_IG , ind.names = Y, col = col.stimu, axes.box = "both")
plot3dIndiv(vac18.splsda.multilevel_mixomics , ind.names = Y, col = col.stimu, axes.box = "both")

# -> ok results seem the same as in tutorial

col.samp <- c("lightgreen", "red", "lightblue", "darkorange",
              "purple", "maroon", "blue", "chocolate", "turquoise",
              "tomato1", "pink2", "aquamarine")
col.stimu2 = unique(col.stimu)
# note for tutorial: need to change name spls.multilevel so splsda.multilevel

source('multilevel_new/pheatmap.multilevel.splsda1fact.R'); require(pheatmap)
pheatmap.multilevel.splsda1fact(vac18.splsda.multilevel_IG, clustering_method = "ward",
                    col_sample = col.samp, col_stimulation = col.stimu2,
                    label_annotation = c("Subject", "Stimulus"),
                    fontsize = 8, border = FALSE, fontsize_row = 3,
                    show_colnames = FALSE,
                    show_rownames = TRUE, 
                    fontsize_col = 2, width = 10)
rm(pheatmap.multilevel.splsda1fact); detach(package:pheatmap)
# error:
# Error in `row.names<-.data.frame`(`*tmp*`, value = value) : invalid 'row.names' length
## This one fails, even after I added the extr class splsda1fact. TO BE FIXED
# I think it is due to the tab.prob.gene in the pheatmap.multilevel.splsda1fact (that way of extract the tab.prob.gene on l53 is really badly coded by the way:
# ### visualisation of the genes instead of the probes
# if(!(is.null(result$tab.prob.gene))) geneX <- result$tab.prob.gene[match(probeX,result$tab.prob.gene[,1]),2]
# )
#=> the pheatmap needs to be changed

pheatmap.multilevel(vac18.splsda.multilevel_mixomics, clustering_method = "ward",
                    col_sample = col.samp, col_stimulation = col.stimu2,
                    label_annotation = c("Subject", "Stimulus"),
                    fontsize = 8, border = FALSE, fontsize_row = 3,
                    show_colnames = FALSE,
                    show_rownames = TRUE, 
                    fontsize_col = 2, width = 10)

# #before, version 5.0-3
# rm(multilevel)
# library(mixOmics)
# vac18.spls.multilevel.BEFORE <- multilevel(X, cond = Y, ncomp = 3,
#                                     tab.prob.gene = vac18$tab.prob.gene,
#                                     sample = vac18$sample, method = "splsda",
#                                     keepX = c(30, 137, 123))
# col.samp <- c("lightgreen", "red", "lightblue", "darkorange",
#               "purple", "maroon", "blue", "chocolate", "turquoise",
#               "tomato1", "pink2", "aquamarine")
# col.stimu2 = unique(col.stimu)
# pheatmap.multilevel(vac18.spls.multilevel.BEFORE, clustering_method = "ward",
#                     col_sample = col.samp, col_stimulation = col.stimu2,
#                     label_annotation = c("Subject", "Stimulus"),
#                     fontsize = 8, border = FALSE, fontsize_row = 3,
#                     show_colnames = FALSE,
#                     #show_rownames = TRUE, 
#                     fontsize_col = 2, width = 10)


# 11.4 multilevel PLS
# back to the update and testing tutorial
rm(list=ls())
source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')

data(liver.toxicity)
repeat.indiv = c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                   6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                   10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                   13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
summary(as.factor(repeat.indiv))

design <- data.frame(sample = repeat.indiv, 
                     stimul = liver.toxicity$treatment$Dose.Group)
# Note for tutorial: would need to be renamed as liver.spls.multillevel!
liver.spls.multilevel_IG <- multilevel(X = liver.toxicity$gene,
                                      Y=liver.toxicity$clinic,
                                      design = design,
                                      ncomp = 3,
                                      keepX = c(50, 50, 50), keepY = c(5, 5, 5),
                                      method = 'spls', mode = 'regression')
rm(multilevel)

liver.spls.multilevel_mixomics <- multilevel(X = liver.toxicity$gene, Y=liver.toxicity$clinic, 
                                             cond = liver.toxicity$treatment$Dose.Group, 
                                             sample = repeat.indiv, ncomp = 3, keepX = c(50, 50, 50), 
                                             keepY = c(5, 5, 5), method = 'spls', mode = 'regression') 


### Dicrepancy found between update and mixomics
# call: ok 
# mode: !!! hard coding !!! in multilevel.spls function
liver.spls.multilevel_mixomics$mode ### always canonical mode is performed => see multilevel.spls function 4th line

### Check done using canonical mode
source("multilevel_new/multilevel.R")
liver.spls.multilevel_IG <- multilevel(X = liver.toxicity$gene,
                                       Y=liver.toxicity$clinic,
                                       design = design,
                                       ncomp = 3,
                                       keepX = c(50, 50, 50), keepY = c(5, 5, 5),
                                       method = 'spls', mode = 'canonical')
rm(multilevel)

liver.spls.multilevel_mixomics <- multilevel(X = liver.toxicity$gene, Y=liver.toxicity$clinic, 
                                             cond = liver.toxicity$treatment$Dose.Group, 
                                             sample = repeat.indiv, ncomp = 3, keepX = c(50, 50, 50), 
                                             keepY = c(5, 5, 5), method = 'spls', mode = 'canonical') 

names(liver.spls.multilevel_IG)
names(liver.spls.multilevel_mixomics)

for (i in names(liver.spls.multilevel_mixomics)[!names(liver.spls.multilevel_mixomics) %in% c("call", "sample", "name.condition")]){
  if (!all.equal(liver.spls.multilevel_IG[[i]], liver.spls.multilevel_mixomics[[i]]) == TRUE) {
    stop(i)
  }
}

### Dicrepancy found between IG update and mixomics
# call: ok 
# sample: ok input / output present only in multilevel (mixOmics)
# name.condition: output missing in multilevel (update_IG) ??? should be added ???

# the following 'runs' but would need to be compared with the before and after update.
plotVar(liver.spls.multilevel_IG, comp = 1:2, X.label = TRUE, Y.label = TRUE, cex = c(0.5, 0.9))
plotVar(liver.spls.multilevel_mixomics, comp = 1:2, X.label = TRUE, Y.label = TRUE, cex = c(0.5, 0.9))

CIM_IG <- cim(liver.spls.multilevel_IG, comp = 1:2, xlab = "genes",
            ylab = "clinic var",
            margins = c(5, 6), zoom = FALSE)

CIM_mixomics <- cim(liver.spls.multilevel_mixomics, comp = 1:2, xlab = "genes",
              ylab = "clinic var",
              margins = c(5, 6), zoom = FALSE)

names(CIM_IG)
names(CIM_mixomics)

for (i in names(CIM_mixomics)){
  if (!all.equal(CIM_IG[[i]], CIM_mixomics[[i]]) == TRUE) {
    stop(i)
  }
}
### CIM_IG = CIM_mixomics

network(liver.spls.multilevel_IG, comp = 1:2, threshold = 0.8,
         Y.names = NULL, keep.var = TRUE,
         color.node = c( "lightcyan","mistyrose"),
         shape.node = c("circle", "rectangle"),
         color.edge = c("red", "green"),
         lty.edge = c("solid", "solid"), lwd.edge = c(1, 1),
          show.edge.labels = FALSE, interactive = FALSE)

network(liver.spls.multilevel_mixomics, comp = 1:2, threshold = 0.8,
        Y.names = NULL, keep.var = TRUE,
        color.node = c( "lightcyan","mistyrose"),
        shape.node = c("circle", "rectangle"),
        color.edge = c("red", "green"),
        lty.edge = c("solid", "solid"), lwd.edge = c(1, 1),
        show.edge.labels = FALSE, interactive = FALSE)


# ----------------------
# CIM
# ---------------------
# new CIM IG
rm(list=ls())
source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')

# --------
## First example: one-factor analysis with sPLS-DA
# --------
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
#make sure sample comes first in the help file?
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)

res_IG <- multilevel(X, ncomp = 3, design = design, method = "splsda", keepX = c(30, 137, 123))


source('multilevel_new/CIM/cim.R')
##source('multilevel_new/CIM/cim.mlspls.R')
source('multilevel_new/CIM/cim.mlsplsda.R')
source('multilevel_new/CIM/color.palettes.R')  # !!! warning message for In grepl("\n", lines, fixed = TRUE) :input string 2 is invalid in this locale
source('multilevel_new/CIM/isColor.R')
source('multilevel_new/CIM/imageMap.R')

stim.col <- c("darkblue", "purple", "green4","red3")
stim.col <- stim.col[as.numeric(design$stimul)]

pdf('testCIMlvl1.pdf')
cim(res_IG, sample.sideColors = stim.col, sample.names = design$stimul, var.names = FALSE)
dev.off()

# -----------
# second example with: two factor analysis with sPLS-DA
# ------------
load('multilevel_new/vac18.simulated.rda')

X <- vac18.simulated$genes
design <- data.frame(samp = vac18.simulated$sample,
                     stim = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

res.2level <- multilevel(X, ncomp = 2, design = design,
                         keepX = c(120, 10), method = 'splsda')

stim.col <- c("darkblue", "purple", "green4","red3")
stim.col <- stim.col[as.numeric(design$stim)]
time.col <- c("orange", "cyan")[as.numeric(design$time)]

pdf('testCIMlvl2.pdf')
cim(res.2level, sample.sideColors = cbind(stim.col, time.col), 
    sample.names = paste(design$time, design$stim, sep = "_"),
    var.names = FALSE)
dev.off()

# \method{cim}{mlspls}(object, comp = 1:object$ncomp, color = NULL,
#                      X.var.names = TRUE, Y.var.names = TRUE, sample.names = TRUE,
#                      x.sideColors = NULL, y.sideColors = NULL, sample.sideColors = NULL, 
#                      mapping = "XY", dist.method = c("euclidean", "euclidean"),
#                      clust.method = c("complete", "complete"), cluster = "both", center = TRUE, 
#                      scale = FALSE, \ldots)


# ----------
# 3rd example with mlspls??
# ---------

# --------------------------------------------------
# lets compare now with a pheatmap output
# -------------------------------------------------
source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')

# --------
## First example: one-factor analysis with sPLS-DA
# --------
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
#make sure sample comes first in the help file?
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)

res_IG <- multilevel(X, ncomp = 3, design = design, method = "splsda", keepX = c(30, 137, 123))


col.samp <- c("lightgreen", "red", "lightblue", "darkorange",
              "purple", "maroon", "blue", "chocolate", "turquoise",
              "tomato1", "pink2", "aquamarine")
# color for plotIndiv
stimu <- as.numeric(vac18$stimulation)

# mixOmics palette (not merged yet)
mixo.gray = gray.colors(1, start = 0.76, gamma = 1)

mixo.col = c('#388ECC', # mixOmics logo blue
             '#F68B33', # mixOmics logo orange
             mixo.gray, # mixOmics logo grey
             '#009E73', # shiny dark green
             '#CC79A7', # shiny purple/pink
             '#F0E442', #shiny yellow
             'black',
             '#D55E00', #shiny dark orange
             '#0072B2', #shiny dark blue
             '#999999'  # shiny grey
             #'#E69F00', # shiny orange
             #'#56B4E9' #Shiny blue
)
col.stimu <- mixo.col[1:4]  ##** need to put the mixOmics colors palette using the acutal function (not merged yet)
col.stimu <- col.stimu[stimu]
col.stimu2 = unique(col.stimu)
# note for tutorial: need to change name spls.multilevel so splsda.multilevel

source('multilevel_new/pheatmap.multilevel.splsda1fact.R'); require(pheatmap)

pdf('example-pheatmaplvl1.pdf')
pheatmap.multilevel.splsda1fact(res_IG, 
                                col_sample = col.samp, col_stimulation = col.stimu2,
                                scale = 'row',
                                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
                                clustering_method = "complete",
                                label_annotation = c("Subject", "Stimulus"),
                                fontsize = 8, border = FALSE, fontsize_row = 3,
                                show_colnames = FALSE,
                                show_rownames = FALSE, 
                                fontsize_col = 2, width = 10)
dev.off()


source('multilevel_new/CIM/cim.R')
##source('multilevel_new/CIM/cim.mlspls.R')
source('multilevel_new/CIM/cim.mlsplsda.R')
source('multilevel_new/CIM/color.palettes.R')  # !!! warning message for In grepl("\n", lines, fixed = TRUE) :input string 2 is invalid in this locale
source('multilevel_new/CIM/isColor.R')
source('multilevel_new/CIM/imageMap.R')
col.stimu.cim <- c("darkblue", "purple", "green4","red3")
col.stimu.cim <- col.stimu.cim[as.numeric(design$stimul)]


# not possible to add the sample information
pdf('testCIMlvl1.pdf')
cim(res_IG, comp = c(1:3), sample.sideColors = col.stimu.cim, sample.names = design$stimul, var.names = FALSE,
   clust.method = c('complete', 'complete'), transpose = TRUE)
dev.off()
# the "ward" method has been renamed to "ward.D"; note new "ward.D2"



rm(pheatmap.multilevel.splsda1fact); detach(package:pheatmap)
# error:
# Error in `row.names<-.data.frame`(`*tmp*`, value = value) : invalid 'row.names' length
## This one fails, even after I added the extr class splsda1fact. TO BE FIXED
# I think it is due to the tab.prob.gene in the pheatmap.multilevel.splsda1fact (that way of extract the tab.prob.gene on l53 is really badly coded by the way:
# ### visualisation of the genes instead of the probes
# if(!(is.null(result$tab.prob.gene))) geneX <- result$tab.prob.gene[match(probeX,result$tab.prob.gene[,1]),2]
# )
#=> the pheatmap needs to be changed

pheatmap.multilevel(vac18.splsda.multilevel_mixomics, clustering_method = "ward",
                    col_sample = col.samp, col_stimulation = col.stimu2,
                    label_annotation = c("Subject", "Stimulus"),
                    fontsize = 8, border = FALSE, fontsize_row = 3,
                    show_colnames = FALSE,
                    show_rownames = TRUE, 
                    fontsize_col = 2, width = 10)











# ---------------------------------
# end of testing for now
# ---------------------------------




# ----- 
# why has the pheatmap function been removed? 
# I want to put it back otherwise it is impossible to provide such code (see below) to users
# Benoit_pheatmap
# -----
library(pheatmap)

select.comp1 = select.var(res, comp = 1)$name
length(select.comp1) # 30
mat.pheatmap = t(res$Xw[, select.comp1])
dim(mat.pheatmap)  #30 42


# colors for sample and stimulation
col_sample = sample(colors())[unique(as.numeric(vac18$sample))]
names(col_sample) = unique(as.numeric(vac18$sample))

col_stimu <- c("darkblue", "purple", "green4","red3")
names(col_stimu) = levels(vac18$stimulation)
# set up annotation colors
annot_col = list(Sample = col_sample, Stimulation = col_stimu)

#set up annotation to match the colors
sample = factor(vac18$sample, levels =  1:length(unique(vac18$sample)))
names(sample) = colnames(mat.pheatmap)
annotation <- data.frame(Sample=sample, Stimulation=as.character(vac18$stimulation))
rownames(annotation) <- colnames(mat.pheatmap) 


## note this does not run anymore (it used to!), I get a palette color issue
# The "ward" method has been renamed to "ward.D"; note new "ward.D2"
# The "ward" method has been renamed to "ward.D"; note new "ward.D2"
# Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : 
#   invalid graphics state

pheatmap(mat.pheatmap,
         annotation = annotation, annotation_colors = annot_col, annotation_legend = TRUE,
         #label_annotation=NULL,
                    border=FALSE,
                    clustering_method="ward",
                    show_colnames = TRUE,
                    show_rownames = TRUE)



# -------------
# other outputs need to be tested, e.g. plotVar
# plot3D, pltVar3D
# ------------

# ==================
## Second example: two-factor analysis with sPLS-DA
# =================
rm(list=ls())
source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')
require(mixOmics)

load('multilevel_new/vac18.simulated.rda')
attributes(vac18.simulated)
# data(vac18.simulated) # simulated data

X <- vac18.simulated$genes
design <- data.frame(sample = vac18.simulated$sample,
                     stimul = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

res.2level_IG <- multilevel(X, ncomp = 2, design = design,
                            keepX = c(200, 200), method = 'splsda')
rm(multilevel)


res.2level_mixomics <- multilevel(X, cond = design[, -1], sample = design[, 1], ncomp = 2,
                                  keepX = c(200, 200), tab.prob.gene = NULL,
                                  method = 'splsda')

names(res.2level_IG)
names(res.2level_mixomics)

for (i in names(res.2level_mixomics)[!names(res.2level_mixomics) %in% c("call", "sample", "name.condition", "name.time")]){
  if (!all.equal(res.2level_IG[[i]], res.2level_mixomics[[i]]) == TRUE) {
    stop(i)
  }
}

### Dicrepancy found between IG update and mixomics
# call: ok 
# sample: ok output present only in multilevel (mixOmics)
# name.condition: output missing in multilevel (update_IG) ??? should be added ???
# name.time: output missing in multilevel (update_IG) ??? should be added ???
# tab.prob.gene: ok output present only in multilevel (mixOmics)

# color and pch for plotIndiv
col.stimul <- as.numeric(design$stimul)
pch.time <- c(20, 4)[as.numeric(design$time)]

plotIndiv(res.2level_IG, col = col.stimul, ind.names = FALSE,pch = pch.time)
legend('bottomright', legend = levels(design$stimul),
       col = unique(col.stimul), pch = 20, cex = 0.8, 
       title = "Stimulation")
legend('topright', col = 'black', legend = levels(design$time),  
       pch = unique(pch.time), cex = 0.8, title = "Time")

plotIndiv(res.2level_mixomics, col = col.stimul, ind.names = FALSE,pch = pch.time)
legend('bottomright', legend = levels(design$stimul),
       col = unique(col.stimul), pch = 20, cex = 0.8, 
       title = "Stimulation")
legend('topright', col = 'black', legend = levels(design$time),  
       pch = unique(pch.time), cex = 0.8, title = "Time")



# -----------------
# need to do the same schmilblick for pheatmap
# I want the function to be added again.
# ----------------

# color for plotIndiv
col.sample = c("lightgreen", "red","lightblue","darkorange","purple","maroon") # 6 samples
col.time = c("pink","lightblue1") # two time points
col.stimu = c('green', 'black', 'red', 'blue') # 4 stimulations
label.stimu = unique(design[, 2])
label.time = unique(design$time)

source('multilevel_new/pheatmap.multilevel.splsda2fact.R'); require(pheatmap)
pheatmap.multilevel.splsda2fact(res.2level_IG, 
                                col_sample=col.sample, 
                                col_stimulation=col.stimu, 
                                col_time=col.time,
                                label_color_stimulation=label.stimu,
                                label_color_time=label.time, 
                                label_annotation=NULL,border=FALSE,
                                clustering_method="ward",
                                show_colnames = FALSE,
                                show_rownames = TRUE,
                                fontsize_row=2)
rm(pheatmap.multilevel.splsda2fact); detach(package:pheatmap)

pheatmap.multilevel(res.2level_mixomics, 
                    col_sample=col.sample, 
                    col_stimulation=col.stimu, 
                    col_time=col.time,
                    label_color_stimulation=label.stimu,
                    label_color_time=label.time, 
                    label_annotation=NULL,border=FALSE,
                    clustering_method="ward",
                    show_colnames = FALSE,
                    show_rownames = TRUE,
                    fontsize_row=2)

# -------------
# other outputs need to be tested, e.g. plotVar
# plot3D, pltVar3D
# ------------





# ==================
## Third example: one-factor integrative analysis with sPLS
# =================
  rm(list=ls())
  source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')
  library(mixOmics)

  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                    6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9, 
                    10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14, 
                    13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
  design <- data.frame(repeat.indiv, liver.toxicity$treatment$Dose.Group)
  
  result.rat_IG <- multilevel(X, Y, design = design, ncomp = 2, keepX = c(50, 50), 
                               keepY = c(5, 5), method = 'spls', mode = "canonical") 
  
  rm(multilevel)
  result.rat_mixomics <- multilevel(X = liver.toxicity$gene, Y=liver.toxicity$clinic, 
                                    cond = liver.toxicity$treatment$Dose.Group, 
                                    sample = repeat.indiv, ncomp = 2, keepX = c(50, 50), 
                                    keepY = c(5, 5), method = 'spls', mode = 'canonical') 

  names(result.rat_IG)
  names(result.rat_mixomics)
  
  for (i in names(result.rat_mixomics)[!names(result.rat_mixomics) %in% c("call", "sample", "name.condition")]){
    if (!all.equal(result.rat_IG[[i]], result.rat_mixomics[[i]]) == TRUE) {
      stop(i)
    }
  }
  ### Dicrepancy found between IG update and mixomics
  # call: ok 
  # sample: ok input / output present only in multilevel (mixOmics)
  # name.condition: output missing in multilevel (update_IG) ??? should be added ???

  # variable plots
  plotVar(result.rat_IG, comp = 1:2, X.label = TRUE, Y.label = TRUE, cex = c(0.5, 0.9)) 
  plotVar(result.rat_mixomics, comp = 1:2, X.label = TRUE, Y.label = TRUE, cex = c(0.5, 0.9)) 
  
  cim_IG = cim(result.rat_IG, comp = 1:2, xlab = "genes", ylab = "clinic var",  margins = c(5, 6))
  cim_mixomics = cim(result.rat_mixomics, comp = 1:2, xlab = "genes", ylab = "clinic var",  margins = c(5, 6))
  
  all.equal(cim_IG, cim_mixomics)
  # cim_IG = cim_mixomics

  network(result.rat_IG, comp = 1:2, threshold = 0.8, 
          Y.names = NULL, keep.var = TRUE,
          color.node = c( "lightcyan","mistyrose"),
          shape.node = c("circle", "rectangle"),
          color.edge = c("red", "green"),
          lty.edge = c("solid", "solid"))
  
  network(result.rat_mixomics, comp = 1:2, threshold = 0.8, 
          Y.names = NULL, keep.var = TRUE,
          color.node = c( "lightcyan","mistyrose"),
          shape.node = c("circle", "rectangle"),
          color.edge = c("red", "green"),
          lty.edge = c("solid", "solid"))

# ==================== TO DO =============================

# why are the functions Split.variation.one.level  and Split.variation.two.level not used any more?
# these functions were used for other purposes (multilevel PCA for example, see paper non submitted.)


# we need to recall it in the multilevel functions, otherwise if we change 1 thing in 1 function we have to change it 2 times.
# we also need to add a detailed help file on these two functions which can be used external to multilevel.R. And add all the warnings


# help file have been added for the split.variation.one.level and the split.variation.two.level

## this below works but is now obsolete since I coded a generic split.variation function (see below)

# # ------------
# ## testing the split.variation.one.level function which has been RE-INCLUDED in the code and improved.
# # ------------
# source('multilevel/multilevel.R')
# 
# # function has been renamed 'split' instead of Split
# 
# ## First example: one-factor analysis with sPLS-DA
# data(vac18.study) # vac18 study
# 
# X <- vac18.study$genes
# dim(X)
# Xw.res = split.variation.one.level(X, rep.measures =  vac18.study$sample)
# dim(Xw.res)
# 
# rownames(Xw.res)
# 
# 
# # ------------
# ## testing the split.variation.two.level function which has been RE-INCLUDED in the code and improved.
# # ------------
# source('multilevel/multilevel.R')
# 
# 
# # function has been renamed 'split' instead of Split
# 
# 
# data(vac18.simulated) # simulated data
# 
# X <- vac18.simulated$genes
# factors <- data.frame(
#   stimul = vac18.simulated$stimulation,
#   time = vac18.simulated$time
#   )
# 
# Xw.res2 = split.variation.two.level(X, factors = factors, rep.measures =  vac18.simulated$sample)
# dim(Xw.res2)


# ================
# I created a generic split.variation function that merges the two
# ================

data(vac18.simulated) # simulated data

source('multilevel/multilevel.R')


## example for one level
X <- vac18.simulated$genes
Xw.res1level = split.variation(X, rep.measures =  vac18.simulated$sample)
dim(Xw.res1level)

# testing few warnings
Xw.res1level = split.variation(X, factors = vac18.simulated$stimulation, rep.measures =  vac18.simulated$sample)
Xw.res1level = split.variation(X, rep.measures =  vac18.simulated$sample[1:10])
Xw.res1level = split.variation(X, rep.measures =  as.character(vac18.simulated$sample))



## example for two level
factors <- data.frame(
  stimul = vac18.simulated$stimulation,
  time = vac18.simulated$time
)
Xw.res2level = split.variation(X, factors =factors, rep.measures =  vac18.simulated$sample)
dim(Xw.res2level)


# ===================== TO DO ============================
# test the tune.multilevel

## First example: one-factor analysis with sPLS-DA
rm(list=ls())
source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')
library(mixOmics)

data(data.simu) # simulated data
design <- data.frame(data.simu$sample, data.simu$stimu)

source('multilevel_new/tune.multilevel.R'); source('multilevel_new/tune.splsdalevel1.R')

# *BG* Remark: could be compared only in loo validation => difference expected in Mfold

result.ex1_IG = tune.multilevel(data.simu$X,
                               design = design,
                               ncomp=2,
                               test.keepX=c(5, 10, 15), 
                               already.tested.X = c(50),
                               method = 'splsda',
                               dist = 'mahalanobis.dist',
                               validation = 'loo') 

result.ex1_IG

rm("tune.multilevel"); rm("tune.splsdalevel1")
result.ex1_mixomics = tune.multilevel(data.simu$X,
                                       cond = data.simu$stimu,
                                       sample = data.simu$sample, 
                                       ncomp=2,
                                       test.keepX=c(5, 10, 15), 
                                       already.tested.X = c(50),
                                       method = 'splsda',
                                       dist = 'mahalanobis.dist',
                                       validation = 'loo') 

result.ex1_mixomics

all.equal(result.ex1_mixomics, result.ex1_IG)

## Second example: two-factor analysis with sPLS-DA
### */* Example 1 */* ###
rm(list=ls())
source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')
library(mixOmics)

data(liver.toxicity)
dose = liver.toxicity$treatment$Dose.Group
time = liver.toxicity$treatment$Time.Group
dose.time = cbind(dose, time)
repeat.indiv = c(1, 2, 1,  2,  1,  2,  1,  2,  3,  3,  4,  
                 3,  4,  3,  4,  4,  5,  6,  5,  5,  6,  5,  6,  7,  7,  
                 8,  6,  7,  8,  7,  8,  8,  9, 10,  9, 10, 11, 9,  9, 
                 10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 
                 14, 13, 14, 15, 16, 15, 16, 15, 16, 15, 16)

design <- data.frame(sample = repeat.indiv, 
                     dose = dose,
                     time = time)

source('multilevel_new/tune.multilevel.R'); source('multilevel_new/tune.splsdalevel2.R')
result.ex2_IG = tune.multilevel(liver.toxicity$gene,
                                design = design, 
                                ncomp=2,
                                test.keepX=c(5, 10, 15), 
                                already.tested.X = c(50),
                                method = 'splsda',
                                dist = 'mahalanobis.dist') 
result.ex2_IG

rm("tune.multilevel"); rm("tune.splsdalevel2")
result.ex2_mixomics = tune.multilevel(liver.toxicity$gene,
                                      cond = dose.time,
                                      sample = repeat.indiv, 
                                      ncomp=2,
                                      test.keepX=c(5, 10, 15), 
                                      already.tested.X = c(50),
                                      method = 'splsda',
                                      dist = 'mahalanobis.dist') 
result.ex2_mixomics

all.equal(result.ex2_mixomics, result.ex2_IG)
# Remark: a slight difference is expected since the correlation is now calculate on centered/scaled data

### */* Example 2 */* ###
rm(list=ls())
source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')
require(mixOmics)

load('multilevel_new/vac18.simulated.rda')
attributes(vac18.simulated)
# data(vac18.simulated) # simulated data

X <- vac18.simulated$genes
design <- data.frame(sample = vac18.simulated$sample,
                     stimul = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

source('multilevel_new/tune.multilevel.R'); source('multilevel_new/tune.splsdalevel2.R')
result.ex2_IG = tune.multilevel(X,
                                design = design, 
                                ncomp=1,
                                test.keepX=c(1, 2, 3), 
                                already.tested.X = NULL,
                                method = 'splsda',
                                dist = 'mahalanobis.dist') 
result.ex2_IG

rm("tune.multilevel"); rm("tune.splsdalevel2")
result.ex2_mixomics = tune.multilevel(X,
                                      cond = cbind(vac18.simulated$stimulation, vac18.simulated$time),
                                      sample = vac18.simulated$sample, 
                                      ncomp=1,
                                      test.keepX=c(1, 2, 3), 
                                      already.tested.X = NULL,
                                      method = 'splsda',
                                      dist = 'mahalanobis.dist') 

result.ex2_mixomics

## Third example: one-factor integrative analysis with sPLS
rm(list=ls())
source('multilevel_new/multilevel.R'); source('multilevel_new/withinVariation.R')
library(mixOmics)

data(liver.toxicity)
dose = liver.toxicity$treatment$Dose.Group
repeat.indiv = c(1,2, 1,  2,  1,  2,  1,  2,  3,  3,  4,  
                 3,  4,  3,  4,  4,  5,  6,  5,  5,  6,  5,  6,  7,  7,  
                 8,  6,  7,  8,  7,  8,  8,  9, 10,  9, 10, 11, 9,  9, 
                 10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 
                 14, 13, 14, 15, 16, 15, 16, 15, 16, 15, 16)

design <- data.frame(sample = repeat.indiv, 
                     dose = dose)

source('multilevel_new/tune.multilevel.R'); source('multilevel_new/tune.splslevel.R')
result.ex3_IG = tune.multilevel(liver.toxicity$gene, liver.toxicity$clinic, 
                                design = design, 
                                ncomp=2,
                                test.keepX=c(5, 10, 15), 
                                test.keepY=c(2,3), 
                                already.tested.X = c(50), already.tested.Y = c(5),
                                method = 'spls') 

result.ex3_IG

rm("tune.multilevel"); rm("tune.splslevel")
result.ex3_mixomics = tune.multilevel(liver.toxicity$gene, liver.toxicity$clinic, 
                                      cond = dose,
                                      sample = repeat.indiv, 
                                      ncomp=2,
                                      test.keepX=c(5, 10, 15), 
                                      test.keepY=c(2,3), 
                                      already.tested.X = c(50), already.tested.Y = c(5),
                                      method = 'spls') 

result.ex3_mixomics

result.ex3_IG$cor.value - result.ex3_mixomics$cor.value

# *BG* remark: difference expected => better if we scale the data

result.ex3_mixomics = tune.multilevel(scale(liver.toxicity$gene), scale(liver.toxicity$clinic), 
                                      cond = dose,
                                      sample = repeat.indiv, 
                                      ncomp=2,
                                      test.keepX=c(5, 10, 15), 
                                      test.keepY=c(2,3), 
                                      already.tested.X = c(50), already.tested.Y = c(5),
                                      method = 'spls') 

result.ex3_IG$cor.value - result.ex3_mixomics$cor.value
# ====================== TO DO ============================

# update the website
# ==================================================


# run it on terminal
R
library(mixOmics, lib='MyR')
