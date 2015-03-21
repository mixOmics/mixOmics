# ----------------------------------------
# testing after the NEW IG update from Jan 2015
# date: 15/02/2015
# latest update: 12/03/2015: testing CIM
# latest update: 13/03/2015 : reformatting for branch multilevel in bitbucket
# latest update: 14/03/2015: testing all examples for R package. compile is ok.
# latest update: 15/03/2015: columns 2 and 3 are no needed in withinVariation (ex for PCA multilevel) and for spls
#                            testing all examples for R package again. compile is ok.
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

# this is a spls (unsupervised analysis) so no need to mention any factor in design
# # we only perform a one level variation split
design <- data.frame(sample = repeat.indiv) 
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
# in design we only need to mention the repeated measurements to split the one level variation
design <- data.frame(sample = vac18$sample)
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
  # note: we made up those data, pretending they are repeated measurements
  repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                    6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                    10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                    13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
  summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
  
  # here we are only interested in a one level variation split since spls is an unsupervised method
  design <- data.frame(sample = repeat.indiv)
  
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

# ================================================================
# end of examples for Rd files testing 
# ================================================================

# R compile is now ok. (15/03/2015)


# ======== that's it for now, KA 15/03/2015 ======================


