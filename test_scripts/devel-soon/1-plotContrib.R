# -----------------------------------------------------------------------------------
# Testing-contrib.R
# Author:    KA Le Cao 
# Date started:  27/07/2015
# Last updated:  
# Objective:
# Latest update: 
# -----------------------------------------------------------------------------------

rm(list=ls())
library(mixOmics)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# Florian
##sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R

# KA
sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R


# ======================
# splsda
# =====================

# the object
# ---------
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))
# ---------


# contribution on comp 1, based on the median. Colors indicate the group in which the median expression is maximal
plotContrib(splsda.liver, comp = 1, method = 'median')
# contribution on comp 2, based on median. Colors indicate the group in which the median expression is maximal
plotContrib(splsda.liver, comp = 2, method = 'median')
# contribution on comp 2, based on median. Colors indicate the group in which the median expression is minimal
plotContrib(splsda.liver, comp = 2, method = 'median', contrib = 'min')

# testing that above ncomp it breaks: ok ERROR output
plotContrib(splsda.liver, comp = 3, method = 'median')

# changing the name to gene names
?liver.toxicity

# if the user input a name.var but names(name.var) is NULL, then a warning will be output and assign names of name.var to colnames(X)
# this is to make sure we can match the name of the selected variables to the contribution plot.
name.var = liver.toxicity$gene.ID[, 'geneBank']
length(name.var)
plotContrib(splsda.liver, comp = 2, method = 'median', name.var = name.var)

# if names are provided: ok, even when NAs
name.var = liver.toxicity$gene.ID[, 'geneBank']
names(name.var) = rownames(liver.toxicity$gene.ID)
plotContrib(splsda.liver, comp = 2, method = 'median', name.var = name.var)

# reduce the cex size of the names
name.var = liver.toxicity$gene.ID[, 'gene.title']
names(name.var) = rownames(liver.toxicity$gene.ID)
plotContrib(splsda.liver, comp = 2, method = 'median', name.var = name.var, cex.name = 0.5)

# if we want to have a look at the contribution (median) for each variable, need to store in an object
# (the plot will still be output)
output.contrib = plotContrib(splsda.liver, comp = 2, method = 'median', name.var = name.var)
output.contrib$contrib

# if we want to change the title of the legend
plotContrib(splsda.liver, comp = 2, method = 'median',legend.title = 'Time')


# we can change the title name
plotContrib(splsda.liver, comp = 2, method = 'median',legend.title = 'Time', title = 'Contribution plot')

# if we dont want to show the legend
plotContrib(splsda.liver, comp = 2, method = 'median', legend = FALSE)

# if we want to change the color of the legend
plotContrib(splsda.liver, comp = 2, method = 'median', legend.color = c(1:4))


# if we are only interested in the numerical outputs and no plot
plot.contrib = plotContrib(splsda.liver, comp = 2, method = 'median', plot = FALSE)
head(plot.contrib$contrib)



# splsda multilevel
# -----------------
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)

# multilevel sPLS-DA model
res.1level <- multilevel(X, ncomp = 3, design = design,
                         method = "splsda", keepX = c(30, 137, 123))


name.var = vac18$tab.prob.gene[, 'Gene']
names(name.var) = colnames(X)
# we can change the cex
plotContrib(res.1level, comp = 1, method = 'median', legend.title = 'Stimu', name.var = name.var, cex.name = 0.5)
plotContrib(res.1level, comp = 2, method = 'median', legend.title = 'Stimu', name.var = name.var, cex.name = 0.2)

# too many transcripts? only output the top ones
plotContrib(res.1level, comp = 2, method = 'median', legend.title = 'Stimu', name.var = name.var, cex.name = 0.5, ndisplay = 60)

# a simple PLSDA
# ----------------
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)

selectVar(plsda.breast, comp = 1)$name

# warning appears and no output
plotContrib(plsda.breast, comp = 1, method = 'median', ndisplay = 60, legend.color = color.mixo(2))




