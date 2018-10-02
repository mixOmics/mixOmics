# ------------------------
# date: 25/04/2015
# author: Kim-Anh Le Cao
# including the tau.estim (from RGCCA/SHIP) in mixOmics to enable regularisation parameters estimation and bypass tune.rcc
# ---------------------


library(mixOmics)
data(nutrimouse)
X = nutrimouse$gene
Y = nutrimouse$lipid

A = list(X, Y)

source('tau.estimate.R')
tau = sapply(A, tau.estimate)
# regularisation parameters for X (first) and then Y
tau

# note: scaling the data below is optional, in our case it seemed to give a good separation of the diets
rcc.result = rcc(scale(X), scale(Y), ncomp = 2, lambda1 = tau[1], lambda2 = tau[2])

plotIndiv(rcc.result, ind.names = nutrimouse$diet, col = as.numeric(nutrimouse$genotype))


## tuning the regularization parameters for 'rcc'
# ----------------------------------------------------
library(mixOmics)

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

# create a list as an input
A = list(X, Y)

source('../mixOmics/R/tau.estim.R')
param = sapply(A, tau.estim)
# regularisation parameters for X (first) and then Y
param

# note: we advise scaling the data as the rcc input in that case
rcc.result = rcc(scale(X), scale(Y), ncomp = 2, lambda1 = param[1], lambda2 = param[2])

# we observe quite a good deparation of the diets
plotIndiv(rcc.result, ind.names = nutrimouse$diet, col = color.mixo(as.numeric(nutrimouse$genotype)))


# this is if we dont do any scaling:
# our feeling (Benoit and I) is thatwhether you scale or not should not affect the rcc result
# however the data are centered and scaled in tau.estim to estimate the crossprod, therefore the regularization estimates are 
# dependent on that scaling.

rcc.result = rcc(X, Y, ncomp = 2, lambda1 = param[1], lambda2 = param[2])

my.colors = color.mixo(as.numeric(nutrimouse$genotype))
plotIndiv(rcc.result, ind.names = nutrimouse$diet, col = my.colors)


