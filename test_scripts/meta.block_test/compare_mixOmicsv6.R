###########################################################################################################################################################
##                  comparison of mixOmicsv6 and mixOmicsv5
##                      expect different results because of different algorithm etc
###########################################################################################################################################################
rm(list=ls())
setwd("/Users/florian/Work/git/package-mixOmics/")
library(mixOmics)


source("mixOmics/R/check_entry.R")
source("mixOmics/R/helpers.R")
source("mixOmics/R/sparse.meta.block.R")

source("mixOmics/R/meta.spls.hybrid.R")
source("mixOmics/R/pls.R")
source("mixOmics/R/plsda.R")
source("mixOmics/R/spls.R")
source("mixOmics/R/splsda.R")
source("mixOmics/R/meta.pls.R")
source("mixOmics/R/meta.plsda.R")
source("mixOmics/R/meta.spls.R")
source("mixOmics/R/meta.splsda.R")

source("mixOmics/R/wrapper.sparse.meta.block.R")
source("mixOmics/R/block.pls.R")
source("mixOmics/R/block.spls.R")
source("mixOmics/R/block.plsda.R")
source("mixOmics/R/block.splsda.R")
source("mixOmics/R/meta.block.pls.R")
source("mixOmics/R/meta.block.spls.R")
source("mixOmics/R/meta.block.plsda.R")
source("mixOmics/R/meta.block.splsda.R")


source("mixOmics/R/wrapper.rgcca.R")
source("mixOmics/R/wrapper.sparse.rgcca.R")

source("mixOmics/R/mixOmics.R")



#
abslist <- function(list = NULL, var = NULL) {
for (i in var){
    if (is.list(list[[which(names(list) == i)]])){
        list[[which(names(list) == i)]] = lapply(list[[which(names(list) == i)]], abs)
    } else {
        list[[which(names(list) == i)]] = abs(list[[which(names(list) == i)]])
    }
}
return(list)
}


#=============================================================================================================================================
#                           PLS
#=============================================================================================================================================

data(nutrimouse)
X <- scale(nutrimouse$lipid)
Y <- scale(nutrimouse$gene)

plsv5=pls(X,Y)
plsv6=wrapper.pls(X,Y)


plsv5.abs = abslist(plsv5, c("variates", "loadings"))
plsv6.abs = abslist(plsv6, c("variates", "loadings"))


checklist = names(plsv5)
checklist = checklist[!checklist %in% c("mat.c","mat.d","mat.e","max.iter","iter","nzv")]

lapply(checklist, function(x){all.equal(plsv5.abs[names(plsv5.abs) == x],
    plsv6.abs[names(plsv6.abs) == x])})


#=============================================================================================================================================
#                           meta.PLS
#=============================================================================================================================================


source("../multi-group/meta.pls.R")
source("../multi-group/meta.spls.R")
source("../multi-group/helpers.R")

study=factor(c(rep(1,15),rep(2,15),rep(3,10)))

plsv5=meta.pls(X,Y,study=study)



source("mixOmics/R/check_entry.R")
source("mixOmics/R/helpers.R")
plsv6=wrapper.meta.pls(X,Y,study=study)


all.equal(plsv5$loadings.global$X,plsv6$loadings$X)
all.equal(plsv5$variates.global$X,plsv6$variates$X)
all.equal(plsv5$loadings.partial$X,plsv6$loadings.partial$X)




#=============================================================================================================================================
#                           PLSDA
#=============================================================================================================================================

load("test_scripts/meta.block_test/Fibro-ESC-iPSC.6exp.167samples.light.Rdata") #load data, type.id, exp among others

#one study
ind=which(exp=="Briggs")
study.light=exp[ind]
type.id.light=type.id[ind]
data.light=data[ind,]
Y.mat.light=unmap(type.id.light)
rownames(Y.mat.light)=rownames(data.light)

plsdav5=plsda(X=data.light,Y=type.id.light,ncomp=3,near.zero.var=FALSE)
plsdav6=wrapper.plsda(X=data.light,Y=type.id.light,ncomp=3,near.zero.var=FALSE)

plsdav5.abs = abslist(plsdav5, c("variates", "loadings"))
plsdav6.abs = abslist(plsdav6, c("variates", "loadings"))


checklist = names(plsdav5)
checklist = checklist[!checklist %in% c("mat.c","mat.d","mat.e","max.iter","iter","nzv")]

lapply(checklist, function(x){all.equal(plsdav5.abs[names(plsdav5.abs) == x],
    plsdav6.abs[names(plsdav6.abs) == x])})
#differences between Y, ind.mat. By design Y is now the factor, to be used in plotIndiv


#=============================================================================================================================================
#                           meta.PLSDA
#=============================================================================================================================================

source("../multi-group/meta.pls.R")
source("../multi-group/meta.spls.R")
source("../multi-group/helpers.R")

Y.mat=unmap(type.id)
rownames(Y.mat)=rownames(data)
plsdav5=meta.pls(data,Y.mat,study=exp)

source("mixOmics/R/check_entry.R")
source("mixOmics/R/helpers.R")
plsdav6=wrapper.meta.plsda(X=data,Y=type.id,ncomp=2,study=exp,near.zero.var=FALSE)


all.equal(plsdav5$loadings.global$X,plsdav6$loadings$X)
all.equal(plsdav5$variates.global$X,plsdav6$variates$X)
all.equal(plsdav5$loadings.partial$X,plsdav6$loadings.partial$X)
all.equal(plsdav5$variates.partial$X,plsdav6$variates.partial$X)




#=============================================================================================================================================
#                           sPLS
#=============================================================================================================================================

data(nutrimouse)
X <- scale(nutrimouse$lipid)
Y <- scale(nutrimouse$gene)

source("mixOmics/R/spls.R")
splsv5=spls(X,Y,ncomp=3,keepX=c(10,10,10))
splsv6=wrapper.spls(X,Y,ncomp=3,keepX=c(10,10,10))


splsv5.abs = abslist(splsv5, c("variates", "loadings"))
splsv6.abs = abslist(splsv6, c("variates", "loadings"))


checklist = names(splsv5)
checklist = checklist[!checklist %in% c("mat.c","mat.d","mat.e","max.iter","iter","nzv")]

lapply(checklist, function(x){all.equal(splsv5.abs[names(splsv5.abs) == x],
    splsv6.abs[names(splsv6.abs) == x])})



#=============================================================================================================================================
#                           meta.sPLS
#=============================================================================================================================================


source("../multi-group/meta.pls.R")
source("../multi-group/meta.spls.R")
source("../multi-group/helpers.R")
study=factor(c(rep(1,15),rep(2,15),rep(3,10)))

splsv5=meta.spls(X,Y,ncomp=3,keepX=c(10,10,10),study=study)

source("mixOmics/R/check_entry.R")
source("mixOmics/R/helpers.R")
source("mixOmics/R/meta.spls.R")
splsv6=wrapper.meta.spls(X,Y,ncomp=3,keepX=c(10,10,10),study=study)


all.equal(splsv5$loadings.global$X,splsv6$loadings$X)
all.equal(splsv5$variates.global$X,splsv6$variates$X)
all.equal(splsv5$loadings.partial$X,splsv6$loadings.partial$X)
all.equal(splsv5$variates.partial$X,splsv6$variates.partial$X)




#=============================================================================================================================================
#                           sPLSDA
#=============================================================================================================================================


#load("test_scripts/meta.block_test/Fibro-ESC-iPSC.6exp.167samples.light.Rdata") #load data, type.id, exp among others

#one study
ind=which(exp=="Briggs")
study.light=exp[ind]
type.id.light=type.id[ind]
data.light=data[ind,]
Y.mat.light=unmap(type.id.light)
rownames(Y.mat.light)=rownames(data.light)

splsdav5=splsda(X=data.light,Y=type.id.light,ncomp=3,near.zero.var=FALSE,keepX=c(10,10,10))
source("mixOmics/R/splsda.R")
splsdav6=wrapper.splsda(X=data.light,Y=type.id.light,ncomp=3,near.zero.var=FALSE,keepX=c(10,10,10))

splsdav5.abs = abslist(splsdav5, c("variates", "loadings"))
splsdav6.abs = abslist(splsdav6, c("variates", "loadings"))


checklist = names(splsdav5)
checklist = checklist[!checklist %in% c("mat.c","mat.d","mat.e","max.iter","iter","nzv")]

lapply(checklist, function(x){all.equal(splsdav5.abs[names(splsdav5.abs) == x],
    splsdav6.abs[names(splsdav6.abs) == x])})
#differences between Y, ind.mat. By design Y is now the factor, to be used in plotIndiv


#=============================================================================================================================================
#                           meta.sPLSDA
#=============================================================================================================================================


source("../multi-group/meta.pls.R")
source("../multi-group/meta.spls.R")
source("../multi-group/helpers.R")

Y.mat=unmap(type.id)
rownames(Y.mat)=rownames(data)
splsdav5=meta.spls(data,Y.mat,study=exp,ncomp=3,keepX=c(10,10,10))

source("mixOmics/R/check_entry.R")
source("mixOmics/R/helpers.R")
source("mixOmics/R/meta.splsda.R")
splsdav6=wrapper.meta.splsda(X=data,Y=type.id,ncomp=3,study=exp,near.zero.var=FALSE,keepX=c(10,10,10))


all.equal(splsdav5$loadings.global$X,splsdav6$loadings$X)
all.equal(splsdav5$variates.global$X,splsdav6$variates$X)
all.equal(splsdav5$loadings.partial$X,splsdav6$loadings.partial$X)
all.equal(splsdav5$variates.partial$X,splsdav6$variates.partial$X)




#=============================================================================================================================================
#                           rGCCA
#=============================================================================================================================================

#difference of algo and outputs

data(nutrimouse)
X <- scale(nutrimouse$lipid)
Y <- scale(nutrimouse$gene)
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
nutri.res$loadings$X

A = list(X = X, Y = Y)


rgccav5=wrapper.rgcca(data=A,tau=c(0.064, 0.008),ncomp=rep(2, length(A)))

source("mixOmics/R/rgcca.R")
rgccav6=rgcca.block(data=A,tau=c(0.064, 0.008),max.iter=1000,ncomp=rep(2, length(A)),tol=1e-26)
#rgccav6=rgcca.block(data=A,tau="optimal",max.iter=1000,ncomp=rep(2, length(A)))

plot(rgccav5$variates[[1]][,1],rgccav6$variates[[1]][,1])

all.equal(rgccav5,rgccav6)





v6=rgcca.block(data=A,max.iter=1000,ncomp=rep(2, length(A)),tol=1e-26) #tau=1, is it supposed to be equal to something else?



