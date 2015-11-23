require(lattice)
require(ellipse)
require(ggplot2)

library(mixOmics)
library(R.utils)

sourceDirectory(path='C:/Users/Yan/Downloads/mixOmics/')


data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

liver.spca <- spca(X, ncomp = 1, keepX = 10)
selectVar(liver.spca)
selectVar(liver.spca,comp=c(1,2))
selectVar(liver.spca,comp=2)



liver.spls = spls(X, Y, ncomp = 2, keepX = c(20, 40),keepY = c(5, 5))
selectVar(liver.spls, comp = 2)

data(srbct)
X = srbct$gene
Y = srbct$class

srbct.splsda = splsda(X, Y, ncomp = 2, keepX = c(5, 10))
select = selectVar(srbct.splsda, comp = 2)
select
srbct$gene.name[substr(select$select, 2,5),]  


data(nutrimouse)

diet = unmap(nutrimouse$diet)
blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = c(3, 2, 2))
selectVar(nutri.sgcca) ###
selectVar(nutri.sgcca,block=c(1,2,3))




