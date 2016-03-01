<<<<<<< HEAD
library(mixOmics)

## default method
#------------------------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

cim(cor(X, Y), cluster = "none")


## CIM representation for objects of class 'rcc'
#------------------------------------------------------------------
nutri.rcc <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

cim(nutri.rcc, xlab = "genes", ylab = "lipids", margins = c(5, 6))

#-- interactive 'zoom' 
\dontrun{
  
  cim(nutri.rcc, xlab = "genes", ylab = "lipids", margins = c(5, 6), 
      zoom = TRUE)
  #-- select the region and "see" the zoom-out region
}

#-- cim from X matrix
diet.col <- palette()[as.numeric(nutrimouse$diet)]
cim(nutri.rcc, mapping = "X", sample.names = nutrimouse$diet,
    sample.sideColors = diet.col, xlab = "lipids",
    clust.method = c("ward", "ward"), margins = c(6, 4))

#-- cim from Y matrix
geno.col = c("orange", "magenta")[as.numeric(nutrimouse$genotype)]
cim(nutri.rcc, mapping = "Y", sample.names = nutrimouse$genotype,
    sample.sideColors = geno.col, xlab = "genes",
    clust.method = c("ward", "ward"))


## CIM representation for objects of class '(s)(i)pca'
#------------------------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene

spca.res <- spca(X, ncomp = 2, keepX = c(30, 30), scale = FALSE)

dose.col <- palette()[as.numeric(as.factor(liver.toxicity$treatment[, 3]))]

cim(spca.res, sample.sideColors = dose.col, var.names = FALSE,
    sample.names = liver.toxicity$treatment[, 3],
    clust.method = c("ward", "ward"))


## CIM representation for objects of class '(s)pls' 
#------------------------------------------------------------------
Y <- liver.toxicity$clinic

toxicity.spls <- spls(X, Y, ncomp = 2, keepX = c(40, 30))
cim(toxicity.spls, transpose = TRUE,   
    clust.method = c("ward", "ward"), margins = c(5, 7))

#-- cim from X matrix
cim(toxicity.spls, mapping = "X", sample.sideColors = dose.col, 
    sample.names = liver.toxicity$treatment[, 3])


## CIM representation for objects of class '(s)plsda' 
#------------------------------------------------------------------
Y <- liver.toxicity$treatment[, 3]

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(40, 30))

cim(splsda.liver, sample.sideColors = dose.col, sample.names = Y)


## CIM representation for objects of class 'multilevel' 
#------------------------------------------------------------------
data(vac18.simulated)
X <- vac18.simulated$genes
design <- data.frame(samp = vac18.simulated$sample,
                     stim = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

res.2level <- multilevel(X, ncomp = 2, design = design,
                         keepX = c(120, 10), method = 'splsda')

stim.col <- c("darkblue", "purple", "green4","red3")
stim.col <- stim.col[as.numeric(design$stim)]
time.col <- c("orange", "cyan")[as.numeric(design$time)]

cim(res.2level, sample.sideColors = cbind(stim.col, time.col), 
    sample.names = paste(design$time, design$stim, sep = "_"),
    var.names = FALSE)
=======
library(mixOmics)
library(R.utils)

sourceDirectory(path='C:/Users/Yan/Downloads/mixOmics1')



set.seed(123)
mat <- matrix(rnorm(200), 20, 10)
mat[ 1:10, seq(1, 10, 2)] <- mat[ 1:10, seq(1, 10, 2)] - 4
mat[11:20, seq(2, 10, 2)] <- mat[11:20, seq(2, 10, 2)] + 2
mat[15:20, seq(2, 10, 2)] <- mat[15:20, seq(2, 10, 2)] + 3
mat <- mat[row(20), ]
colnames(mat) <- paste("row", 1:10, sep = " ")
rownames(mat) <- paste("gene", 1:20, sep = "")


obj.cim=cim(mat)

#test color
col=palette(rainbow(51))

cim(mat, color = color.GreenRed(51))
cim(mat, color = col)

#test col.sideColors
cond=rep(c("cond 1","cond 2"),5)
cond.col <- c("cond 1" = "darkviolet", "cond 2" = "darkorange")

cim(mat, col.sideColors = 1:10)
cim(mat, col.sideColors = cond.col[cond])

#test row.sideColors
ddr <- as.hclust(obj.cim$ddr)
cl <- cutree(ddr, k = 3)
gene.col <- c("up" = "red", "down" = "green", "null" = "black")

cim(mat, row.sideColors = 1:20)
cim(mat, row.sideColors = gene.col[cl])


#test legend 
cim(mat, col.sideColors = cond.col[cond],
    legend=list( legend = names(cond.col), col = cond.col,title = "Condition", cex = 0.7))
cim(mat, col.sideColors = cond.col[cond],
    legend=list( legend = names(cond.col), col = cond.col,title = "Condition", cex = 1.5))
cim(mat, col.sideColors = cond.col[cond],row.sideColors = gene.col[cl],
    legend=list( legend = names(gene.col), col = gene.col,title = "Condition", cex = 0.7))
cim(mat, col.sideColors = cond.col[cond],row.sideColors = gene.col[cl],
    legend=list( legend = names(gene.col), col = gene.col))


#test zoom
cim(mat,zoom=TRUE)  


#test cluster
cim(mat,cluster="column")
cim(mat,cluster="none") 
cim(mat,cluster="row")


#test transpose
cim(mat,transpose=TRUE)
cim(mat,cluster="column",transpose=TRUE)
cim(mat,cluster="none",transpose=TRUE) 
cim(mat,cluster="row",transpose=TRUE)


#test row.cex, col.cex
cim(mat,row.cex=0.8,col.cex=1.2)
cim(mat,row.cex=1.2,col.cex=0.8)
cim(mat,row.cex=1.2,col.cex=0.8, transpose=TRUE)


#test dist.method
cim(mat,dist.method=c("correlation","maximum"))
cim(mat,dist.method=c("manhattan","canberra"))
cim(mat,dist.method=c("binary","minkowski"))

#test cut.tree
cim(mat,cut.tree=c(0.5,0.2))
cim(mat,cut.tree=c(0.99,0))

#test simkey, keysize
cim(mat,symkey=FALSE)
cim(mat,color = color.GreenRed(51), symkey=FALSE)
cim(mat,keysize=c(5,1))
cim(mat,keysize=c(1,5),symkey=FALSE)

#test lhei, lwid
cim(mat,lhei=c(3,4),lwid=c(2,5))
cim(mat,lhei=c(1,2))
cim(mat,lwid=c(1,2))


#test main, xlab, ylab
cim(mat,main="plot",xlab="row",ylab="gene")

#test margins
cim(mat,main="plot",xlab="row",ylab="gene",margins=c(7,7))
cim(mat,main="Cim.test",xlab="row",ylab="gene",margins=c(5.5,5))





data(liver.toxicity)
X <- liver.toxicity$gene

liver.spca <- spca(X, ncomp = 3, keepX = rep(30, 3), scale = FALSE)

cim(liver.spca)

#test row.names
cim(liver.spca,row.names = liver.toxicity$treatment[, 3])
cim(liver.spca,row.names=rep(c("alpha","beta"),each=32))

#test row.sideColors
cim(liver.spca,row.names = liver.toxicity$treatment[, 3],row.sideColors=color.jet(4)[factor(liver.toxicity$treatment[, 3])],
    legend = list( legend = levels(factor(liver.toxicity$treatment[, 3])), col = color.jet(4)))
cim(liver.spca,row.sideColors=rep(c("blue","green"),each=32))


#test clust.method
cim(liver.spca, clust.method = c("single", "mcquitty"))
cim(liver.spca, clust.method = c("average", "centroid"))
cim(liver.spca, clust.method = c("median", "complete"))

#test cutoff
cim(liver.spca,cutoff=0.5)
cim(liver.spca,cutoff=0.8)
cim(liver.spca,cutoff=0.8,cluster="column")
cim(liver.spca,cutoff=0.8,cluster="row")


#test col.names, row.names
cim(liver.spca,col.names=FALSE)
cim(liver.spca,row.names=FALSE)
cim(liver.spca,col.names=FALSE,row.names=FALSE)



data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
diet <- nutrimouse$diet
geno <- nutrimouse$genotype

nutri.rcc <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

cim(nutri.rcc)


#test mapping
cim(nutri.rcc,mapping="X")
cim(nutri.rcc,mapping="Y")
cim(nutri.rcc,mapping="XY")

#test comp
cim(nutri.rcc,comp=1)
cim(nutri.rcc,comp=3)
cim(nutri.rcc,comp=c(1,2))
cim(nutri.rcc,comp=c(2,1))
cim(nutri.rcc,comp=c(1,2,2))



#test row.sideColors
col=palette(rainbow(21))
cim(nutri.rcc, mapping = "X", col.sideColors = col)
cim(nutri.rcc, mapping = "XY", row.sideColors = col,zoom=T)
cim(nutri.rcc, mapping = "Y", row.sideColors = color.jet(40))

#test col.sideColors
col=color.jet(120)
cim(nutri.rcc, mapping = "X", row.sideColors = color.jet(40))
cim(nutri.rcc, mapping = "XY", col.sideColors = col)
cim(nutri.rcc, mapping = "Y", col.sideColors = col)

#test cutoff
cim(nutri.rcc,cutoff=0.5, mapping="X")
cim(nutri.rcc,cutoff=0.8, mapping="X")
cim(nutri.rcc,cutoff=0.6,cluster="column", mapping="X")
cim(nutri.rcc,cutoff=0.6,cluster="row", mapping="X")
cim(nutri.rcc,cutoff=0.5, mapping="Y")
cim(nutri.rcc,cutoff=0.8, mapping="Y")
cim(nutri.rcc,cutoff=0.6,cluster="column", mapping="Y")
cim(nutri.rcc,cutoff=0.6,cluster="row", mapping="Y")


data(vac18.simulated)
X <- vac18.simulated$genes
design <- data.frame(samp = vac18.simulated$sample,
                     stim = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

res.2level <- multilevel(X, ncomp = 2, design = design,
                         keepX = c(120, 10), method = 'splsda')

stim.col <- c("darkblue", "purple", "green4","red3")
stim.col <- stim.col[as.numeric(design$stim)]
time.col <- c("orange", "cyan")[as.numeric(design$time)]

cim(res.2level, row.sideColors = cbind(stim.col, time.col), 
    row.names = paste(design$time, design$stim, sep = "_"),
    col.names = FALSE,
    legend=list(legend=unique(design$stim),col=unique(stim.col),title="Stimulations"))

cim(res.2level, row.sideColors = cbind(stim.col, time.col), 
    row.names = paste(design$time, design$stim, sep = "_"),
    col.names = FALSE,
    legend=list(title="Stimulations"))



data(liver.toxicity)
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)

design <- data.frame(row = repeat.indiv) 
res.spls.1level <- multilevel(X = liver.toxicity$gene,
                              Y=liver.toxicity$clinic,
                              design = design,
                              ncomp = 3,
                              keepX = c(50, 50, 50), keepY = c(5, 5, 5),
                              method = 'spls', mode = 'canonical')

col <- c("darkblue", "purple", "green4","red3")
cim(res.spls.1level,mapping="Y",
    row.names=liver.toxicity$treatment[,3],
    row.sideColors = col[factor(liver.toxicity$treatment[,3])],
    legend=list(legend=unique(liver.toxicity$treatment[,3]),col=col))



data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

cim(toxicity.spls,mapping="X",comp=c(1,2,2,1))
cim(toxicity.spls,mapping="X",comp=c(1,2,2,3))
cim(toxicity.spls,mapping="Y",comp=c(1,2,3,1))
cim(toxicity.spls,mapping="Y",comp=1,row.names = 1:64)
cim(toxicity.spls,mapping="Y",comp=3, col.names = 1:10)


#test cutoff
cim(toxicity.spls,cutoff=0.5, mapping="X")
cim(toxicity.spls,cutoff=0.8, mapping="X")
cim(toxicity.spls,cutoff=0.6,cluster="column", mapping="X")
cim(toxicity.spls,cutoff=0.6,cluster="row", mapping="X")
cim(toxicity.spls,cutoff=0.5, mapping="Y")
cim(toxicity.spls,cutoff=0.8, mapping="Y")
cim(toxicity.spls,cutoff=0.6,cluster="column", mapping="Y")
cim(toxicity.spls,cutoff=0.6,cluster="row", mapping="Y")
>>>>>>> b52c17368bdc9075f68bea999239c3cddcd49838
