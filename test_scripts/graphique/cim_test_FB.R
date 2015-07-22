library(mixOmics)
library(R.utils)



set.seed(123)
mat <- matrix(rnorm(200), 20, 10)
mat[ 1:10, seq(1, 10, 2)] <- mat[ 1:10, seq(1, 10, 2)] - 4
mat[11:20, seq(2, 10, 2)] <- mat[11:20, seq(2, 10, 2)] + 2
mat[15:20, seq(2, 10, 2)] <- mat[15:20, seq(2, 10, 2)] + 3
mat <- mat[sample(20), ]
colnames(mat) <- paste("sample", 1:10, sep = " ")
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

#test legend ####
cim(mat, col.sideColors = cond.col[cond],
    legend=list( legend = names(cond.col), col = cond.col,title = "Condition", cex = 0.7))
cim(mat, col.sideColors = cond.col[cond],
    legend=list( legend = names(cond.col), col = cond.col,title = "Condition", cex = 1.5))
cim(mat, col.sideColors = cond.col[cond],row.sideColors = gene.col[cl],
    legend=list( legend = names(gene.col), col = gene.col,title = "Condition", cex = 0.7))


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
cim(mat,main="plot",xlab="sample",ylab="gene")

#test margins
cim(mat,main="plot",xlab="sample",ylab="gene",margins=c(7,7))
cim(mat,main="Cim.test",xlab="sample",ylab="gene",margins=c(5.5,5))





data(liver.toxicity)
X <- liver.toxicity$gene

liver.spca <- spca(X, ncomp = 3, keepX = rep(30, 3), scale = FALSE)

cim(liver.spca)

#test sample.names
cim(liver.spca,sample.names = liver.toxicity$treatment[, 3])
cim(liver.spca,sample.names=rep(c("alpha","beta"),each=32))

#test sample.sideColors
cim(liver.spca,sample.names = liver.toxicity$treatment[, 3],sample.sideColors=color.jet(4)[factor(liver.toxicity$treatment[, 3])],
    legend = list( legend = levels(factor(liver.toxicity$treatment[, 3])), col = color.jet(4)))
cim(liver.spca,sample.sideColors=rep(c("blue","green"),each=32))


#test clust.method
cim(liver.spca, clust.method = c("single", "mcquitty"))
cim(liver.spca, clust.method = c("average", "centroid"))
cim(liver.spca, clust.method = c("median", "complete"))



#test var.names, sample.names
cim(liver.spca,var.names=FALSE)
cim(liver.spca,sample.names=FALSE)
cim(liver.spca,var.names=FALSE,sample.names=FALSE)



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




#test x.sideColors
col=palette(rainbow(21))
cim(nutri.rcc, mapping = "X", x.sideColors = col)
cim(nutri.rcc, mapping = "XY", x.sideColors = col,zoom=T)
cim(nutri.rcc, mapping = "Y", x.sideColors = col)

#test y.sideColors
col=palette(rainbow(120))
cim(nutri.rcc, mapping = "X", y.sideColors = col)
cim(nutri.rcc, mapping = "XY", y.sideColors = col)
cim(nutri.rcc, mapping = "Y", y.sideColors = col)



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
    var.names = FALSE,
    legend=list(legend=unique(design$stim),col=unique(stim.col),title="Stimulations",cex=0.9))




data(liver.toxicity)
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)

design <- data.frame(sample = repeat.indiv) 
res.spls.1level <- multilevel(X = liver.toxicity$gene,
                              Y=liver.toxicity$clinic,
                              design = design,
                              ncomp = 3,
                              keepX = c(50, 50, 50), keepY = c(5, 5, 5),
                              method = 'spls', mode = 'canonical')

stim.col <- c("darkblue", "purple", "green4","red3")
cim(res.spls.1level,mapping="Y",
    sample.sideColors = stim.col[factor(liver.toxicity$treatment[,3])],
    legend=list(legend=unique(liver.toxicity$treatment[,3]),col=stim.col,cex=0.9))