library(mixOmics)
# ##library(R.utils)


rm(list=ls())

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}


# FB
##sourceDirectory(path='C:/Users/Yan/Downloads/mixOmics1')


# Florian
sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE) #load all the functions inside mixOmics/R

# KA
sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside mixOmics/R


# =======================================================
# simple CIM example on a matrix as input, by default distance = euclidian and mclust method = complete
# =======================================================

# Simulate some data:
# -----------------
set.seed(123)
mat <- matrix(rnorm(200), 20, 10)
mat[ 1:10, seq(1, 10, 2)] <- mat[ 1:10, seq(1, 10, 2)] - 4
mat[11:20, seq(2, 10, 2)] <- mat[11:20, seq(2, 10, 2)] + 2
mat[15:20, seq(2, 10, 2)] <- mat[15:20, seq(2, 10, 2)] + 3
mat <- mat[sample(20), ]
colnames(mat) <- paste("sample", 1:10, sep = " ")
rownames(mat) <- paste("gene", 1:20, sep = "")
dim(mat)


# The samples are grouped into 2 conditions, which we color accordingly
cond <- rep("cond 1", 10)
cond[seq(2, 10, 2)] <- "cond 2"
cond.col <- c("cond 1" = "darkviolet", "cond 2" = "darkorange")
cond.col

# this option shows horizontal sample bars and legend
obj.cim=cim(mat, color = color.GreenRed(51), col.sideColors = cond.col[cond],
            legend=list( legend = names(cond.col), col = cond.col, title = "Condition", cex = 0.9))


# Based on this cluster, we can extract the dendrogram ddr:
ddr <- as.hclust(obj.cim$ddr)
# and cut the trees to obtain 3 clusters
cl <- cutree(ddr, k = 3)
# we give a color to the genes, so that we can assign each gene to one of the 3 colors clusters, up, down and null:
gene.col <- c("up" = "red", "down" = "green", "null" = "black")
gene.col
# and now add a vertical side bar for the genes, colors are set to green red
cim(mat, color = color.GreenRed(51), col.sideColors = rep(c("red","blue"),each=5),
    row.sideColors=gene.col[cl],legend=list(legend = names(gene.col), col = gene.col,
                                            title = "Regulation", cex = 0.9))  

# Input matrix is transposed, legend is added, colors by default
cim(mat, transpose=TRUE,
    legend=list(legend = names(gene.col), col = gene.col,
                cex=0.8,title = "Regulation")) 

# row names and col names were changed, horizontal and vertical colored side bars added
cim(mat, 
    row.names=names(gene.col)[cl], row.sideColors=gene.col[cl],
    col.sideColors = cond.col[cond], col.names=names(cond.col)[rep(c(1,2),5)],
    legend=list(legend = names(gene.col), col = gene.col, title = "Regulation", cex = 0.8)) 


# data transposed, no side bar, no legend
cim(mat, row.cex=0.8, col.cex=1.2, transpose=TRUE) 

# data are transposed and the cluster is only performed on the columns
cim(mat,cluster="column", 
    legend=list(legend = names(gene.col), col = gene.col,
                title = "Regulation", cex = 0.9))

# data are transposed and the cluster is only performed on the columns
cim(mat,cluster="row") 


# distances are set to correlation distance for the rows, and max distance for the columns
cim(mat, dist.method=c("correlation","maximum"))


# cutting tree at different levels. That is pretty cool!
cim(mat,cut.tree=c(0.5,0.2))

# adding a title to the plot, label axes, adjusting margin
cim(mat, symkey=TRUE, keysize=c(1,1), main = "My plot", xlab="sample", ylab="gene", margins=c(7,7))

# Increasing the height of the dendrogram, results in an increase of the color key
cim(mat,lhei=c(3,4),lwid=c(2,5))

# Increasing the height of the dendrogram, results in an increase of the color key
# here the symkey is set to asymetric, see ??cim
cim(mat,lhei=c(3,4),lwid=c(2,5), symkey = FALSE)



# ========================================================
# CIM using mixOmics on a single data set (example with sPCA)
# ========================================================
data(liver.toxicity)
X <- liver.toxicity$gene

# 20 genes are selected on each component
liver.spca <- spca(X, ncomp = 3, keepX = rep(20, 3))

# Playing with the different available clustering methods, indicating the group colors of the samples based on acetaminophen dose
# note that because it is a sparse method, only the genes selected on the 3 components will appear
cim(liver.spca,
    sample.names = liver.toxicity$treatment[, 3], sample.sideColors=color.mixo(factor(liver.toxicity$treatment[, 3])),
    clust.method = c("single", "mcquitty"))

cim(liver.spca,
    sample.names = liver.toxicity$treatment[, 3], sample.sideColors=color.mixo(factor(liver.toxicity$treatment[, 3])),
    clust.method = c("average", "centroid"))

cim(liver.spca,
    sample.names = liver.toxicity$treatment[, 3],sample.sideColors=color.mixo(factor(liver.toxicity$treatment[, 3])),
    clust.method = c("median", "complete"))

# other names
cim(liver.spca, sample.names=rep(c("alpha","beta"),each=32))

# or no var.names, no sample.names
cim(liver.spca, var.names=FALSE, sample.names=FALSE)

# with a sparse method, show only the variables selected on specific components
# FB: le nom des variables (colonnes) est tjs le meme ce qui peut preter a confusion
# peut on extraire plutot le nom des variables?
##reparer erreur tres legere du ? object$names
# FB: peut etre que tu as oublie de livrer: mais ca ne change rien.

# ici les varaibles sont numerotees de 1:20
cim(liver.spca, comp = 1)
# ici les variables sont numerotees de 1:20, c est ca qui prete a confusion par rapport a au dessus. un utilisateur lambda croiait
# que ce sont les meme variables
cim(liver.spca, comp = 2)
cim(liver.spca, comp = c(1,3))

# ==================================================================
## CIM representation for objects of class 'sPLS-DA' 
#===================================================================
# Setting up the Y outcome first
Y <- liver.toxicity$treatment[, 3]

# 1 - sPLS-DA analysis, where 40 and 30 genes are selected on each component
splsda.liver1 <- splsda(X, Y, ncomp = 3, keepX = c(40, 30, 50))

# CIM default representation includes the total of 70 genes selected, with the dose color
dose.col <- color.mixo(as.numeric(as.factor(liver.toxicity$treatment[, 3])))
cim(splsda.liver1, sample.sideColors = dose.col, sample.names = Y)

# 2 - sPLS-DA analysis, where 40 genes are selected on the first component
splsda.liver2 <- splsda(X, Y, ncomp = 1, keepX = c(40))
# The default CIM will show only those 40 genes selected
cim(splsda.liver2, sample.sideColors = dose.col, sample.names = Y)

# with a sparse method, show only the variables selected on specific components
# FB: le nom des variables (colonnes) est tjs le meme ce qui peut preter a confusion
# peut on extraire plutot le nom des variables?
cim(splsda.liver1, comp = 1)
cim(splsda.liver1, comp = 2)
cim(splsda.liver1, comp = c(1,2))
cim(splsda.liver1, comp = c(1,2,2))
cim(splsda.liver1, comp = c(2,1))

# =============================================
## The CIM default method to visualise cross-correlations between two data sets 
# =============================================
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

# cross correlation between X and Y, by default this is Pearson correlation
cim(cor(X, Y), cluster = "none")

# cross correlation between X and Y, with Spearman correlation
cim(cor(X, Y, method = 'spearman'), cluster = "none")

# or one data set (just the matrix as an input)
cim(as.matrix(X), cluster = "none")


# other options: for X matrix
dose.col <- color.mixo(as.numeric(as.factor(liver.toxicity$treatment[, 3])))
cim(toxicity.spls, mapping = "X", sample.sideColors = dose.col, 
    sample.names = liver.toxicity$treatment[, 3])

# with a transpose and a Ward clustering method, changing the margins to allow more space for the row and column names
cim(toxicity.spls, transpose = TRUE,   
    clust.method = c("ward", "ward")), margins = c(5, 7))


# ========================================================
# CIM using mixOmics for integrative analysis, example using rCCA
# ========================================================

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
diet <- nutrimouse$diet
geno <- nutrimouse$genotype

# --------------------------
# with rCCA
# -------------------------
nutri.rcc <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# plotting 1 dataset, X
cim(nutri.rcc, mapping="X")

# plotting 1 dataset, X (lipids), changing colors
col.lipid= color.jet(n = ncol(X))     
cim(nutri.rcc, mapping = "X", xlab = "lipids", 
    ylab = "samples", x.sideColors = col.lipid, margins = c(6, 5))

# plotting 1 dataset, Y (genes), changing colors
col.sample=color.jet(n = ncol(Y)) 
cim(nutri.rcc, mapping = "Y", xlab = "genes",
ylab = "samples", y.sideColors = col.sample, margins = c(6, 5))

# plotting two data sets as a cross-correlation, see our Data Mining paper for more details on how
# the similarity matrix is computed
cim(nutri.rcc, mapping="XY")

# other options: adding label axes
cim(nutri.rcc, xlab = "genes", ylab = "lipids", margins = c(5, 6))

# In a Windows and Linux environment (not RStudio Mac), you can select the region and "see" the zoom-out region
# you need to click on 'stop' to exit
# FB: une fois que je lance ca, je sais vraiment pas comment en sortir. j ai pas de stop ni exit button, c est ou?
cim(nutri.rcc, xlab = "genes", ylab = "lipids", margins = c(5, 6), zoom =TRUE)


# By default the clustering method is set to Complete for both rows and columns
# -----------------------
# here we set it to Ward, dispaying only the X data set
# changing the side color with respect to Diet
diet.col <- color.mixo(as.numeric(nutrimouse$diet))
cim(nutri.rcc, mapping = "X", sample.names = nutrimouse$diet,
    sample.sideColors = diet.col, xlab = "lipids",
    clust.method = c("ward", "ward"), margins = c(6, 4))

# For the Y data set
geno.col = c("orange", "magenta")[as.numeric(nutrimouse$genotype)]
cim(nutri.rcc, mapping = "Y", sample.names = nutrimouse$genotype,
    sample.sideColors = geno.col, xlab = "genes",
    clust.method = c("ward", "ward"))


# only calculates the similarity matrix based on the variate from comp = 1
cim(nutri.rcc, comp = 1)

# only calculates the similarity matrix based on the variate from comp = 1
cim(nutri.rcc, comp = 2)

# ----------------------------
# with sPLS
# ---------------------------
data(liver.toxicity)

X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
liver.spls <- spls(X, Y, ncomp = 3,
                      keepX = c(20, 50, 50), keepY = c(10, 10, 10))


# default
cim(liver.spls)

# transpose matrix, choose clustering method
cim(liver.spls, transpose = TRUE,   
    clust.method = c("ward", "ward"), margins = c(5, 7))


# Here we visualise only the X variables selected 
cim(liver.spls, mapping="X")

# Here we should visualise only the Y variables selected
cim(liver.spls, mapping="Y") 

# Here we only visualise the similarity matrix between the variables by spls  
cim(liver.spls, cluster="none")

# plotting two data sets with the similarity matrix as input in the funciton (see our Data Mining paper for more details)
# Only the variables selected by the sPLS model in X and Y are represented
cim(liver.spls, mapping="XY")

# on the X matrix only, side col var to indicate dose
dose.col <- color.mixo(as.numeric(as.factor(liver.toxicity$treatment[, 3])))
cim(liver.spls, mapping = "X", sample.sideColors = dose.col, 
    sample.names = liver.toxicity$treatment[, 3])

# CIM default representation includes the total of 120 genes selected, with the dose color
# with a sparse method, show only the variables selected on specific components
cim(liver.spls, comp = 1)
cim(liver.spls, comp = 2)
cim(liver.spls, comp = c(1,2))
cim(liver.spls, comp = c(1,3))






# ----------------------------------------------------------------
## CIM representation for objects of class splsda 'multilevel' 
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


# The row side bar indicates the two levels of the facteor, stimulation and time.
# the sample names have been motified on the plot.
cim(res.2level, sample.sideColors = cbind(stim.col, time.col), 
    sample.names = paste(design$time, design$stim, sep = "_"),
    var.names = FALSE,
  #setting up legend:
    legend=list(legend = unique(design$stim), col = unique(stim.col), title = "Condition", cex = 0.9)
)


# ----------------------------------------------------------------
## CIM representation for objects of class spls 'multilevel' 
#------------------------------------------------------------------

data(liver.toxicity)
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)

# sPLS is a non supervised technique, and so we only indicate the sample repetitions in the design (1 factor only here, sample)
# sPLS takes as an input 2 data sets, and the variables selected
design <- data.frame(sample = repeat.indiv) 
res.spls.1level <- multilevel(X = liver.toxicity$gene,
                              Y=liver.toxicity$clinic,
                              design = design,
                              ncomp = 2,
                              keepX = c(50, 50), keepY = c(5, 5),
                              method = 'spls', mode = 'canonical')

stim.col <- c("darkblue", "purple", "green4","red3")

# showing only the Y variables, and only those selected in comp 1 
cim(res.spls.1level, mapping="Y",
    sample.sideColors = stim.col[factor(liver.toxicity$treatment[,3])], comp = 1,
    #setting up legend:
    legend=list(legend = unique(liver.toxicity$treatment[,3]), col=stim.col, title = "Dose", cex=0.9))


# showing only the X variables, for all selected on comp 1 and 2 
cim(res.spls.1level, mapping="X",
    sample.sideColors = stim.col[factor(liver.toxicity$treatment[,3])], 
    #setting up legend:
    legend=list(legend = unique(liver.toxicity$treatment[,3]), col=stim.col, title = "Dose", cex=0.9))


# These are the cross correlations between the variables selected in X and Y.
# The similarity matrix is obtained as in our paper in Data Mining
cim(res.spls.1level, mapping="XY")





