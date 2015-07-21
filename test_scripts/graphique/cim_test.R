library(mixOmics)
##library(R.utils)

##sourceDirectory(path='C:/Users/Yan/Downloads/mixOmics1')

source("../../mixOmics/R/cim.R"); 
## KA: function renamed as CIM

set.seed(123)
mat <- matrix(rnorm(200), 20, 10)
mat[ 1:10, seq(1, 10, 2)] <- mat[ 1:10, seq(1, 10, 2)] - 4
mat[11:20, seq(2, 10, 2)] <- mat[11:20, seq(2, 10, 2)] + 2
mat[15:20, seq(2, 10, 2)] <- mat[15:20, seq(2, 10, 2)] + 3
mat <- mat[sample(20), ]
colnames(mat) <- paste("sample", 1:10, sep = " ")
rownames(mat) <- paste("gene", 1:20, sep = "")
dim(mat)

cond <- rep("cond 1", 10)
cond[seq(2, 10, 2)] <- "cond 2"
cond.col <- c("cond 1" = "darkviolet", "cond 2" = "darkorange")

obj.cim=cim(mat, color = color.GreenRed(51), col.sideColors = cond.col[cond],
            legend=list( legend = names(cond.col), col = cond.col,title = "Condition", cex = 0.9))

ddr <- as.hclust(obj.cim$ddr)
cl <- cutree(ddr, k = 3)
gene.col <- c("up" = "red", "down" = "green", "null" = "black")

cim(mat, color = color.GreenRed(51), col.sideColors = rep(c("red","blue"),each=5),
    row.sideColors=gene.col[cl],legend=list(legend = names(gene.col), col = gene.col,
                                            title = "Regulation", cex = 0.9))  

# testing transpose matrix

## FB: title in legend does not appear
#corrigé
cim(mat,transpose=TRUE,
    legend=list(legend = names(gene.col), col = gene.col,
                cex=1.5,title = "Regulation")) 

cim(mat,row.names=names(gene.col)[cl],row.sideColors=gene.col[cl],
    col.sideColors = cond.col[cond],col.names=names(cond.col)[rep(c(1,2),5)],
    legend=list(legend = names(gene.col), col = gene.col,
                title = "Regulation", cex = 0.5)) 


cim(mat,row.cex=0.8,col.cex=1.2,transpose=TRUE,
    legend=list(legend = names(gene.col), col = gene.col,
                title = "Regulation", cex = 0.9)) 

cim(mat,cluster="column",transpose=TRUE,
    legend=list(legend = names(gene.col), col = gene.col,
                title = "Regulation", cex = 0.9))

cim(mat,cluster="row",transpose=TRUE,
    legend=list(legend = names(gene.col), col = gene.col,
                title = "Regulation", cex = 0.9)) 



cim(mat,dist.method=c("correlation","maximum"))


# cutting tree at different levels. Interesting!
cim(mat,cut.tree=c(0.5,0.2))

# adding title of plot
cim(mat,symkey=TRUE,keysize=c(1,1),main = "My plot",xlab="sample",ylab="gene",margins=c(7,7),
    legend=list(legend = names(gene.col), col = gene.col))

cim(mat,lhei=c(3,4),lwid=c(2,5))


# ---------------
# CIM on a single data set
# --------------
data(liver.toxicity)
X <- liver.toxicity$gene

liver.spca <- spca(X, ncomp = 3, keepX = rep(30, 3), scale = FALSE)

cim(liver.spca,
    sample.names = liver.toxicity$treatment[, 3],sample.sideColors=c(1,2,3,4)[factor(liver.toxicity$treatment[, 3])],
    clust.method = c("single", "mcquitty"))
cim(liver.spca,
    sample.names = liver.toxicity$treatment[, 3],sample.sideColors=c(1,2,3,4)[factor(liver.toxicity$treatment[, 3])],
    clust.method = c("average", "centroid"))
cim(liver.spca,
    sample.names = liver.toxicity$treatment[, 3],sample.sideColors=c(1,2,3,4)[factor(liver.toxicity$treatment[, 3])],
    clust.method = c("median", "complete"))


cim(liver.spca,X.var.names=FALSE,scale=TRUE,comp=c(1,2,2))

# --------------------
# CIM for integrative analysis
# ------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
diet <- nutrimouse$diet
geno <- nutrimouse$genotype

# with rCCA
# ---------
nutri.rcc <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# plotting 1 dataset, X
cim(nutri.rcc, mapping="X")

# FB: below: dont understand what the vector col is supposed to do?
##x.sidecolor ne marche pas lorsque le mapping est Y erreur de ma part.
# plotting 1 dataset, X, changing colors
col=palette(rainbow(n = 21))
cim(nutri.rcc, mapping = "X", xlab = "lipids", 
    ylab = "samples", x.sideColors = col, margins = c(6, 5))

col=palette(rainbow(n = 120))
cim(nutri.rcc, mapping = "Y", xlab = "lipids",
ylab = "samples", y.sideColors = col, margins = c(6, 5))

# plotting two data sets
cim(nutri.rcc, mapping="XY")

# other options
cim(nutri.rcc, xlab = "genes", ylab = "lipids", margins = c(5, 6))

cim(nutri.rcc, xlab = "genes", ylab = "lipids", margins = c(5, 6) )

#-- select the region and "see" the zoom-out region -> FB: does that work? it does not in my RStudio. SHoudnt it be when interactive = TRUE
##ca marche sur windows et linux par contre pour quitter il y a un bouton arreter qu'il faut cliquer avant de quitter

#-- cim from X matrix
diet.col <- color.mixo(as.numeric(nutrimouse$diet))
cim(nutri.rcc, mapping = "X", sample.names = nutrimouse$diet,
    sample.sideColors = diet.col, xlab = "lipids",
    clust.method = c("ward", "ward"), margins = c(6, 4))

#-- cim from Y matrix
geno.col = c("orange", "magenta")[as.numeric(nutrimouse$genotype)]
cim(nutri.rcc, mapping = "Y", sample.names = nutrimouse$genotype,
    sample.sideColors = geno.col, xlab = "genes",
    clust.method = c("ward", "ward"))


# with sPLS
# ---------
data(liver.toxicity)

X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3,
                      keepX = c(20, 50, 50), keepY = c(10, 10, 10))

# FB: le code ci dessous ne devrait montrer que 20 variables.
# par defaut la CIM sample plotter les 120 variables selectionnees sur les 3 comp, meme lorsque l on demande comp = 1
# plot 1 data set, X
#Désolé mais la j'ai pas compris. J'ai réutilisé la meme formule que dans l'ancien cim.
#est que je dois faire une selection sur les variables?
cim(toxicity.spls,mapping="X",comp=c(1))

# et cela n a pas de sens de mettre comp = c(2,1) (ca veut dire quoi? normalement on commence par 1, puis 2 etc puisque l info decroissante)
##erreur d'écritutre dsl

cim(toxicity.spls,mapping="Y",comp=c(2,1)) 

cim(toxicity.spls,cluster="none")

# plotting two data sets, FB make sure there are only 20 genes in there for comp = 1
cim(toxicity.spls, mapping="XY", comp = 1)

# other options


# ----------------------------
## CIM default method
#----------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

# cross correlation between X and Y. FB: which correlation is this? (or which formula)
##C'est pour donner une matrice en argument par défaut on a juste besoin d'une matrice de corrélation
cim(cor(X, Y), cluster = "none")

# or one data set
cim(as.matrix(X), cluster = "none")


# other options: for X matrix
dose.col <- color.mixo(as.numeric(as.factor(liver.toxicity$treatment[, 3])))
cim(toxicity.spls, mapping = "X", sample.sideColors = dose.col, 
    sample.names = liver.toxicity$treatment[, 3])
# with a transpose
cim(toxicity.spls, transpose = TRUE,   
    clust.method = c("ward", "ward"), margins = c(5, 7))

# -----------------------------------------------------------------
## CIM representation for objects of class '(s)(i)pca'
#------------------------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene

spca.res <- spca(X, ncomp = 2, keepX = c(30, 30), scale = FALSE)

# FB: je sais pas pk cette option la ne donne que du rouge. A retester
#dose.col <- palette()[as.numeric(as.factor(liver.toxicity$treatment[, 3]))]
#testé pas vu de problème
dose.col <- color.mixo(as.numeric(as.factor(liver.toxicity$treatment[, 3])))

cim(spca.res, sample.sideColors = dose.col, var.names = FALSE,
    sample.names = liver.toxicity$treatment[, 3],
    clust.method = c("ward", "ward"))


# ----------------------------------------------------------------
## CIM representation for objects of class '(s)plsda' 
#------------------------------------------------------------------
Y <- liver.toxicity$treatment[, 3]

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(40, 30))

cim(splsda.liver, sample.sideColors = dose.col, sample.names = Y)


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


# FB: est ce que tu peux me rajouter une legende ici? de mon experience c est le plus difficile a mettre en place sur la cim
cim(res.2level, sample.sideColors = cbind(stim.col, time.col), 
    sample.names = paste(design$time, design$stim, sep = "_"),
    var.names = FALSE
    )


# ----------------------------------------------------------------
## CIM representation for objects of class spls 'multilevel' 
#------------------------------------------------------------------

# FB est ce que tu peux me rajouter un exemple pr un object de type spls multilevel? Merci!






