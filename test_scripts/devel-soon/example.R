rm(list=ls())

sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}

sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R


### RGCCA / SGCCA with mixOmics package version 5
require(mixOmics)
data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutrimouse.sgcca1 <- wrapper.rgcca(data = data,
                                   design = design1,
                                   tau = "optimal",
                                   ncomp = c(2, 2, 1),
                                   scheme = "centroid",
                                   verbose = FALSE, 
                                   bias = FALSE)
head(nutrimouse.sgcca1$loadings[[1]])
head(nutrimouse.sgcca1$loadings[[2]])
head(nutrimouse.sgcca1$loadings[[3]])

nutrimouse.sgcca1$tau

design2 = matrix(c(0,0,1,0,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
nutrimouse.sgcca2 <- wrapper.rgcca(data = data,
                                   design = design2,
                                   tau = "optimal",
                                   ncomp = c(2, 2, 1),
                                   scheme = "centroid",
                                   verbose = FALSE,
                                   bias = FALSE)
head(nutrimouse.sgcca2$loadings[[1]])
head(nutrimouse.sgcca2$loadings[[2]])
head(nutrimouse.sgcca2$loadings[[3]])

nutrimouse.sgcca2$tau

### RGCCA / SGCCA with new scripts 
#(source r script helpers.R and wrapper.sgcca.v6)
source("./R scripts/helpers.R"); source("./R scripts/wrapper.sgcca.v6.R"); 
source("./R scripts/wrappers.R"); source("./R scripts/plotIndiv.R")

nutrimouse.sgcca1.update <- wrapper.rgcca(blocks = data,
                                         design = design1,
                                         tau = "optimal",
                                         ncomp = c(2, 2, 1),
                                         scheme = "centroid",
                                         verbose = FALSE)

plotIndiv(nutrimouse.sgcca1.update, blocks = NULL)
plotIndiv(nutrimouse.sgcca1.update, blocks = c(1,2))

head(nutrimouse.sgcca1.update$loadings[[1]])
head(nutrimouse.sgcca1.update$loadings[[2]])
head(nutrimouse.sgcca1.update$loadings[[3]])

#-- Comparison with mixOmics package

all.equal(lapply(nutrimouse.sgcca1.update$loadings, abs), lapply(nutrimouse.sgcca1$loadings, abs), check.attributes = FALSE)
all.equal(lapply(nutrimouse.sgcca1.update$loadings.star, abs), lapply(nutrimouse.sgcca1$loadings.star, abs), check.attributes = FALSE)
all.equal(lapply(nutrimouse.sgcca1.update$variates, abs), lapply(nutrimouse.sgcca1$variates, abs), check.attributes = FALSE)


rm(list=ls())
data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
nutrimouse.sgcca1 <- wrapper.sgcca(data = data,
                                   design = design1,
                                   penalty = c(0.3, 0.5, 1),
                                   ncomp = c(2, 2, 1),
                                   scheme = "centroid",
                                   verbose = FALSE, 
                                   bias = FALSE)

head(nutrimouse.sgcca1$loadings[[1]])
head(nutrimouse.sgcca1$loadings[[2]])
head(nutrimouse.sgcca1$loadings[[3]])


source("./R scripts/helpers.R"); source("./R scripts/wrapper.sgcca.v6.R"); 
source("./R scripts/wrappers.R"); source("./R scripts/plotIndiv.R")

nutrimouse.sgcca1.update <- wrapper.sgcca(blocks = data,
                                   design = design1,
                                   penalty = c(0.3, 0.5, 1),
                                   ncomp = c(2, 2, 1),
                                   scheme = "centroid",
                                   verbose = FALSE, 
                                   bias = FALSE)

head(nutrimouse.sgcca1.update$loadings[[1]])
head(nutrimouse.sgcca1.update$loadings[[2]])
head(nutrimouse.sgcca1.update$loadings[[3]])

all.equal(lapply(nutrimouse.sgcca1.update$loadings, abs), lapply(nutrimouse.sgcca1$loadings, abs), check.attributes = FALSE)
all.equal(lapply(nutrimouse.sgcca1.update$loadings.star, abs), lapply(nutrimouse.sgcca1$loadings.star, abs), check.attributes = FALSE)
all.equal(lapply(nutrimouse.sgcca1.update$variates, abs), lapply(nutrimouse.sgcca1$variates, abs), check.attributes = FALSE)

rm(list=ls())
data(breast.tumors)
X <- breast.tumors$gene.exp
X <- t(na.omit(t(X)))
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)
palette(c("red", "blue"))
col.breast <- as.numeric(as.factor(Y))
plotIndiv(plsda.breast, ind.names = TRUE, col = col.breast)
legend('bottomleft', c("After", "Before"), pch = c(16, 16), 
       col = unique(col.breast), cex = 1, pt.cex = c(1.2, 1.2), 
       title = "Treatment")
palette("default")

source("./R scripts/helpers.R"); source("./R scripts/wrapper.sgcca.v6.R"); 
source("./R scripts/wrappers.R"); source("./R scripts/plotIndiv.R")
sgccda.breast <- wrapper.sgccda(blocks = list(X = X),
                                Y = Y,
                                ncomp = c(2),
                                verbose = FALSE,
                                bias = FALSE)


all.equal(lapply(sgccda.breast$loadings, abs), lapply(plsda.breast$loadings, abs), check.attributes = FALSE)
all.equal(lapply(sgccda.breast$variates, abs), lapply(plsda.breast$variates, abs), check.attributes = FALSE)


plotIndiv(sgccda.breast, blocks = "X")
plotIndiv(sgccda.breast, blocks = "X", style = "lattice")
plotIndiv(sgccda.breast, blocks = "X", style = "graphics")



data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

source("./R scripts/helpers.R"); source("./R scripts/wrapper.sgcca.v6.R"); 
source("./R scripts/wrappers.R"); source("./R scripts/plotIndiv.R")

nutrimouse.sgccda <- wrapper.sgccda(blocks = data,
                                    Y = Y,
                                    design = design1,
                                    ncomp = c(2, 2, 1),
                                    scheme = "centroid",
                                    verbose = FALSE,
                                    bias = FALSE)

head(nutrimouse.sgccda$loadings[[1]])
head(nutrimouse.sgccda$loadings[[2]])
head(nutrimouse.sgccda$loadings[[3]])

nutrimouse.sgccda <- wrapper.sgccda(blocks = data,
                                    Y = Y,
                                    design = design1,
                                    ncomp = c(2, 2),
                                    scheme = "centroid",
                                    verbose = FALSE,
                                    bias = FALSE)

head(nutrimouse.sgccda$loadings[[1]])
head(nutrimouse.sgccda$loadings[[2]])
head(nutrimouse.sgccda$loadings[[3]])

nutrimouse.sgccda <- wrapper.sgccda(blocks = data,
                                    Y = Y,
                                    design = design1,
                                    ncomp = c(2, 2),
                                    keep.blocks = list(c(10,10), c(15,15), c(5,5)),
                                    scheme = "centroid",
                                    verbose = FALSE,
                                    bias = FALSE)

head(nutrimouse.sgccda$loadings[[1]])
head(nutrimouse.sgccda$loadings[[2]])
head(nutrimouse.sgccda$loadings[[3]])

nutrimouse.sgccda <- wrapper.sgccda(blocks = data,
                                    Y = Y,
                                    design = design1,
                                    ncomp = c(2, 2),
                                    keep.blocks = list(c(10,10), c(15,15)),
                                    scheme = "centroid",
                                    verbose = FALSE,
                                    bias = FALSE)

head(nutrimouse.sgccda$loadings[[1]])
head(nutrimouse.sgccda$loadings[[2]])
head(nutrimouse.sgccda$loadings[[3]])




