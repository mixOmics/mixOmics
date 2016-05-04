#########Meeting 13/10

########Issue 41

####remove playwith in plot.perf ; ok

######Example cim where similarity > 1

data(nutrimouse) 
X <- nutrimouse$gene 
Y <- nutrimouse$lipid[,1:5]
res.pls1 <- pls(X, Y, ncomp = 10, mode = "regression") 
cim(res.pls1)

######plotVar: add plot= TRUE argument 

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

plotVar(nutri.res)
plotVar(nutri.res,plot=FALSE)

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

plotVar(toxicity.spls)
plotVar(toxicity.spls,plot=FALSE)

data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(nutrimouse$gene, nutrimouse$lipid,Y)

design = matrix(c(0,1,1,
                  1,0,1,
                  1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


wrap.result.sgcca = wrapper.sgcca(blocks = data, design = design, penalty = c(.3,.3, 1),
                                  ncomp = c(2, 2, 1),
                                  scheme = "centroid", verbose = FALSE)

plotVar(wrap.result.sgcca, comp = c(1,2), 
        blocks = c(1,2), comp.select = c(1,1), 
        overlap = FALSE,
        main = 'Variables selected on component 1 only')

plotVar(wrap.result.sgcca, comp = c(1,2), 
        blocks = c(1,2), comp.select = c(1,1), 
        overlap = FALSE,plot=FALSE,
        main = 'Variables selected on component 1 only')


####check that by default comp.select = comp and that by default comp = c(1,2) : ok

#####by default pls, rcc and gcca should show all spaces 

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

plotIndiv(nutri.res) 
plotIndiv(nutri.res,rep.space = "X-variate") 

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50), 
                      keepY = c(10, 10, 10))

plotIndiv(toxicity.spls)
plotIndiv(toxicity.spls,rep.space = "X-variate") 


data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
nutrimouse.sgcca <- wrapper.sgcca(blocks = data,
                                  design = design1,
                                  penalty = c(0.3, 0.5, 1),
                                  ncomp = c(2, 2, 3),
                                  scheme = "centroid",
                                  verbose = FALSE, 
                                  bias = FALSE)

plotIndiv(nutrimouse.sgcca)


data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

splsda.breast <- splsda(X, Y,keepX=c(10,10),ncomp=2)
plotIndiv(splsda.breast)


######remove that title in all single omics methods

data(multidrug)
pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, 
          group = as.numeric(as.factor(multidrug$cell.line$Class)),main="titel")

######object == sgccda then do not plot the Y block

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)

nutrimouse.sgccda1 <- wrapper.sgccda(blocks = data,
                                     Y = Y,
                                     design = design1,
                                     ncomp = c(3, 3),
                                     keep = list(c(10,10,10), c(15,15,15)),
                                     scheme = "centroid",
                                     verbose = FALSE,
                                     bias = FALSE)

plotIndiv(nutrimouse.sgccda1)

#####Florian but useless

plotIndiv(nutrimouse.sgccda1,style='3d',plot.ellipse = TRUE)
plotIndiv(nutrimouse.sgccda1,style='3d',plot.ellipse = TRUE,alpha=0.9)
plotIndiv(nutrimouse.sgccda1,style='3d',plot.ellipse = TRUE,alpha=0.1)

############plotArrow

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, mode = "canonical", ncomp = 3, 
                      keepX = c(50, 50, 50), keepY = c(10, 10, 10))

color.toxicity <- as.numeric(liver.toxicity$treatment[, 3])
label.toxicity <- liver.toxicity$treatment[, 1]
plotArrow(toxicity.spls$variates$X[, c(1, 2)], 
          toxicity.spls$variates$Y[, c(1, 2)], arrows =T,points=T, label=label.toxicity,
          cex.lab = 0.5, color =  color.toxicity)

plotArrow(toxicity.spls$variates$X[, c(1, 2)], 
          toxicity.spls$variates$Y[, c(1, 2)], points = T,arrows=T,
          cex.lab = 0.5, color =  color.toxicity,group=as.factor(liver.toxicity$treatment[, 3]))

plotArrow(toxicity.spls$variates$X[, c(1, 2)], 
          toxicity.spls$variates$Y[, c(1, 2)], star=T,points = T,arrows=T,
          cex.lab = 0.5, color =  color.toxicity,group=as.factor(liver.toxicity$treatment[, 3]))


