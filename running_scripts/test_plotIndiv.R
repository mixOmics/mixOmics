#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

# default, only in the X space
plotIndiv(nutri.res)
plotIndiv(nutri.res,layout=c(2,1))
plotIndiv(nutri.res,layout=c(2,1),style="graphics")

par(opar)

# ellipse with respect to genotype in the XY space,
# names indicate genotype
plotIndiv(nutri.res, rep.space= 'XY-variate',
ellipse = TRUE, ellipse.level = 0.9,
group = nutrimouse$genotype, ind.names = nutrimouse$genotype)

# ellipse with respect to genotype in the XY space, with legend
plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype,
legend = TRUE)


# lattice style
plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype,
legend = TRUE, style = 'lattice')

# classic style, in the Y space
plotIndiv(nutri.res, rep.space= 'Y-variate', group = nutrimouse$genotype,
legend = TRUE, style = 'graphics')

# indicating centroid
plotIndiv(nutri.res, rep.space= 'X-variate', ind.names = FALSE,
          group = nutrimouse$genotype, centroid = TRUE)

# indicating the star and centroid
plotIndiv(nutri.res, rep.space= 'X-variate', ind.names = FALSE,
          group = nutrimouse$genotype, centroid = TRUE, star = TRUE)


# indicating the star and ellipse
plotIndiv(nutri.res, rep.space= 'X-variate', ind.names = FALSE,
          group = nutrimouse$genotype, centroid = TRUE, star = TRUE, ellipse = TRUE)


## plot of individuals for objects of class 'pls' or 'spls'
# ----------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
keepY = c(10, 10, 10))

#default
plotIndiv(toxicity.spls)

# in the Y space, colors indicate time of necropsy, text is the dose
plotIndiv(toxicity.spls, rep.space= 'Y-variate',
group = liver.toxicity$treatment[, 'Time.Group'],
ind.names = liver.toxicity$treatment[, 'Dose.Group'],
legend = TRUE)


# in the X space, with graphics style
plotIndiv(toxicity.spls, rep.space= 'X-variate',
group = liver.toxicity$treatment[, 'Time.Group'],
ind.names = liver.toxicity$treatment[, 'Dose.Group'],
legend = TRUE,title="",style="graphics")


## plot of individuals for objects of class 'plsda' or 'splsda'
# ----------------------------------------------------
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

splsda.breast <- splsda(X, Y,keepX=c(10,10),ncomp=2)

# default option: note the outcome color is included by default!
plotIndiv(splsda.breast)

# default option with no ind name: pch and color are set automatically, with legend
plotIndiv(splsda.breast, ind.names = FALSE, comp = c(1, 2), legend = TRUE)

# trying the different styles
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), 
          ellipse = TRUE, style = "ggplot2", cex = c(1, 1))
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), 
          ellipse = TRUE, style = "graphics", cex = c(1, 1))
plotIndiv(splsda.breast, ind.names = TRUE, comp = c(1, 2), 
          ellipse = TRUE, style = "lattice", cex = c(1, 1))

# with centroid and stars
plotIndiv(splsda.breast, ind.names = TRUE, centroid = TRUE, star = TRUE)



## variable representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------
data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
design1 = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
nutrimouse.sgcca <- wrapper.sgcca(X = data,
design = design1,
penalty = c(0.3, 0.5, 1),
ncomp = c(2, 2, 3),
scheme = "centroid",
verbose = FALSE,
bias = FALSE)

# default style: one panel for each block
plotIndiv(nutrimouse.sgcca)

# for the block 'lipid' with ellipse plots and legend, different styles
plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, legend =TRUE, 
          ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", title = 'my plot')
plotIndiv(nutrimouse.sgcca, style = "lattice", group = nutrimouse$diet, 
          legend = TRUE, ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid",
title = 'my plot')
plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, 
          legend = TRUE, ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid",
          title = 'my plot')

# stars and centroids
plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, 
          legend = FALSE, ellipse = TRUE, ellipse.level = 0.5, blocks = c("gene", "lipid"),
          title = 'my plot', centroid = TRUE, star = TRUE)

# stars and centroids
plotIndiv(nutrimouse.sgcca, style = "graphics", group = nutrimouse$diet, ind.names = nutrimouse$diet, 
          legend = TRUE, ellipse = TRUE, ellipse.level = 0.5, blocks = c("gene", "lipid"),
          title = 'my plot', centroid = TRUE, star = TRUE)



## variable representation for objects of class 'sgccda'
# ----------------------------------------------------
# Note: the code differs from above as we use a 'supervised' GCCA analysis
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)

nutrimouse.sgccda1 <- wrapper.sgccda(X = data,
Y = Y,
design = design1,
ncomp = c(2, 2),
keepX = list(gene=c(10,10), lipid = c(15,15)),
scheme = "centroid",mode="regression",
verbose = FALSE,
bias = TRUE,scale=TRUE,tol=1e-6)
nutrimouse.sgccda1$loadings$Y


# plotIndiv
# ----------
# displaying all blocks. bu default colors correspond to outcome Y
plotIndiv(nutrimouse.sgccda1)

# displaying only 2 blocks
plotIndiv(nutrimouse.sgccda1, blocks = c(1,2), group = nutrimouse$diet)

# with some ellipse, legend and title
plotIndiv(nutrimouse.sgccda1, blocks = c(1,2), group = nutrimouse$diet,
ellipse = TRUE, legend = TRUE, title = 'my sample plot')

# stars and centroids
plotIndiv(nutrimouse.sgccda1, style = "graphics", group = nutrimouse$diet, 
          legend = TRUE, ellipse = TRUE, ellipse.level = 0.5, 
          title = 'my plot', centroid = TRUE, star = TRUE)



#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################
### splsda, plsda

if(additional.test==TRUE)
{
    
    ### PLS, sPLS
    
    data(liver.toxicity)
    
    X <- liver.toxicity$gene
    Y <- liver.toxicity$clinic
    toxicity.spls <- spls(X, Y, ncomp = 3,
    keepX = c(50, 50, 50), keepY = c(10, 10, 10))
    
    
    col <- c("blue", "red", "darkgreen", "darkviolet")
    pch <- c(15, 16, 17, 18)
    
    cex <- c(1, 1.5, 2, 2.5)
    
    ##test ind.names, rep.spec=XY-variate, abline and group
    
    plotIndiv(toxicity.spls, comp = c(1,2), ind.names = rep(c("a","b","c","d"),each=16),
    rep.space = "XY-variate",abline=TRUE,group=rep(c(1,2,3,4), each=16))
    
    # FB: Error: Incompatible lengths for set aesthetics: shape, colour, size
    # is that supposed to be the case, if so PLEASE put COMMENT HERE.
    #KA: nothing error.
    ##test rep.spec=X-variate, X.label,pch and group
    
    
    plotIndiv(toxicity.spls, comp = c(2,1),
    rep.space = "X-variate", X.label="X1",ind.names=FALSE,pch=pch,
    group=rep(c(1,2,3,4), each=16))
    
    ##test rep.spec=Y-variate, cex,lattice, ellipse, group
    
    plotIndiv(toxicity.spls,rep.space = "Y-variate",cex=c(1, 1.5, 2, 2.5),
    ind.names=FALSE,style="lattice", ellipse=TRUE,
    group=rep(c(1,2,3,4), each=16))
    
    
    ##test rep.spec=XY-variate, lattice
    
    plotIndiv(toxicity.spls, comp = c(1,3),style="lattice",
    rep.space="XY-variate")
    
    group=rep(c("group 1","group 2","group 3","group 4"),each = 16)
    
    ##test  legend,abline,group
    
    plotIndiv(toxicity.spls,legend=TRUE,abline = TRUE,group=group)
    
    ##test  legend,col,lattice,group
    
    plotIndiv(toxicity.spls,legend=TRUE,ind.names=FALSE,col=col,style="lattice",group=group)
    
    ##test ind.names, legend,lattice and graphics,group
    
    plotIndiv(toxicity.spls,legend=TRUE,ind.names=TRUE,style="lattice",group=group)
    plotIndiv(toxicity.spls,legend=TRUE,ind.names=TRUE,style="graphics",group=group) #####
    
    ##test legend,lattice and graphics,group
    
    plotIndiv(toxicity.spls,legend=TRUE,ind.names=FALSE,group=group)
    
    ##test  graphics and group with add after points and abline
    
    plotIndiv(toxicity.spls,ind.names=FALSE,group=group,style="graphics")
    abline(0,0)
    points(5,2)
    
    ### RCC
    
    data(nutrimouse)
    
    X <- nutrimouse$lipid
    Y <- nutrimouse$gene
    
    col=c("red","blue")
    pch=c(0,1)
    cex=c(1,2)
    
    nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)
    
    ##test group, col
    
    plotIndiv(nutri.res,group=rep(c("group 1", "group 2"),each=20))
    plotIndiv(nutri.res,col=col,group=rep(c("group 1", "group 2"),each=20))
    plotIndiv(nutri.res,col=col,group=rep(c("group 1", "group 2"),each=20),ind.names = F)
    
    ##test lattice, group,col,legend and ellipse
    
    nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 10)
    plotIndiv(nutri.res,col=col,style="lattice",legend=TRUE,ellipse=TRUE,
    group=rep(c("group 1", "group 2"),each=20)) ### ellipse
    
    ##test title,Y.label
    
    nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 10, lambda2 = 0.008)
    plotIndiv(nutri.res,title="plot",Y.label="Y1",comp=c(2,3))
    
    ##test graphics,legend with add after points and abline
    # working properly if layout is done, otherwise par(mfrow) is reset to 1,1, which messes up with the abline and points. no problem when plotting only one block
    plotIndiv(nutri.res,style="graphics",layout=c(1,2))
    abline(0,0)
    plotIndiv(nutri.res,style="graphics",legend=T)
    points(0,0)
    plotIndiv(nutri.res,style="graphics")
    abline(0,0)
    
    par(mfrow=c(1,1))
    ###PLSDA
    
    data(breast.tumors)
    X = breast.tumors$gene.exp
    Y = breast.tumors$sample$treatment
    
    ##test lattice and ind.names
    
    breast.plsda = plsda(X, Y, ncomp = 3)
    plotIndiv(breast.plsda, ind.names = TRUE,style="lattice",legend=FALSE)
    
    col=unique(color.mixo(as.numeric(map(breast.plsda$ind.mat))))
    
    ##test ind.names,col,graphics,ellipse and cex
    
    plotIndiv(breast.plsda,cex=2,col=col, ind.names = TRUE,style="graphics",ellipse=TRUE)
    
    ##test ellipse and ellipse.level
    
    plotIndiv(breast.plsda, ind.names = FALSE,ellipse.level=0.5,ellipse=TRUE)
    
    ###BLOCKS
    
    diet = unmap(nutrimouse$diet)
    blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    
    
    # nutri.sgcca <- wrapper.sgcca(blocks, ncomp = c(3, 2, 1))####
    nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = c(3, 3, 3))
    
    ##test 1 block
    
    #blocks
    plotIndiv(nutri.sgcca) ####
    plotIndiv(nutri.sgcca, blocks = "lipid", ind.names = TRUE)
    
    ##test 1 block, vector ind.names, group and legend
    
    plotIndiv(nutri.sgcca, blocks = "gene", ind.names = nutrimouse$diet,
    group=nutrimouse$diet,legend=TRUE)
    
    ##test all options
    
    plotIndiv(nutri.sgcca, blocks = "gene", ind.names = FALSE,legend=TRUE,
    style="graphics",Y.label="Dim")
    plotIndiv(nutri.sgcca, blocks = "gene", ind.names = 1:40,
    style="lattice",col="blue",X.label="Dim",abline=TRUE)
    plotIndiv(nutri.sgcca, blocks = c("gene","lipid"),
    col="blue",X.label="Dim",abline=TRUE,legend=TRUE)
    
    ###SPCA
    
    liver.spca= spca(liver.toxicity$gene,ncomp=3)
    
    ##test group,ellipse, lattice and graphics, legend
    
    plotIndiv(liver.spca, ind.names=FALSE,
    group=liver.toxicity$treatment[, 3],ellipse=TRUE,style="lattice")
    plotIndiv(liver.spca, ind.names=FALSE,legend=TRUE,ellipse=TRUE,
    group=liver.toxicity$treatment[, 3],style="graphics") ###
    
    ##test group
    
    plotIndiv(liver.spca, ind.names = FALSE,
    group=liver.toxicity$treatment[, 3])
    
    ###PCA
    
    data(multidrug)
    pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)
    
    ##test X.label, Y-label and Z-label
    
    plotIndiv(pca.res)
    plotIndiv(pca.res,X.label="X")
    plotIndiv(pca.res,Y.label="Y")
    plotIndiv(pca.res,Z.label="Z",style="3d")
    
    ##test xlim,ylim, different style
    
    plotIndiv(toxicity.spls,ylim=list(c(-10,5),c(1,5)),style="graphics")
    plotIndiv(toxicity.spls,ylim=list(c(-10,5),c(1,5)),style="ggplot2")#not working
    plotIndiv(toxicity.spls,xlim=list(c(-10,5),c(1,5)),style="lattice")
    plotIndiv(nutri.sgcca,ylim=list(c(1,5),c(0,5),c(0,1)),style="graphics")
    plotIndiv(nutri.sgcca,ylim=list(c(1,5),c(0,5),c(0,1)),style="lattice")
    plotIndiv(nutri.sgcca,ylim=list(c(1,5),c(4,5),c(0,1)),style="ggplot2") #not working

    
    # testing star and centroid
    # if group not provided should throw a warning
    plotIndiv(nutrimouse.sgcca, legend =TRUE, ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", 
              title = 'my plot', star = TRUE)
    # star
    plotIndiv(nutrimouse.sgcca, group = nutrimouse$diet, legend =TRUE, ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", title = 'my plot', star = TRUE)

    # star and centroid
    plotIndiv(nutri.sgcca, group = nutrimouse$diet, legend =TRUE, ellipse = TRUE, ellipse.level = 0.5, blocks = "lipid", title = 'my plot', star = TRUE, centroid = TRUE)

    # centroid only
    plotIndiv(nutri.sgcca, group = nutrimouse$diet, legend =TRUE, ellipse = TRUE, 
              ellipse.level = 0.5, blocks = "lipid", title = 'my plot', star = FALSE, centroid = TRUE)
    
    #####MULTI
    
    data(vac18)
    X <- vac18$genes
    Y <- vac18$stimulation
    # sample indicates the repeated measurements
    design <- data.frame(sample = vac18$sample, 
                         stimul = vac18$stimulation)
    
    # multilevel sPLS-DA model
    res.1level <- splsda(X, Y = design[,2], ncomp = 3, multilevel = design[,1],
                             keepX = c(30, 137, 123))
    
    # set up colors for plotIndiv
    col.stim <- c("darkblue", "purple", "green4","red3")
    plotIndiv(res.1level, ind.names = Y, col.per.group = col.stim)
    
    data(vac18.simulated) # simulated data
    
    X <- vac18.simulated$genes
    design <- data.frame(sample = vac18.simulated$sample,
                         stimu = vac18.simulated$stimulation,
                         time = vac18.simulated$time)
    
    res.2level <- splsda(X, Y = design[,2], ncomp = 2, multilevel = design[,1],
                             keepX = c(200, 200))
    
    plotIndiv(res.2level, group = design$stimu, ind.names = vac18.simulated$time,
              legend = TRUE, style = 'lattice')
    
    data(liver.toxicity)
    # note: we made up those data, pretending they are repeated measurements
    repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                      6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                      10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                      13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
    summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each
    
    # this is a spls (unsupervised analysis) so no need to mention any factor in design
    # we only perform a one level variation split
    design <- data.frame(sample = repeat.indiv)
    res.spls.1level <- spls(X = liver.toxicity$gene,
                                  Y=liver.toxicity$clinic,
                                  multilevel = design,
                                  ncomp = 3,
                                  keepX = c(50, 50, 50), keepY = c(5, 5, 5),
                                  mode = 'canonical')
    
    # set up colors and pch for plotIndiv
    col.stimu <- 1:nlevels(design$stimu)
    
    plotIndiv(res.spls.1level, rep.space = 'X-variate', ind.names = FALSE, 
              group = liver.toxicity$treatment$Dose.Group,
              pch = 20, title = 'Gene expression subspace',
              legend = TRUE)
    
    
    plotIndiv(res.spls.1level, rep.space = 'Y-variate', ind.names = FALSE,
              group = liver.toxicity$treatment$Dose.Group,
              pch = 20, title = 'Clinical measurements ssubpace',
              legend = TRUE)
    
    plotIndiv(res.spls.1level, rep.space = 'XY-variate', ind.names = FALSE,
              group = liver.toxicity$treatment$Dose.Group,
              pch = 20, title = 'Both Gene expression and Clinical subspaces',
              legend = TRUE)

    par(mfrow=c(1,1))
}
par(opar)

