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

plotVar(nutri.res) #(default)

plotVar(nutri.res, comp = c(1,3), cutoff = 0.5)

## variable representation for objects of class 'pls' or 'spls'
# ----------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
keepY = c(10, 10, 10))

plotVar(toxicity.spls, cex = c(1,0.8))

## variable representation for objects of class 'splsda'
# ----------------------------------------------------

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- as.factor(liver.toxicity$treatment[, 4])

ncomp <- 3
keepX <- rep(20, ncomp)

splsda.liver <- splsda(X, Y, ncomp = ncomp, keepX = keepX)
plotVar(splsda.liver)


## variable representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------
## see example in ??wrapper.sgcca
data(nutrimouse)
# need to unmap the Y factor diet
Y = unmap(nutrimouse$diet)
# set up the data as list
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)

# set up the design matrix:
# with this design, gene expression and lipids are connected to the diet factor
# design = matrix(c(0,0,1,
#                   0,0,1,
#                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
1,0,1,
1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


#note: the penalty parameters will need to be tuned
wrap.result.sgcca = wrapper.sgcca(X = data, design = design, penalty = c(.3,.3, 1),
ncomp = 2,
scheme = "centroid")
wrap.result.sgcca

#variables selected on component 1 for each block
selectVar(wrap.result.sgcca, comp = 1, block = c(1,2))$'gene'$name
selectVar(wrap.result.sgcca, comp = 1, block = c(1,2))$'lipid'$name

#variables selected on component 2 for each block
selectVar(wrap.result.sgcca, comp = 2, block = c(1,2))$'gene'$name
selectVar(wrap.result.sgcca, comp = 2, block = c(1,2))$'lipid'$name

plotVar(wrap.result.sgcca, comp = c(1,2), block = c(1,2), comp.select = c(1,1),title = c('Variables selected on component 1 only'))
plotVar(wrap.result.sgcca, comp = c(1,2), block = c(1,2), comp.select = c(2,2),title = c('Variables selected on component 2 only'))

# -> this one shows the variables selected on both components
plotVar(wrap.result.sgcca, comp = c(1,2), block = c(1,2), title = c('Variables selected on components 1 and 2'))

## variable representation for objects of class 'rgcca'
# ----------------------------------------------------
data(nutrimouse)
# need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
Y = unmap(nutrimouse$diet)

data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
byrow = TRUE, dimnames = list(names(data), names(data)))

nutrimouse.rgcca <- wrapper.rgcca(X = data,
design = design,
tau = "optimal",
ncomp = 2,
scheme = "centroid")

plotVar(nutrimouse.rgcca, comp = c(1,2), block = c(1,2), cex = c(1.5, 1.5))
plotVar(nutrimouse.rgcca, comp = c(1,2), block = c(1,2))


# set up the data as list
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid,Y = Y)
# with this design, gene expression and lipids are connected to the diet factor
# design = matrix(c(0,0,1,
#                   0,0,1,
#                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

# with this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
1,0,1,
1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#note: the tau parameter is the regularization parameter
wrap.result.rgcca = wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0),
ncomp = 2,
scheme = "centroid")
#wrap.result.rgcca
plotVar(wrap.result.rgcca, comp = c(1,2), block = c(1,2))

#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    
    ####### SPLS
    
    data(liver.toxicity)
    
    X <- liver.toxicity$gene
    Y <- liver.toxicity$clinic
    toxicity.spls <- spls(X, Y, ncomp = 3,
    keepX = c(30, 30, 30), keepY = c(10, 10, 10))
    
    
    plotVar(toxicity.spls)
    
    # test comp
    plotVar(toxicity.spls,comp=c(1,3))
    plotVar(toxicity.spls,comp=c(3,1))
    plotVar(toxicity.spls,comp=c(2,3))
    
    # test var.names
    plotVar(toxicity.spls,var.names = c(FALSE,TRUE))
    plotVar(toxicity.spls,var.names = c(TRUE,FALSE))
    plotVar(toxicity.spls,var.names = c(FALSE,FALSE))
    
    # test xlab, ylab
    plotVar(toxicity.spls,X.label="alpha")
    plotVar(toxicity.spls,Y.label="beta")
    plotVar(toxicity.spls,X.label="alpha",Y.label="beta")
    
    # test col
    plotVar(toxicity.spls,col=c("green","red"))
    col.x=rep(c("blue","red"),each=1558)
    col.y=rep("green",10)
    plotVar(toxicity.spls,col=list(col.x,col.y))
    plotVar(toxicity.spls,col=list(col.x,col.y),var.names = c(TRUE,FALSE))
    plotVar(toxicity.spls,col=list(col.x,col.y),var.names = c(FALSE,TRUE))
    
    # test cex
    plotVar(toxicity.spls,cex=c(3.5,5.5))
    cex.x=rep(3,3116)
    cex.y=rep(5,10)
    plotVar(toxicity.spls,cex=list(cex.x,cex.y))
    plotVar(toxicity.spls,col=list(col.x,col.y),cex=list(cex.x,cex.y))
    
    # test font
    plotVar(toxicity.spls,font=c(1,2))
    plotVar(toxicity.spls,font=c(1,3))
    plotVar(toxicity.spls,font=c(2,2))
    plotVar(toxicity.spls,font=c(4,2))
    
    # test rad.in
    plotVar(toxicity.spls,rad.in=0.1)
    plotVar(toxicity.spls,rad.in=0.9)
    plotVar(toxicity.spls,rad.in=0)
    plotVar(toxicity.spls,rad.in=1)
    
    # test cutoff
    plotVar(toxicity.spls,cutoff=0.6)
    plotVar(toxicity.spls,cutoff=0.2)
    
    #test overlap
    plotVar(toxicity.spls,overlap = F)
    plotVar(toxicity.spls,overlap = F, legend=TRUE)


    # test legend ######
    plotVar(toxicity.spls,legend=FALSE)
    plotVar(toxicity.spls,legend=TRUE)
    
    
    
    ########### RCC
    
    data(linnerud)
    X <- linnerud$exercise
    Y <- linnerud$physiological
    linn.res <- rcc(X, Y)
    
    plotVar(linn.res)
    
    # test comp
    plotVar(linn.res,comp=c(1,2))
    plotVar(linn.res,comp=c(2,1))
    
    # test X.var.names, Y.var.names
    plotVar(linn.res,var.names = c(TRUE,FALSE))
    plotVar(linn.res,var.names = c(FALSE,TRUE))
    
    # test col
    plotVar(linn.res,col=c("green","red"))
    col.x=rep("blue",3)
    col.y=rep("green",3)
    plotVar(linn.res,col=list(col.x,col.y))
    plotVar(linn.res,col=list(col.x,col.y),var.names = c(FALSE,FALSE))
    
    
    # test cex
    plotVar(linn.res,cex=c(3.5,5.5))
    cex.x=rep(c(2.8,4.2),c(2,1))
    cex.y=rep(5.5,3)
    plotVar(linn.res,cex=list(cex.x,cex.y))
    plotVar(linn.res,col=list(col.x,col.y),cex=list(cex.x,cex.y))
    
    # test font
    plotVar(linn.res,font=c(1,2))
    plotVar(linn.res,font=c(1,3))
    plotVar(linn.res,font=c(2,2))
    plotVar(linn.res,font=c(4,2))
    
    # test rad.in
    plotVar(linn.res,rad.in=0.1)
    plotVar(linn.res,rad.in=0.9)
    plotVar(linn.res,rad.in=0)
    plotVar(linn.res,rad.in=1)
    
    
    # test legend
    plotVar(linn.res,legend=FALSE)
    plotVar(linn.res,legend=TRUE)
    
    
    ######### SPCA
    
    data(liver.toxicity)
    rat.spca <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))
    
    plotVar(rat.spca)
    plotVar(rat.spca,col=list(rep(c("green","red"),each=1558)))
    ######### SGCCA
    
    data(nutrimouse)
    diet = unmap(nutrimouse$diet)
    blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = 3)
    
    plotVar(nutri.sgcca)
    
    # test comp
    plotVar(nutri.sgcca,comp=c(1,3))
    plotVar(nutri.sgcca,comp=c(2,3))
    plotVar(nutri.sgcca,comp=c(3,2))
    
    #test blocks
    plotVar(nutri.sgcca,comp=c(1,2), blocks = c(1,3))
    plotVar(nutri.sgcca,comp=c(1,3), blocks= 2:1)
    plotVar(nutri.sgcca,comp=c(1,2), blocks= c("diet","lipid"))
    plotVar(nutri.sgcca,comp=c(1,2),blocks=2:3)
    plotVar(nutri.sgcca,blocks=1:3)
    
    #test var.names
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,FALSE))
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE),blocks=2:3)
    plotVar(nutri.sgcca,var.names = TRUE) ######
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,TRUE),blocks=1:3)
    plotVar(nutri.sgcca,var.names = c(FALSE,FALSE,TRUE),blocks=1:3)
    
    #test col, cex, pch
    plotVar(nutri.sgcca,col=c("yellow","pink"),blocks=1:2)
    plotVar(nutri.sgcca,col=c("yellow","pink","chartreuse"))
    
    plotVar(nutri.sgcca,blocks=1:3,cex=c(3.5,4.5,6.5))
    plotVar(nutri.sgcca,blocks=1:3,pch=15:17)
    plotVar(nutri.sgcca,blocks=1:3,pch=15:17,cex=c(3.5,4.5,6.5))
    
    #test legend
    plotVar(nutri.sgcca,legend=TRUE)
    
    
    
    
    
    ############# lattice
    
    ####### SPLS
    
    
    plotVar(toxicity.spls, style = "lattice")
    
    # test comp
    plotVar(toxicity.spls,comp=c(1,3), style = "lattice")
    plotVar(toxicity.spls,comp=c(3,1), style = "lattice")
    plotVar(toxicity.spls,comp=c(2,3), style = "lattice")
    
    # test X.var.names, Y.var.names
    plotVar(toxicity.spls,var.names = c(TRUE,FALSE), style = "lattice")
    
    
    # test xlab, ylab
    plotVar(toxicity.spls,X.label="alpha",Y.label="beta", style = "lattice")
    
    # test col
    plotVar(toxicity.spls,col=c("green","red"), style = "lattice")
    col.x=rep(c("blue","red"),each=1558)
    col.y=rep("green",10)
    plotVar(toxicity.spls,col=list(col.x,col.y))
    plotVar(toxicity.spls,col=list(col.x,col.y),var.names = TRUE, style = "lattice")
    plotVar(toxicity.spls,col=list(col.x,col.y),var.names = FALSE, style = "lattice")
    
    # test cex
    plotVar(toxicity.spls,cex=c(0.5,1.5), style = "lattice")
    cex.x=rep(c(0.8,1.2),1558)
    cex.y=rep(1.5,10)
    plotVar(toxicity.spls,cex=list(cex.x,cex.y), style = "lattice")
    plotVar(toxicity.spls,col=list(col.x,col.y),cex=list(cex.x,cex.y), style = "lattice")
    
    # test font
    plotVar(toxicity.spls,font=c(1,2), style = "lattice")
    plotVar(toxicity.spls,font=c(1,3), style = "lattice")
    plotVar(toxicity.spls,font=c(2,2), style = "lattice")
    plotVar(toxicity.spls,font=c(4,2), style = "lattice")
    # test rad.in
    plotVar(toxicity.spls,rad.in=0.1, style = "lattice")
    plotVar(toxicity.spls,rad.in=0.9, style = "lattice")
    plotVar(toxicity.spls,rad.in=0, style = "lattice")
    plotVar(toxicity.spls,rad.in=1, style = "lattice")
    
    
    plotVar(toxicity.spls,overlap = FALSE, style = "lattice")
    plotVar(toxicity.spls,overlap = FALSE, style = "lattice", legend=TRUE)
    
    # test legend
    plotVar(toxicity.spls,legend=FALSE, style = "lattice")
    plotVar(toxicity.spls,legend=TRUE, style = "lattice")
    
    
    ########### RCC
    
    plotVar(linn.res, style = "lattice")
    
    # test comp
    plotVar(linn.res,comp=c(1,2), style = "lattice")
    plotVar(linn.res,comp=c(2,1), style = "lattice")
    
    # test col
    plotVar(linn.res,col=c("green","red"), style = "lattice")
    col.x=rep(c("blue","red"),c(1,2))
    col.y=rep("green",3)
    plotVar(linn.res,col=list(col.x,col.y), style = "lattice")
    plotVar(linn.res,col=list(col.x,col.y),var.names = c(TRUE,FALSE), style = "lattice")
    plotVar(linn.res,col=list(col.x,col.y),var.names = c(FALSE,TRUE), style = "lattice")
    
    
    
    ######### SPCA
    
    data(liver.toxicity)
    rat.spca <- spca(liver.toxicity$gene, ncomp = 3, keepX = rep(50, 3))
    
    plotVar(rat.spca, style = "lattice")
    plotVar(rat.spca,legend = FALSE, style = "lattice")
    plotVar(rat.spca,col=list(rep(c("green","red"),each=1558)),style="lattice")
    
    
    ######### SGCCA
    
    plotVar(nutri.sgcca, style = "lattice")
    
    # test comp
    plotVar(nutri.sgcca,comp=c(1,3), style = "lattice")
    plotVar(nutri.sgcca,comp=c(2,3), style = "lattice")
    plotVar(nutri.sgcca,comp=c(3,2), style = "lattice")
    
    #test blocks
    plotVar(nutri.sgcca,comp=c(1,2), blocks = c(1,3), style = "lattice")
    plotVar(nutri.sgcca,comp=c(1,3), blocks= 2:1, style = "lattice")
    plotVar(nutri.sgcca,comp=c(1,2), blocks= c("diet","lipid"), style = "lattice")
    plotVar(nutri.sgcca,comp=c(1,2),blocks=2:3, style = "lattice")
    plotVar(nutri.sgcca,blocks=1:3, style = "lattice")
    
    #test var.names
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,FALSE), style = "lattice")
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE),blocks=2:3, style = "lattice")
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,TRUE),blocks=1:3, style = "lattice")
    
    
    ############# graphics
    
    ####### SPLS
    
    
    plotVar(toxicity.spls, style = "graphics")
    
    # test comp
    plotVar(toxicity.spls,comp=c(1,3), style = "graphics")
    plotVar(toxicity.spls,comp=c(3,1), style = "graphics")
    plotVar(toxicity.spls,comp=c(2,3), style = "graphics")
    
    # test X.var.names, Y.var.names
    plotVar(toxicity.spls,var.names = TRUE, style = "graphics")
    plotVar(toxicity.spls,var.names = FALSE, style = "graphics")
    
    # test xlab, ylab
    plotVar(toxicity.spls,X.label="alpha", style = "graphics")
    plotVar(toxicity.spls,Y.label="beta", style = "graphics")
    
    
    # test col
    plotVar(toxicity.spls,col=c("green","red"), style = "graphics")
    col.x=rep(c("blue","red"),each=1558)
    col.y=rep("green",10)
    plotVar(toxicity.spls,col=list(col.x,col.y), style = "graphics")
    plotVar(toxicity.spls,col=list(col.x,col.y),var.names = TRUE, style = "graphics")
    plotVar(toxicity.spls,col=list(col.x,col.y),var.names = FALSE, style = "graphics")
    
    # test cex
    plotVar(toxicity.spls,cex=c(0.5,1.5), style = "graphics")
    cex.x=rep(c(0.5,1.2),1558)
    cex.y=rep(2,10)
    plotVar(toxicity.spls,cex=list(cex.x,cex.y), style = "graphics")
    plotVar(toxicity.spls,col=list(col.x,col.y),cex=list(cex.x,cex.y), style = "graphics")
    
    # test font
    plotVar(toxicity.spls,font=c(1,2), style = "graphics")
    plotVar(toxicity.spls,font=c(1,3), style = "graphics")
    plotVar(toxicity.spls,font=c(2,2), style = "graphics")
    plotVar(toxicity.spls,font=c(4,2), style = "graphics")
    # test rad.in
    plotVar(toxicity.spls,rad.in=0.1, style = "graphics")
    plotVar(toxicity.spls,rad.in=0.9, style = "graphics")
    plotVar(toxicity.spls,rad.in=0, style = "graphics")
    plotVar(toxicity.spls,rad.in=1, style = "graphics")
    
    plotVar(toxicity.spls,overlap = FALSE, style = "graphics")
    plotVar(toxicity.spls,overlap = FALSE, style = "graphics", legend=TRUE)

    # test legend
    plotVar(toxicity.spls,legend=FALSE, style = "graphics")
    plotVar(toxicity.spls,legend=TRUE, style = "graphics")
    
    
    ########### RCC
    
    data(linnerud)
    X <- linnerud$exercise
    Y <- linnerud$physiological
    linn.res <- rcc(X, Y)
    
    plotVar(linn.res, style = "graphics")
    
    # test comp
    plotVar(linn.res,comp=c(1,2), style = "graphics")
    plotVar(linn.res,comp=c(2,1), style = "graphics")
    
    # test col
    plotVar(linn.res,col=c("green","red"), style = "graphics")
    col.x=rep(c("blue","red"),c(1,2))
    col.y=rep("green",3)
    plotVar(linn.res,col=list(col.x,col.y), style = "graphics")
    plotVar(linn.res,col=list(col.x,col.y),var.names = TRUE, style = "graphics")
    plotVar(linn.res,col=list(col.x,col.y),var.names = FALSE, style = "graphics")
    
    
    ######### SGCCA
    
    plotVar(nutri.sgcca, style = "graphics")
    
    # test comp
    plotVar(nutri.sgcca,comp=c(1,3), style = "graphics")
    plotVar(nutri.sgcca,comp=c(2,3), style = "graphics")
    plotVar(nutri.sgcca,comp=c(3,2), style = "graphics")
    
    #test blocks
    plotVar(nutri.sgcca,comp=c(1,2), blocks = c(1,3), style = "graphics")
    plotVar(nutri.sgcca,comp=c(1,3), blocks= 2:1, style = "graphics")
    plotVar(nutri.sgcca,comp=c(1,2), blocks= c("diet","lipid"), style = "graphics")
    plotVar(nutri.sgcca,comp=c(1,2),blocks=2:3, style = "graphics")
    plotVar(nutri.sgcca,blocks=1:3, style = "graphics")
    
    #test var.names
    plotVar(nutri.sgcca,var.names = TRUE, style = "graphics")
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,FALSE), style = "graphics")
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE),blocks=2:3, style = "graphics")
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,TRUE),blocks=1:3, style = "graphics")
    
    #test font
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,TRUE), font= 1:3, style = "graphics")
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,TRUE), font= 2:4, style = "graphics")
    plotVar(nutri.sgcca,var.names = c(FALSE,TRUE,TRUE), font= c(1,2,1), style = "graphics")
    plotVar(nutri.sgcca,var.names = c(TRUE,FALSE,FALSE), font= c(1,2,2), style = "graphics")
    plotVar(nutri.sgcca,var.names = c(TRUE,FALSE,TRUE), font= c(4,2,1), style = "graphics")
    
    #test legend
    plotVar(nutri.sgcca,legend=FALSE)
}
par(opar)
