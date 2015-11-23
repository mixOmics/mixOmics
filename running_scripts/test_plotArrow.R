#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################

##!! NOTE from KA: the help file not updated yet.

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

plotArrow(nutri.res)

# names indicate genotype
plotArrow(nutri.res,
group = nutrimouse$genotype, ind.names  =   nutrimouse$genotype)

plotArrow(nutri.res, group = nutrimouse$genotype,
add.legend = TRUE)



## plot of individuals for objects of class 'pls' or 'spls'
# ----------------------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
keepY = c(10, 10, 10))

#default
plotArrow(toxicity.spls)

# colors indicate time of necropsy, text is the dose
plotArrow(toxicity.spls,  group = liver.toxicity$treatment[, 'Time.Group'],
ind.names  = liver.toxicity$treatment[, 'Dose.Group'],
add.legend = TRUE)

# colors indicate time of necropsy, text is the dose, label at start of arrow

plotArrow(toxicity.spls,  group = liver.toxicity$treatment[, 'Time.Group'], ind.names  = liver.toxicity$treatment[, 'Dose.Group'],
add.legend = TRUE, position.names = 'start')

## variable representation for objects of class 'sgcca' (or 'rgcca')
# ----------------------------------------------------
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

# default style: same color for all samples
plotArrow(nutrimouse.sgcca)

plotArrow(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE,
          main = 'my plot')

# ind.names to visualise the unique individuals
plotArrow(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE,
          main = 'my plot', ind.names = TRUE)

# ind.names to visualise the unique individuals
plotArrow(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE,
          main = 'my plot', ind.names = TRUE,position.names   = 'start')

plotArrow(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE,
          main = 'my plot', ind.names = TRUE,position.names   = 'end')

# ind.names indicates the diet
plotArrow(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE,
          main = 'my plot', ind.names = nutrimouse$diet, position.names= 'start')

# ind.names to visualise the unique individuals, start position
plotArrow(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE,
          main = 'my plot', ind.names = TRUE, position.names   = 'start')

# end position
plotArrow(nutrimouse.sgcca, group = nutrimouse$diet, add.legend =TRUE,
          main = 'my plot', ind.names = TRUE, position.names   = 'end')


## variable representation for objects of class 'sgccda'
# ----------------------------------------------------
# Note: the code differs from above as we use a 'supervised' GCCA analysis
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)

nutrimouse.sgccda1 <- wrapper.sgccda(blocks = data,
Y = Y,
design = design1,
ncomp = c(2, 2),
keep = list(c(10,10), c(15,15)),
scheme = "centroid",
verbose = FALSE,
bias = FALSE)


#  default colors correspond to outcome Y
plotArrow(nutrimouse.sgccda1)

# with legend and title and indiv ID
plotArrow(nutrimouse.sgccda1,  add.legend = TRUE, main = 'my sample plot', ind.names = TRUE, position.names = 'start')


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    
    
    library(mixOmics)
    
    
    ### PLS, sPLS
    
    data(liver.toxicity)
    
    X <- liver.toxicity$gene
    Y <- liver.toxicity$clinic
    toxicity.spls <- spls(X, Y, ncomp = 3,
    keepX = c(50, 50, 50), keepY = c(10, 10, 10))
    
    
    col <- c("blue", "red", "darkgreen", "darkviolet")
    pch <- c(15, 16, 17, 18)
    
    cex <- c(1, 1.5, 2, 2.5)
    
    
    ###test spls with vector ind.names and abline
    
    plotArrow(toxicity.spls, comp = c(1,2), ind.names  = rep(c("a","b","c","d"),each=16),abline.line=TRUE)
    
    ###test spls with ind.names=T,ind.names start arrows, pch and group
    
    plotArrow(toxicity.spls, comp = c(2,1),ind.names  =TRUE,position.names ='start',pch=pch,
    group=rep(c(1,2,3,4), each=16))
    
    ###test spls with ind.names=T,ind.names end arrows, cex and group
    
    plotArrow(toxicity.spls,cex=c(1, 1.5, 2, 2.5),
    position.names='end',ind.names=T,
    group=rep(c(1,2,3,4), each=16))
    
    ###change comp
    
    plotArrow(toxicity.spls, comp = c(1,3))
    
    ###test group and col with add.legend
    
    group=rep(c("group 1","group 2","group 3","group 4"),each = 16)
    
    plotArrow(toxicity.spls,add.legend=TRUE,abline.line=TRUE,group=group)
    
    plotArrow(toxicity.spls,add.legend=TRUE,ind.names =TRUE,col=col,group=group)
    
    ###test group with add.legend and position end
    
    plotArrow(toxicity.spls,add.legend=TRUE,ind.names=TRUE,position.names='end',group=group)
    
    ###test group with add.legend and position start
    
    plotArrow(toxicity.spls,ind.names=T,position.names='start',group=group,abline.line = TRUE)
    abline(0,0)
    points(5,2)
    
    ###rcc
    
    data(nutrimouse)
    
    X <- nutrimouse$lipid
    Y <- nutrimouse$gene
    
    col=c("red","blue")
    pch=c(0,1)
    cex=c(1,2)
    
    nutri.res <- rcc(X, Y, ncomp = 5, lambda1 = 0.064, lambda2 = 0.008)
    
    ###test with group, col and ind.names
    
    plotArrow(nutri.res,group=rep(c("group 1", "group 2"),each=20))
    plotArrow(nutri.res,col=col,group=rep(c("group 1", "group 2"),each=20))
    plotArrow(nutri.res,col=col,group=rep(c("group 1", "group 2"),each=20),ind.names   = T)
    
    ###test with group, col and add.legend
    
    nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 10)
    plotArrow(nutri.res,col=col,add.legend=TRUE,
    group=rep(c("group 1", "group 2"),each=20))
    
    nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 10, lambda2 = 0.008)
    plotArrow(nutri.res,main="plot",comp=c(2,3))
    
    ###test add after plot
    
    plotArrow(nutri.res)
    abline(0,0)
    plotArrow(nutri.res,add.legend=T)
    points(0,0)
    
    
    #blocks
    diet = unmap(nutrimouse$diet)
    blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = c(3, 3, 2))
    
    ###test blocks with ind.names and position end
    
    plotArrow(nutri.sgcca) ####
    plotArrow(nutri.sgcca,  ind.names=T,position.names  = 'end')
    
    ###test blocks with ind.names vector and group
    
    plotArrow(nutri.sgcca, ind.names= nutrimouse$diet,
    group=nutrimouse$diet,add.legend=TRUE)
    
    ###test ind.names  with TRUE and vector, position start, unique col , add.legend and abline
    
    plotArrow(nutri.sgcca, ind.names=T,position.names  ='start',add.legend=TRUE)
    plotArrow(nutri.sgcca, ind.names= 1:40,col="blue",abline.line=TRUE)
    plotArrow(nutri.sgcca, col="blue",abline.line=TRUE,add.legend=TRUE)
    
    
    ###test all with xlim and ylim
    
    plotArrow(toxicity.spls,ylim=c(1,5))
    plotArrow(nutri.sgcca,xlim=c(1,5))
    
}