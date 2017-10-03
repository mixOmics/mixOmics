#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

## RGCCA
# --------------
#library(mixOmics)
data(nutrimouse)
# need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
byrow = TRUE, dimnames = list(names(data), names(data)))

nutrimouse.rgcca <- wrapper.rgcca(X = data,
#design = design,
tau = "optimal",
ncomp = 2,
scheme = "centroid")

# blocks should specify the block data set where the sample plot can be performed
# (ideally when there are >= 2 components!)
# we indicate the diet variable colors.
plotIndiv(nutrimouse.rgcca)

plotIndiv(nutrimouse.rgcca, blocks = c(1,2), group = nutrimouse$diet, ellipse = TRUE)

plotIndiv(nutrimouse.rgcca, blocks = c(1), group = nutrimouse$diet, ellipse = TRUE)
plotIndiv(nutrimouse.rgcca, blocks = c(1), group = nutrimouse$diet, ellipse = TRUE,star=TRUE,centroid=TRUE)


plotIndiv(nutrimouse.rgcca, blocks = c(1,2), group = nutrimouse$diet, ellipse = TRUE)
plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, ellipse = TRUE)

#without layout, mfrow is reset
plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, ellipse = TRUE,style="graphics")
plot(1:10,1:10)

#with layout, we can combine
plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, ellipse = TRUE,style="graphics",layout=c(2,2))
plot(1:10,1:10)

plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, ellipse = TRUE,style="graphics",add.legend=TRUE)

# have a look at the looadings
head(nutrimouse.rgcca$loadings[[1]])
head(nutrimouse.rgcca$loadings[[2]])
head(nutrimouse.rgcca$loadings[[3]])




## sparse RGCCA
# --------------
data(nutrimouse)
# need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
byrow = TRUE, dimnames = list(names(data), names(data)))

nutrimouse.srgcca <- wrapper.rgcca(X = data,
#design = design,
tau = "optimal",
ncomp = 2,
scheme = "centroid")

# blocks should specify the block data set where the sample plot can be performed
# (ideally when there are >= 2 components!)
# we indicate the diet variable colors.
#plotIndiv(nutrimouse.srgcca, blocks = c(1,2), group = nutrimouse$diet, ellipse = TRUE)

# have a look at the looadings
head(nutrimouse.srgcca$loadings[[1]])
head(nutrimouse.srgcca$loadings[[2]])
head(nutrimouse.srgcca$loadings[[3]])




nutrimouse.srgcca2 <- wrapper.rgcca(X = data,
#design = design,
tau = "optimal",
ncomp = 2,
keepX=NULL,
scheme = "centroid")
# have a look at the looadings
head(nutrimouse.srgcca2$loadings[[1]])
head(nutrimouse.srgcca2$loadings[[2]])
head(nutrimouse.srgcca2$loadings[[3]])


## sGCCA
# -------------
# same data as above but sparse approach

# version 1 using the penalisation penalty criterion
# ---
nutrimouse.sgcca <- wrapper.sgcca(X = data,
design = design,
penalty = c(0.3, 0.5, 1),
ncomp = 3,
scheme = "centroid",
bias = FALSE)

# In plotIndiv we indicate the diet variable colors and the blocks to be plotted
# (only blocks with comp  >=2!)
plotIndiv(nutrimouse.sgcca, blocks = c(1,2), group = nutrimouse$diet,
ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgcca, comp = 1, block = 1)
selectVar(nutrimouse.sgcca, comp = 1, block = 2)

# variable plot on the selected variables
plotVar(nutrimouse.sgcca, col = color.mixo(1:3), cex = c(2,5,3))

# version 2 using the keep penalty criterion (number of variables to keep)
# it is a list per block and per component, need to specify all variables for the
# Y 'outcome' here
# (see below for sgccda code, which is more appropriate)
# ----
nutrimouse.sgcca <- wrapper.sgcca(X = data,
design = design,
ncomp = 2,
# for keep: each element of the list corresponds to a block
# and is of length the # comp per block
keepX = list(gene = c(10,10), lipid = c(15,15), Y = c(ncol(Y))),
scheme = "centroid",
bias = FALSE)


# In plotIndiv we indicate the diet variable colors and the blocks to be plotted
# (only blocks with comp  >=2!)
plotIndiv(nutrimouse.sgcca, blocks = c(1,2), group = nutrimouse$diet,
ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgcca, comp = 1, block = 1)
selectVar(nutrimouse.sgcca, comp = 1, block = 2)


## sGCC-DA
# -------------
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

nutrimouse.sgccda <- wrapper.sgccda(X = data,
Y = Y,
design = design,
keepX = list(gene = c(10,10), lipid = c(15,15)),
ncomp = 3,
scheme = "centroid",
bias = FALSE)



plotIndiv(nutrimouse.sgccda,style="lattice") # Amrit function
plotIndiv(nutrimouse.sgccda) # Amrit function
#plotIndiv(nutrimouse.sgccda, blocks = c(1,2), group = nutrimouse$diet,
#ellipse = TRUE)

# which variables are selected on a given component?
selectVar(nutrimouse.sgccda, comp = 1, block = 1)
selectVar(nutrimouse.sgccda, comp = 1, block = 2)

# variable plot on the selected variables
plotVar(nutrimouse.sgccda, col = color.mixo(1:2), cex = c(2,2))


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    cat("no additional tests")
    
}
par(opar)

