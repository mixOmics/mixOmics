#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

## RGCCA
# --------------
library(mixOmicsv6)
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
ncomp = c(2, 2, 2),
scheme = "centroid",
verbose = FALSE)

# blocks should specify the block data set where the sample plot can be performed
# (ideally when there are >= 2 components!)
# we indicate the diet variable colors.
plotIndiv(nutrimouse.rgcca, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)

plotIndiv(nutrimouse.rgcca, blocks = c(1), group = nutrimouse$diet, plot.ellipse = TRUE)
plotIndiv(nutrimouse.rgcca, blocks = c(1), group = nutrimouse$diet, plot.ellipse = TRUE,plot.star=TRUE,plot.centroid=TRUE)


plotIndiv(nutrimouse.rgcca, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)
plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, plot.ellipse = TRUE)
plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, plot.ellipse = TRUE,style="graphics")
plot(1:10,1:10)

plotIndiv(nutrimouse.rgcca, group = nutrimouse$diet, plot.ellipse = TRUE,style="graphics",add.legend=TRUE)
plot(1:10,1:10)

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
ncomp = c(2, 2, 2),
scheme = "centroid",
verbose = FALSE)

# blocks should specify the block data set where the sample plot can be performed
# (ideally when there are >= 2 components!)
# we indicate the diet variable colors.
#plotIndiv(nutrimouse.srgcca, blocks = c(1,2), group = nutrimouse$diet, plot.ellipse = TRUE)

# have a look at the looadings
head(nutrimouse.srgcca$loadings[[1]])
head(nutrimouse.srgcca$loadings[[2]])
head(nutrimouse.srgcca$loadings[[3]])




nutrimouse.srgcca2 <- wrapper.rgcca(X = data,
#design = design,
tau = "optimal",
ncomp = c(2, 2, 2),
keepX=NULL,
scheme = "centroid",
verbose = FALSE)
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
ncomp = c(3, 3, 3),
scheme = "centroid",
verbose = FALSE,
bias = FALSE)

# In plotIndiv we indicate the diet variable colors and the blocks to be plotted
# (only blocks with comp  >=2!)
plotIndiv(nutrimouse.sgcca, blocks = c(1,2), group = nutrimouse$diet,
plot.ellipse = TRUE)

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
ncomp = c(2, 2, 2),
# for keep: each element of the list corresponds to a block
# and is of length the # comp per block
keepX = list(c(10,10), c(15,15), c(ncol(Y))),
scheme = "centroid",
verbose = FALSE,
bias = FALSE)


# In plotIndiv we indicate the diet variable colors and the blocks to be plotted
# (only blocks with comp  >=2!)
plotIndiv(nutrimouse.sgcca, blocks = c(1,2), group = nutrimouse$diet,
plot.ellipse = TRUE)

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
keepX = list(c(10,10), c(15,15)),
ncomp = c(3, 3),
scheme = "centroid",
verbose = FALSE,
bias = FALSE)



plotIndiv(nutrimouse.sgccda,style="lattice") # Amrit function
plotIndiv(nutrimouse.sgccda) # Amrit function
#plotIndiv(nutrimouse.sgccda, blocks = c(1,2), group = nutrimouse$diet,
#plot.ellipse = TRUE)

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

