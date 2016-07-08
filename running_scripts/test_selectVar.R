
#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
# last modified: 02-03-2016
opar <- par(no.readonly = TRUE)

data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic

# example with sPCA
# ------------------
liver.spca <- spca(X, ncomp = 1, keepX = 10)
selectVar(liver.spca, comp = 1)$name
selectVar(liver.spca, comp = 1)$value

#example with sIPCA
# -----------------

liver.sipca <- sipca(X, ncomp = 3, keepX = rep(10, 3))
selectVar(liver.sipca, comp = 1)

# example with sPLS
# -----------------
liver.spls = spls(X, Y, ncomp = 2, keepX = c(20, 40),keepY = c(5, 5))
selectVar(liver.spls, comp = 2)

# example with sPLS-DA
data(srbct)   # an example with no gene name in the data
X = srbct$gene
Y = srbct$class

srbct.splsda = splsda(X, Y, ncomp = 2, keepX = c(5, 10))
select = selectVar(srbct.splsda, comp = 2)
select
# this is a very specific case where a data set has no rownames.
srbct$gene.name[substr(select$select, 2,5),]


# example with sGCCA
# -----------------
data(nutrimouse)

# ! need to unmap the Y factor
Y = unmap(nutrimouse$diet)
data = list(block1=nutrimouse$gene, block2=nutrimouse$lipid,Y=Y)
# in this design, gene expression and lipids are connected to the diet factor
# and gene expression and lipids are also connected
design = matrix(c(0,1,1,
1,0,1,
1,1,0), ncol = 3, nrow = 3, byrow = T)
#note: the penalty parameters need to be tuned
wrap.result.sgcca = wrapper.sgcca(X = data, design = design, penalty = c(.3,.3, 1),
ncomp = 2,
scheme = "centroid", verbose = FALSE)

#variables selected and loadings values on component 1 for the two blocs
selectVar(wrap.result.sgcca, comp = 1, block = c(1,2))

#variables selected on component 1 for each block
selectVar(wrap.result.sgcca, comp = 1, block = c(1,2))$'block1'$name
selectVar(wrap.result.sgcca, comp = 1, block = c(1,2))$'block2'$name

#variables selected on component 2 for each block
selectVar(wrap.result.sgcca, comp = 2, block = c(1,2))$'block1'$name
selectVar(wrap.result.sgcca, comp = 2, block = c(1,2))$'block2'$name

# loading value of the variables selected on the first block
selectVar(wrap.result.sgcca, comp = 1, block = 1)$'block1'$value


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################
if(additional.test==TRUE)
{
    data(liver.toxicity)
    X = liver.toxicity$gene
    Y = liver.toxicity$clinic

    liver.spca <- spca(X, ncomp = 1, keepX = 10)
    selectVar(liver.spca)
    #try(selectVar(liver.spca,comp=c(1,2))) # error comp>ncomp
    #try(selectVar(liver.spca,comp=2)) # error comp>ncomp



    liver.spls = spls(X, Y, ncomp = 2, keepX = c(20, 40),keepY = c(5, 5))
    selectVar(liver.spls, comp = 2)

    data(srbct)
    X = srbct$gene
    Y = srbct$class

    srbct.splsda = splsda(X, Y, ncomp = 2, keepX = c(5, 10))
    select = selectVar(srbct.splsda, comp = 2)
    select
    srbct$gene.name[substr(select$name, 2,5),]


    data(nutrimouse)

    diet = unmap(nutrimouse$diet)
    blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)

    nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = 3)
    selectVar(nutri.sgcca) ###
    selectVar(nutri.sgcca,block=c(1,2,3))
    selectVar(nutri.sgcca,block=c(3))

}
par(opar)

