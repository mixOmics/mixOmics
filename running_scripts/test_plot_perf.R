# created on 12/03/15
# last modified: 02-03-2016
# Author: F.Rohart
#purpose: test the pls/plsda/spls/splsda function


#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

#source("mixOmics/R/plot.perf.R")
#source("mixOmics/R/internal_graphic.perf.R")


## validation for objects of class 'pls' (classification)
# ----------------------------------------
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic


liver.pls <- pls(X, Y, ncomp = 10)

liver.val <- perf(liver.pls, validation = "Mfold", folds = 5)

plot(liver.val)



## validation for objects of class 'plsda' (classification)
# ----------------------------------------
data(srbct)
X <- srbct$gene
Y <- srbct$class

ncomp = 5

srbct.plsda <- plsda(X, Y, ncomp = ncomp)

# with Mfold
# ---------
set.seed(45)
error <- perf(srbct.plsda, validation = "Mfold", folds = 8, dist = "all", progressBar = FALSE)

for(di in c("all","max.dist","centroids.dist","mahalanobis.dist"))
{
    for(mea in c("all", "BER" , "overall"))
    {
        for(overla in c("all","measure"))
        {
            for(leg in c("vertical","horizontal"))
            {
                #quartz()
                print(paste(di,mea,overla,leg))
                plot(error,
                dist =di,
                measure = mea,
                xlab = NULL,
                ylab = NULL,
                overlay=overla,
                legend.position=leg)
            }
        }
    }
}

print(paste(di,mea,overla,leg))
plot(error, dist ="all", measure = "all", overlay="all", legend="vertical")

## validation for objects of class 'splsda' (classification)
# ----------------------------------------
data(srbct)
X <- srbct$gene
Y <- srbct$class

ncomp = 5

srbct.splsda <- splsda(X, Y, ncomp = ncomp, keepX = rep(10, ncomp))

# with Mfold
# ---------
set.seed(45)
error <- perf(srbct.splsda, validation = "Mfold", folds = 8, dist = "all", progressBar = TRUE, nrepeat =4)

for(di in c("all","max.dist","centroids.dist","mahalanobis.dist"))
{
    for(mea in c("all", "BER" , "overall"))
    {
        for(overla in c("all","measure"))
        {
            for(leg in c("vertical","horizontal"))
            {
                
                print(paste(di,mea,overla,leg))
                plot(error,
                dist =di,
                measure = mea,
                xlab = NULL,
                ylab = NULL,
                overlay=overla,
                legend.position=leg)
            }
        }
    }
}

## validation for objects of class 'mint.splsda' (classification)
# ----------------------------------------


data(stemcells)
res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3, keepX = c(10, 5, 15),
study = stemcells$study)

out = perf(res)
for(di in c("all","max.dist","centroids.dist","mahalanobis.dist"))
{
    for(mea in c("all", "BER" , "overall"))
    {
        for(overla in c("all","measure"))
        {
            for(leg in c("vertical","horizontal"))
            {
                
                print(paste(di,mea,overla,leg))
                plot(out,
                dist =di,
                measure = mea,
                xlab = NULL,
                ylab = NULL,
                overlay=overla,
                legend.position=leg)
            }
        }
    }
}

## validation for objects of class 'sgccda' (classification)
# ----------------------------------------

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- block.splsda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,
scheme = "centroid",
verbose = FALSE,
bias = FALSE)

perf = perf(nutrimouse.sgccda, nrepeat = 5 , folds = 5)
for(di in c("all","max.dist","centroids.dist","mahalanobis.dist"))
{
    for(mea in c("all", "BER" , "overall"))
    {
        for(overla in c("all","measure"))
        {
            for(leg in c("vertical","horizontal"))
            {
                for(sd in c("TRUE","FALSE"))
                {
                    print(paste(di,mea,overla,leg))
                    plot(perf,
                    dist =di,
                    measure = mea,
                    xlab = NULL,
                    ylab = NULL,
                    overlay=overla,
                    legend.position=leg, sd =sd)
                }
            }
        }
    }
}

par(opar)

