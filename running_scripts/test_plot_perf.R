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
for(typ in c("p","l","b","c","o","h","s","S"))
{
  for(overla in c("all","measure"))
{
  for(leg in c("vertical","horizontal"))
{
    
    print(paste(di,mea,typ,overla,leg))
plot(error,
            dist =di,
            measure = mea,
            type=typ,
            xlab = NULL,
            ylab = NULL,
            overlay=overla,
            legend=leg)
}
}
}
}
}
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
error <- perf(srbct.splsda, validation = "Mfold", folds = 8, dist = "all", progressBar = FALSE)

par(opar)
