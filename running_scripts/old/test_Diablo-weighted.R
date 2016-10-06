#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

get.confusion_matrix = function(Y.learn,Y.test,pred)
{
    ClassifResult = array(0,c(nlevels(factor(Y.learn)),nlevels(factor(Y.learn))))
    rownames(ClassifResult) = levels(factor(Y.learn))
    colnames(ClassifResult) = paste("predicted.as.",levels(factor(Y.learn)),sep = "")
    #--------record of the classification accuracy for each level of Y
    for(i in 1:nlevels(factor(Y.learn)))
    {
        ind.i = which(Y.test == levels(factor(Y.learn))[i])
        for(ij in 1:nlevels(factor(Y.learn)))
        {
            ClassifResult[i,ij] = sum(pred[ind.i] == levels(Y.learn)[ij])
            
        }
    }
    ClassifResult
}

get.BER = function(X)
{
    if(!is.numeric(X)| !is.matrix(X) | length(dim(X)) != 2 | nrow(X)!=ncol(X))
    stop("'X' must be a square numeric matrix")
    
    nlev = nrow(X)
    #calculation of the BER
    ClassifResult.temp = X
    diag(ClassifResult.temp) = 0
    BER = sum(apply(ClassifResult.temp,1,sum,na.rm = TRUE)/apply(X,1,sum,na.rm = TRUE),na.rm = TRUE)/nlev
    return(BER)
}


get.result = function(X.train, Y.train, ncomp, keepX, design, X.test, Y.test)
{
    
    model = block.splsda(X = X.train, Y = Y.train, ncomp = ncomp, keepX = keepX, design = design)

    pred.train <- predict(object = model, newdata = X.train, method = "all")
    pred.test <- predict(object = model, newdata = X.test, method = "all")
    
    # pred$class gives the predicted class for each dataset
    # pred$vote gives the majority vote.
    
    # correlation of each dataset to the variates of Y
    x.xList <- list()
    compt = 1
    for(comp in 1:ncomp)
    {
        for(i in 1:length(model$variates)){
            corDat <- rep(0, length(model$variates))
            names(corDat) <- paste("cor", names(model$variates)[i], names(model$variates), sep = "_")
            for(j in (i):length(model$variates)){
                corDat[j] <- as.numeric(cor(model$variates[[i]][,comp], model$variates[[j]][,comp]))
            }
            x.xList[[compt]] <- corDat
            compt = compt +1
        }
    }
    corMat.diablo <- do.call(rbind, x.xList)
    rownames(corMat.diablo) <- paste(names(model$variates),".comp",rep(1:ncomp,each=length(model$variates)),sep="")
    colnames(corMat.diablo) <- names(model$variates)
    
    temp = matrix(corMat.diablo[,"Y"],ncol=2)
    correlation = apply(temp, 1, function(x){mean(abs(x))})[1:length(X.train)]
    names(correlation) = names(X.train)
    
    
    ##### on the training set
    class = sapply(pred.train$class$max.dist,function(x) x[,ncomp])
    
    
    #aggregate(correlation,list(class["A0T1",]),sum)
    
    pred.weighted.vote = apply(class,1,function(x){
        temp = aggregate(correlation,list(x),sum)
        ind = which(temp[,2]== max (temp[,2]))
        if(length(ind) == 1)
        {
            res = temp[ind, 1]
        } else {
            res = NA
        }
    })
    
    #differences between pred.train$vote and pred.weighted.vote
    
    temp = cbind(pred.train$vote$max.dist[,ncomp],pred.weighted.vote)
    temp = replace(temp,which(is.na(temp)),"NA")
    
    sum(apply(temp,1,function(x){as.character(x[1])!=as.character(x[2])}),na.rm=TRUE)
    
    # comparison pred.weighted.vote with the truth Y
    pred.train$vote$max.dist[,ncomp] = replace(pred.train$vote$max.dist[,ncomp],which(is.na(pred.train$vote$max.dist[,ncomp])),"NA")
    
    get.confusion_matrix(Y,Y,pred.weighted.vote)
    BER.train.weighted = get.BER(get.confusion_matrix(Y,Y,pred.weighted.vote))
    
    get.confusion_matrix(Y,Y,pred.train$vote$max.dist[,ncomp])
    BER.train = get.BER(get.confusion_matrix(Y,Y,pred.train$vote$max.dist[,ncomp]))
    
    
    
    ##### on the test set
    class = sapply(pred.test$class$max.dist,function(x) x[,ncomp])
    
    correlation.test = correlation[which(names(correlation)%in%colnames(class))]
    pred.weighted.vote = apply(class,1,function(x){
        temp = aggregate(correlation.test,list(x),sum)
        ind = which(temp[,2]== max (temp[,2]))
        if(length(ind) == 1)
        {
            res = temp[ind, 1]
        } else {
            res = NA
        }
    })
    
    #differences between pred.train$vote and pred.weighted.vote
    
    temp = cbind(pred.test$vote$max.dist[,ncomp],pred.weighted.vote)
    temp = replace(temp,which(is.na(temp)),"NA")
    
    sum(apply(temp,1,function(x){as.character(x[1])!=as.character(x[2])}),na.rm=TRUE)
    
    
    # comparison pred.weighted.vote with the truth Y
    pred.test$vote$max.dist[,ncomp] = replace(pred.test$vote$max.dist[,ncomp],which(is.na(pred.test$vote$max.dist[,ncomp])),"NA")
    
    get.confusion_matrix(Y,Y.test,pred.weighted.vote)
    BER.test.weighted = get.BER(get.confusion_matrix(Y,Y.test,pred.weighted.vote))
    
    get.confusion_matrix(Y,Y.test,pred.test$vote$max.dist[,ncomp])
    BER.test = get.BER(get.confusion_matrix(Y,Y.test,pred.test$vote$max.dist[,ncomp]))
    
    
    BER = cbind(c(BER.train.weighted,BER.test.weighted), c(BER.train,BER.test))
    colnames(BER) = c("weighted", "non-weighted")
    rownames(BER) = c("training", "test")
    
    return(list(BER = BER, design = model$design))
}


## sGCC-DA
# -------------
library(mixOmics)
#source("mixOmics/R/plotIndiv.R")
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- wrapper.sgccda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = 2,#c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE)

set.seed(43)
a=perf(nutrimouse.sgccda)

set.seed(43)
a2=perf(nutrimouse.sgccda,cpus=4)

plot(a)


load("/Users/florian/Work/git/multi-group/trainTestDatasetsNormalized.RDATA")

Y <- droplevels(pam50Train0$Call)
names(Y) <- rownames(pam50Train0)
X <- list(mRNA = mrnaTrain0, miRNA = mirnaTrain0, CpGs = methTrain0, Proteins = protTrain0)
all(names(Y) == rownames(X[[1]]))
all(names(Y) == rownames(X[[2]]))
all(names(Y) == rownames(X[[3]]))
all(names(Y) == rownames(X[[4]]))
dim(X[[1]]); dim(X[[2]]); dim(X[[3]]); dim(X[[4]]);

ncomp = 3

keep = 3
keepX = list(mRNA=rep(keep,ncomp), miRNA=rep(keep,ncomp), CpGs=rep(keep,ncomp), Proteins=rep(keep,ncomp))

## Validate in Test set
Y.test0 <- factor(pam50Test0[, "Call"])
names(Y.test0) <- rownames(pam50Test0)
Y.test <- droplevels(Y.test0[Y.test0 != "Normal"])
X.test <- list(mRNA = mrnaTest0, miRNA = mirnaTest0,
CpGs = methTest0)



# design from the paper
design <- matrix(c(0, 1, 1, 1,
1, 0, 0, 0,
1, 0, 0, 0,
1, 0, 0, 0), nrow = 4, ncol = 4, dimnames = list(names(X), names(X)))


## weighting with the pls correlation from the first component
x.xList <- list()
for(i in 1:length(X)){
    corDat <- rep(0, length(X))
    names(corDat) <- paste("cor", names(X)[i], names(X), sep = "_")
    for(j in 1:length(X)){
        result <- pls(X = X[[i]], Y = X[[j]], ncomp = 1)
        corDat[j] <- as.numeric(cor(result$variates$X, result$variates$Y))
    }
    x.xList[[i]] <- corDat
}
correlation.design.comp1 <- do.call(rbind, x.xList)
rownames(correlation.design.comp1) <- colnames(correlation.design.comp1) <- c("mRNA", "miRNA", "CpGs", "Proteins")

## weighting with the pls correlation from all component
x.xList <- list()
for(i in 1:length(X)){
    corDat <- rep(0, length(X))
    names(corDat) <- paste("cor", names(X)[i], names(X), sep = "_")
    for(j in 1:length(X)){
        result <- pls(X = X[[i]], Y = X[[j]], ncomp = ncomp)
        corDat[j] <- as.numeric(mean(diag(cor(result$variates$X, result$variates$Y)))) # average the correlation over the ncomp components
    }
    x.xList[[i]] <- corDat
}
correlation.design.allcomp <- do.call(rbind, x.xList)
rownames(correlation.design.allcomp) <- colnames(correlation.design.allcomp) <- c("mRNA", "miRNA", "CpGs", "Proteins")


# full design
full.design = matrix(1,nrow=4,ncol=4)-diag(4)



table.full.design = get.result(X.train = X, Y.train = Y, ncomp = 4, keepX = keepX, design = full.design, X.test = X.test, Y.test = Y.test)
table.full.design$BER

table.design = get.result(X.train = X, Y.train = Y, ncomp = ncomp, keepX = keepX, design = design, X.test = X.test, Y.test = Y.test)
table.design$BER

table.correlation.design.comp1 = get.result(X.train = X, Y.train = Y, ncomp = 4, keepX = keepX, design = correlation.design.comp1, X.test = X.test, Y.test = Y.test)
table.correlation.design.comp1$BER

table.correlation.design.allcomp = get.result(X.train = X, Y.train = Y, ncomp = 4, keepX = keepX, design = correlation.design.allcomp, X.test = X.test, Y.test = Y.test)
table.correlation.design.allcomp$BER










## weighting with the block.splsda correlation
comp = 1
x.xList <- list()
for(i in 1:length(X)){
    corDat <- rep(0, length(X))
    names(corDat) <- paste("cor", names(X)[i], names(X), sep = "_")
    for(j in (i):length(X)){
        corDat[j] <- as.numeric(cor(tcga$variates[[i]][,comp], tcga$variates[[j]][,comp]))
    }
    x.xList[[i]] <- corDat
}
corMat.diablo.comp1 <- do.call(rbind, x.xList)
rownames(corMat.diablo.comp1) <- colnames(corMat.diablo.comp1) <- c("mRNA", "miRNA", "CpGs", "Proteins")



comp = 2
x.xList <- list()
for(i in 1:length(tcga$variates)){
    corDat <- rep(0, length(tcga$variates))
    names(corDat) <- paste("cor", names(tcga$variates)[i], names(tcga$variates), sep = "_")
    for(j in (i):length(tcga$variates)){
        corDat[j] <- as.numeric(cor(tcga$variates[[i]][,comp], tcga$variates[[j]][,comp]))
    }
    x.xList[[i]] <- corDat
}


corMat.diablo.comp2 <- do.call(rbind, x.xList)
rownames(corMat.diablo.comp2) <- colnames(corMat.diablo.comp2) <- names(tcga$variates)

# correlation to the answer? how?

# explained variance? No because explained variance is not linked to discriminative power
# lapply(tcga$explained_variance,sum)



library(dplyr)
library(tidyr)



load("/Users/florian/Work/git/multi-group/trainTestDatasetsNormalized.RDATA")

Y <- droplevels(pam50Train0$Call)
names(Y) <- rownames(pam50Train0)
X <- list(mRNA = mrnaTrain0, miRNA = mirnaTrain0, CpGs = methTrain0, Proteins = protTrain0)
all(names(Y) == rownames(X[[1]]))
all(names(Y) == rownames(X[[2]]))
all(names(Y) == rownames(X[[3]]))
all(names(Y) == rownames(X[[4]]))
dim(X[[1]]); dim(X[[2]]); dim(X[[3]]); dim(X[[4]]);


###### From Amrit
ncomp <- 3
iter <- 10
folds <- 5
cpus <-  8

design <- matrix(c(0, 1, 1, 1,
1, 0, 0, 0,
1, 0, 0, 0,
1, 0, 0, 0), nrow = 4, ncol = 4)

design2 = matrix(c(0, 0, 0, 1,
0, 0, 0, 1,
0, 0, 0, 1,
1, 1, 1,0), nrow = 4, ncol = 4)

## select variables (set grid)
keep <- c(1, 2, 3, seq(5, 15, by = 5), 50, 100)
errList <- list()

for(i in 1 : length(keep)){
    print(paste("keepX",keep[i]))
    keepX = list(mRNA=rep(keep[i], 3), miRNA=rep(keep[i], 3), CpGs=rep(keep[i], 3), Proteins=rep(keep[i], 3))
    result = block.splsda(X = X, Y = Y, ncomp = ncomp,
    keepX = keepX, design = design,
    mode = "regression", bias = FALSE, init="svd",scheme="centroid", tol=1e-30, max.iter=500)
    #lapply(result$loadings, function(x)
    #apply(x, 2, function(i) names(i)[which(i != 0)]))
    
    #result2 = block.splsda(X = c(X[-1],X[1]), Y = Y, ncomp = ncomp,
    #keepX = keepX, design = design2,
    #mode = "regression", bias = FALSE, init="svd",scheme="centroid", tol=1e-30,max.iter=50000)
    #lapply(result2$loadings, function(x)
    #apply(x, 2, function(i) names(i)[which(i != 0)]))

    feat1 <- lapply(result$loadings, function(x)
    apply(x, 2, function(i) names(i)[which(i != 0)]))
    
    cvPerf = lapply(1 : iter, function(u){perf(result, validation = "Mfold", folds = folds, method.predict = "all")})
    ## Majority Vote
    cvPerf2 = lapply(1 : iter, function(x){unlist(cvPerf[[x]][names(cvPerf[[x]]) == "MajorityClass.error.rate"], recursive = FALSE)})
    cvPerf2 = lapply(cvPerf2, function(x){lapply(x, function(y) {t(y)})})
    cv.err <- lapply(1:3, function(x) {
        error <- as.data.frame(cbind(apply(simplify2array(lapply(cvPerf2, function(y) y[[x]])), c(1,2), mean),
        apply(simplify2array(lapply(cvPerf2, function(y) y[[x]])), c(1,2), sd)))
        colnames(error) <- paste(colnames(error), rep(c("mean", "sd"), each=(nlevels(Y)+2)), sep=".")
        return(error=error)
    })
    ErrSelected0 <- as.data.frame(do.call(rbind, cv.err))
    ErrSelected0$Method <- factor(rep(c("max.dist", "centroids.dist", "mahalanobis.dist"), e=max(ncomp)), levels=c("max.dist", "centroids.dist", "mahalanobis.dist"))
    ErrSelected0$Comp <- factor(rep(paste0("Comp1-", 1:max(ncomp)), 3), levels = paste0("Comp1-", 1:max(ncomp)))
    ErrSelected0$Group <- rep(1:max(ncomp), 3)
    colNames <- c("Basal", "Her2", "LumA", "LumB", "ER", "BER", "Method", "Comp", "Group")
    ErrSelected <- rbind(setNames(ErrSelected0[, c("Basal.mean", "Her2.mean", "LumA.mean", "LumB.mean", "Overall.ER.mean", "Overall.BER.mean", "Method", "Comp", "Group")], colNames),
    setNames(ErrSelected0[, c("Basal.sd", "Her2.sd", "LumA.sd", "LumB.sd", "Overall.ER.sd", "Overall.BER.sd", "Method", "Comp", "Group")], colNames))
    ErrSelected$Mean.SD <- rep(c("Mean", "SD"), each = nrow(ErrSelected0))
    
    diabloTrainErr <- ErrSelected %>% gather(ErrorRate, Value, Basal:BER) %>% tbl_df %>%
    spread(Mean.SD, Value)
    errList[[i]] <- diabloTrainErr
}


## plot error rates
errDat <- do.call(rbind, lapply(1:length(errList), function(i){
    mat <- subset(errList[[i]], Method == "centroids.dist" & ErrorRate == "BER")
    mat$keep <- keep[i]
    mat
}))

WhereAmI <- "~/Dropbox/Manuscript/mixOmics.org:DIABLO/"
source(paste0(WhereAmI, "functions/functions.R"))

ggplot(errDat, aes(x = keep, y = Mean, color = Comp)) + geom_point() +
geom_line() + scale_x_log10() +
geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD)) +
customTheme(sizeStripFont=15, xAngle=0, hjust=0.5, vjust=0.5,
xSize=10, ySize=10, xAxisSize=10, yAxisSize=10) +
xlab("Number of variables per component for each dataset") + ylab("Balanced Error Rate (BER, 10x5-fold cross-validation)") +
annotate("text", x = 50, y = 0.11, label = "p = 15 (5 per comp)") +
annotate("text", x = 50, y = 0.1, label = "BER = 13.2%")








