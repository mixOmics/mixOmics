---
title: "sPCA"
author: "KA Le Cao, Sebastien Dejean, Xin-Yi Chua, Danielle Davenport"
date: "19 January 2017"
output: html_document
---



# Principal Component Analysis (PCA)

Principal Component Analysis (Jolliffe, 2005) is primarily used to explore one single type of ‘omics data (e.g. transcriptomics, proteomics, metabolomics, etc) and identify the largest sources of variation. PCA is a mathematical procedure that uses orthogonal linear transformation of data from possibly correlated variables into uncorrelated principal components (PCs). The first principal component explains as much of the variability in the data as possible, and each following PC explains as much of the remaining variability as possible. Only the PCs which explain the most variance are retained. This is why choosing the number of dimensions or components **(ncomp)** is crucial (see the function **tune.pca**, below).

In **mixOmics**, PCA is numerically solved in two ways:

**1.** With singular value decomposition (SVD) of the data matrix,which is the most computationally efficient way and is also adopted by most softwares and the R function *prcomp* in the stat package. 

**2.** With the Non-linear Iterative Partial Least Squares (NIPALS) in the case of missing values, which uses an iterative power method. See [Methods: Missing values](http://mixomics.org/methods/missing-values/). 

Both methods are embedded in the **mixOmics** *pca* function and will be used accordingly. 

Input data should be centered *(center = TRUE)* and possibly (sometimes preferably) scaled so that all variables have a unit variance. This is especially advised in the case where the variance is not homogeneous across variables *(scale = TRUE)*. By default, the variables are centered and scaled in the function, but the user is free to choose other options.


```r
library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene # Using one data set only
```

## Choosing the optimal parameters

We can obtain as many dimensions (i.e. number of PCs) as the minimum between the number of samples and variables. However, the goal is to reduce the complexity of the data and therefore summarize the data in fewer underlying dimension. 

The number of principal Components to retain (also called the number of dimensions) is therefore crucial when performing PCA. The function **tune.pca** will plot the barplot of the proportion of explained variance for min(*n*, *p*)principal components, where *n* is the number of samples, and *p* the number of variables. 







