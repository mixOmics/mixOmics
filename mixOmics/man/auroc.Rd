\name{auroc}
\encoding{latin1}
\alias{auroc}
\alias{auroc.plsda}
\alias{auroc.splsda}
\alias{auroc.mint.plsda}
\alias{auroc.mint.splsda}
\alias{auroc.sgccda}


\title{Area Under the Curve (AUC) and Receiver Operating Characteristic (ROC) curves for supervised classification}

\description{
Calculates the AUC and plots ROC for supervised objects from s/plsda, mint.s/plsda and block.spslda or wrapper.sgccda.
}

\usage{
auroc.plsda(object, newdata = object$X, outcome.test = as.factor(object$Y), plot = TRUE, ncomp = 1, ...)
auroc.splsda(object, newdata = object$X, outcome.test = as.factor(object$Y), plot = TRUE, ncomp = 1, ...)
auroc.mint.plsda(object, newdata = object$X, outcome.test = as.factor(object$Y), study.test = object$study, plot = TRUE, ncomp = 1, ...)
auroc.mint.splsda(object, newdata = object$X, outcome.test = as.factor(object$Y), study.test = object$study, plot = TRUE, ncomp = 1, ...)
auroc.sgccda(object, newdata = object$X, outcome.test = as.factor(object$Y), plot = TRUE, block = 1, ncomp = 1, ...)

}	

\arguments{
\item{object}{Object from one of the following supervised analysis class: "plsda", "splsda", "mint.plsda", "mint.splsda", "block.splsda" or "wrapper.sgccda"}
\item{newdata}{numeric matrix of predictors, by default set to the training data set (see details).}
\item{outcome.test}{Either a factor or a class vector for the discrete outcome, by default set to the outcome vector from the training set (see details). }
\item{plot}{Whether the ROC curves should be plotted, by default set to TRUE (see details).}
\item{ncomp}{ROC plots will be output for that specific component from the multivariate model, default to 1.}
\item{block}{For a block.plsda or block.splsda object, either the block number (integer) or the name of the block (set of characters) to specify which ROC plot to be output, default to 1.}
\item{...}{external optional arguments for plotting}
}

\details{
For more than two classes in the categorical outcome Y, the AUC is calculated as one class vs. the other and the ROC curves one class vs. the others are output.

The ROC and AUC are calculated based on the predicted scores obtained from the multivariate methods. Be mindful that our multivariate supervised methods already use a prediction threshold based on distances (see \code{predict}) that optimally determine class membership of the samples tested. As such AUC and ROC are not needed to estimate the performance of the model (see \code{perf}, \code{tune} that report classification error rates). We provide those outputs as complementary performance measures. 

The pvalue is from a Wilcoxon test between the predicted scores between one class vs the others. 

External independent data set could be input with the same number of predictors than in the training data set. In that case the multivariate model would predict on the test data (based on the trained model stored in object) and compare the prediction to the avaible and known outcome vector from the external independent data set
}

\value{
Depending on the type of object used, a list that contains:
The AUC and Wilcoxon test pvalue for each 'one vs other' classes comparison performed, either per component (splsda, plsda, mint.plsda, mint.splsda), or per block and per component (wrapper.sgccda, block.plsda, blocksplsda).
}


\author{Benoit Gautier, Francois Bartolo, Florian Rohart}

\seealso{\code{\link{tune}}, \code{\link{perf}}, and http://www.mixOmics.org for more details.}

\examples{
## First example with PLSDA, 2 classes
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment

plsda.breast <- plsda(X, Y, ncomp = 2)
auc.plsda.breast = auroc(plsda.breast, ncomp = 1)

## second example with sPLSDA
splsda.breast <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))
auroc(plsda.breast, plot = FALSE)

## third example with sPLSDA with 4 classes
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))
auc.splsda.liver = auroc(splsda.liver, ncomp = 1)

## fourth example with block.plsda
data(nutrimouse)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
# with this design, all blocks are connected
design = matrix(c(0,1,1,0), ncol = 2, nrow = 2,
byrow = TRUE, dimnames = list(names(data), names(data)))

block.plsda.nutri = block.plsda(X = data, Y = nutrimouse$diet)
auc.block.plsda.nutri = auroc(block.plsda.nutri, block = 'lipid')

## fith example with block.splsda
list.keepX = list(gene = rep(10, 2), lipid = rep(5,2))
block.splsda.nutri = block.splsda(X = data, Y = nutrimouse$diet, keepX = list.keepX)
auc.block.splsda.nutri = auroc(block.splsda.nutri, block = 1)


}

\keyword{regression}
\keyword{multivariate}
