#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

## Hilbert matrix
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
X.na <- X <- hilbert(9)[, 1:6]

## Hilbert matrix with missing data
idx.na <- matrix(sample(c(0, 1, 1, 1, 1), 36, replace = TRUE), ncol = 6)
X.na[idx.na == 0] <- NA
X.rec <- nipals(X.na, reconst = TRUE)$rec
round(X, 2)
round(X.rec, 2)


#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    
    
    X <- matrix(sample(c(0, 1, 1, 1, 1), 36, replace = TRUE), ncol = 3)
    
    #test ncomp
    nipals(X,ncomp=1)
    nipals(X,ncomp=2)
    nipals(X,ncomp=3)
    
    #test reconst
    nipals(X,reconst=TRUE)$rec
    nipals(X,reconst=FALSE)$rec
    
    #test max.iter
    nipals(X,max.iter=10,ncomp=2)$eig
    nipals(X,max.iter=20,ncomp=3)$eig
    nipals(X,max.iter=1000)$eig
    
    #test tol
    nipals(X,tol=1e-1)$eig
    nipals(X,tol=1e-2)$eig
    nipals(X,tol=1e-3)$eig
    
}
par(opar)
