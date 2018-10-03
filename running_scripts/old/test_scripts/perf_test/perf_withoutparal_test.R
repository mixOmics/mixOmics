## cheking results from Tenenhaus ##
####################################
rm(list = ls())
require(R.utils); require(mixOmics)
sourceDirectory("/Users/bgautier/bitbucket/mixomics/mixOmics/R/", modifiedOnly = FALSE)


#-- Section 9.7 --# 
#-----------------#
data("linnerud")
X <- linnerud$exercise
Y <- linnerud$physiological

linn.pls <- pls(X, Y, ncomp = 3, mode = "regression", tol = 1e-25)
linn.spls <- spls(X, Y, ncomp = 3, mode = "regression", tol = 1e-25)

## loadings and variates
linn.pls$loadings # pag. 142-143
linn.pls$variates # pag. 143-144

## prediction
indiv <- c(200, 40, 60)
pred <- predict(linn.pls, indiv)

pred$variates     # pag. 148
pred$predict[, , 2] # pag. 149

## Q2 criterion
res_pls <- perf(linn.pls, validation = "loo")
res_spls <- perf(linn.spls, validation = "loo")

all.equal(res_pls, res_spls)
res_pls$Q2 # pag. 142
res_spls$Q2 # pag. 142

#-- Section 7.2 --# 
#-----------------#
load("/Users/bgautier/bitbucket/perfwithoutparal_test/cornell.RData")
X <- cornell[, -8]
# Y <- cornell[, 8]
Y <- as.matrix(c(98.7, 97.8, 96.6, 92.0, 86.6, 91.2, 81.9, 83.1, 82.4, 83.2, 81.4, 88.1))

corn.pls <- pls(X, Y, ncomp = 4, mode = "regression", tol = 1e-25)
corn.spls <- spls(X, Y, ncomp = 4, mode = "regression")

## variates
cbind(corn.pls$variates$X, comp_Y = corn.pls$variates$Y[, 1])  # pag. 84
apply(cbind(corn.pls$variates$X, comp_Y = corn.pls$variates$Y[, 1]), 2, mean)  # pag. 84
apply(cbind(corn.pls$variates$X, comp_Y = corn.pls$variates$Y[, 1]), 2, var)  # pag. 84

## prediction
pred <- predict(corn.pls, X)
pred$predict[[3]] # pag. 85

## Q2 criterion
res_pls <- perf.pls(corn.pls, validation = "loo")
res_spls <- perf.spls(corn.spls, validation = "loo")

all.equal(res_pls, res_spls)
res_pls$Q2 # pag. 83
res_spls$Q2 # pag. 83

#-- Section 10.3 --# 
#------------------#
the <- read.csv("/Users/bgautier/bitbucket/perfwithoutparal_test/the_data.csv")

X <- cbind(unmap(the[, 1]), unmap(the[, 2]), unmap(the[, 3]), unmap(the[, 4]))
colnames(X) <- c("chaud", "ti?de", "froid", 
                 "sucre0", "sucre1", "sucre2", 
                 "fort", "moyen", "l?ger",
                 "citron1", "citron0")

Y <- the[, 5:10]
Y = 19-Y

the.pls <- pls(X, Y, ncomp = 5, mode = "regression", tol = 1e-25)
the.spls <- spls(X, Y, ncomp = 5, mode = "regression", tol = 1e-25)

for (i in 1 : 4){
  print(mean(sapply(1 : 6, function(x){cor(the.pls$variates$X[, i], Y[, x])^2})))
}

for (i in 1 : 4){
  print(mean(sapply(1 : 11, function(x){cor(the.pls$variates$X[, i], X[, x])^2})))
}

## Q2 criterion
res_pls <- perf.pls(the.pls, validation = "loo")
res_spls <- perf.spls(the.spls, validation = "loo")

all.equal(res_pls, res_spls)

res_pls$Q2 # pag. 179
res_spls$Q2 # pag. 179

plot(the.pls$mat.c[, 1], the.pls$mat.c[, 2], cex = 0.1, xlim = c(-0.7, 0.8), ylim = c(-0.6, 0.7))
text(the.pls$mat.c[, 1], the.pls$mat.c[, 2], labels = row.names(the.pls$mat.c))
text(the.pls$loadings$Y[, 1], the.pls$loadings$Y[, 2], labels = row.names(the.pls$loadings$Y))
for (i in 1 : 6){
  for (j in 1 : 10){
    arrows(0,0,the.pls$loadings$Y[i, 1], the.pls$loadings$Y[i, 2], length = 0.2, angle = j)
  }
}
abline(h = 0, v = 0, col = "red")

pred = predict(the.pls, newdata = X)
plot(pred$predict[, , 4][, 6], Y[, 6], cex = 0.1)
text(pred$predict[, , 4][, 6], Y[, 6], label = c(1:18))

plot(pred$predict[, , 4][, 4], Y[, 4], cex = 0.1)
text(pred$predict[, , 4][, 4], Y[, 4], label = c(1:18))

#-- Processionnaire du pin --#
#----------------------------#
proc.pin = read.csv("/Users/bgautier/bitbucket/perfwithoutparal_test/proc.pin.csv", row.names = 1)
proc.pin2 = proc.pin[1 : 33,]
proc.pin2[, 11] = log(proc.pin2[, 11])

cor(proc.pin2) 

proc.pin3 = proc.pin[-(1 : 33),]
proc.pin3[, 11] = log(proc.pin3[, 11])

cor(proc.pin3)

X = proc.pin2[, -11]
Y = proc.pin2[, 11, drop = FALSE]

pin.pls <- pls(X, Y, ncomp = 10, mode = "regression")
res_pls <- perf.pls(pin.pls, validation = "loo")
res_pls$Q2

