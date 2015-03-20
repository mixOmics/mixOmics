u <- matrix(runif(1000 * 1000, min = 0, max = 10), nrow = 1000, ncol = 1000, 
            dimnames = list(paste0("ind", 1:1000), paste0("col", 1 : 1000)))

mymean = function(x){
  temp = 0
  for (i in 1 : length(x)){temp = temp + x[i]}
  return(temp/length(x))
}

source("paral.options.R")
library(parallel)
cluster = paral.options(parallel = makeCluster(4, type = "PSOCK"))$parallel

### Comparison CPU time on Apply Functions using/or not Clusters
### apply vs parApply
test1 = apply(u, 2, mymean)
system.time(apply(u, 2, mymean))

test1.p = applyFunc(cluster, u, 2, mymean)
system.time(applyFunc(cluster, u, 2, mymean))

all.equal(test1, test1.p)

### sapply vs parSapply
test2 = sapply(1 : ncol(u), function(x){mymean(u[, x])})
system.time(sapply(1 : ncol(u), function(x){mymean(u[, x])}))

clusterExport(cluster, c("mymean", "u"))
test2.p = sapplyFunc(cluster, 1: ncol(u), function(x){mymean(u[, x])})
system.time(sapplyFunc(cluster, 1: ncol(u), function(x){mymean(u[, x])}))

all.equal(test2, test2.p)

### lapply vs parLapply
test3 = lapply(1 : ncol(u), function(x){mymean(u[, x])})
system.time(lapply(1 : ncol(u), function(x){mymean(u[, x])}))

clusterExport(cluster, c("mymean", "u"))
test3.p = lapplyFunc(cluster, 1: ncol(u), function(x){mymean(u[, x])})
system.time(lapplyFunc(cluster, 1: ncol(u), function(x){mymean(u[, x])}))

all.equal(test3, test3.p)

stopCluster(paral.options()$parallel)
