# mixOmics input data 

# created on 06-06-2016
# last modified: 08-06-2016
# Author: K-A Lê Cao
# purpose: reduce size of some data sets


library(mixOmicsv6)


# -------------------
# smaller data set: for package
# ------------------

# normalised data set
load('../Data/mixDIABLO/trainTestDatasetsNormalized.RDATA')

#             training data
# ----------------------------------
# remove lumB
remove.lumB.train = which(pam50Train0$Call == 'LumB')
set.seed(24)
# subset variables
keep.mrna  = sample(1:ncol(mrnaTrain0), 200, replace = FALSE)
##keep.methyl  = sample(1:ncol(methTrain0), 200, replace = FALSE)

# renaming training data
data.train = list(
  #clinical = clinTrain0[-c(remove.lumB.train),],
  #methylation = methTrain0[-c(remove.lumB.train), keep.methyl],
  mirna = mirnaTrain0[-c(remove.lumB.train),],
  mrna = mrnaTrain0[-c(remove.lumB.train), keep.mrna],
  protein = protTrain0[-c(remove.lumB.train),],
  subtype = droplevels(pam50Train0[-c(remove.lumB.train),]$Call)  # change this as subtype
)
lapply(data.train, dim)

# but still too big
# only keep a sample of the samples, LumA: 50%, Her2: 20%, Basal = 30%
set.seed(24)
keep.sample.her2 = sample(which(data.train$subtype == 'Her2'), round(150*0.2))
keep.sample.basal = sample(which(data.train$subtype == 'Basal'), round(150*0.3))
keep.sample.lumA = sample(which(data.train$subtype == 'LumA'), 150 - length(keep.sample.her2) - length(keep.sample.basal))
#keep.sample.train = sample(1:nrow(data.train$mirna), 150, replace = FALSE)
summary(data.train$subtype[c(keep.sample.basal, keep.sample.her2, keep.sample.lumA)])

keep.sample.train = c(keep.sample.basal, keep.sample.her2, keep.sample.lumA)


data.train.small = list(
  mirna = data.train$mirna[keep.sample.train,],
  mrna = data.train$mrna[keep.sample.train, ],
  protein = data.train$protein[keep.sample.train,],
  subtype = droplevels(data.train$subtype[keep.sample.train])  # change this as subtype
)
lapply(data.train.small, dim)
# check
summary(data.train.small$subtype)


#          test data
# ---------------------------
remove.lumB.test = which(pam50Test0$Call == 'LumB')

data.test = list(
  #clinical = clinTest0[-c(remove.lumB.test),],
  #methylation = methTest0[-c(remove.lumB.test), keep.methyl],
  mirna = mirnaTest0[-c(remove.lumB.test),],
  mrna = mrnaTest0[-c(remove.lumB.test), keep.mrna],
  protein = protAnnotation0[-c(remove.lumB.test),],
  subtype = droplevels(pam50Test0[-c(remove.lumB.test),]$Call)  # change this as subtype
)
lapply(data.test, dim)


# only keep a sample of the samples, LumA: 50%, Her2: 20%, Basal = 30%
set.seed(24)
keep.sample.her2 = sample(which(data.test$subtype == 'Her2'), round(70*0.2))
keep.sample.basal = sample(which(data.test$subtype == 'Basal'), round(70*0.3))
keep.sample.lumA = sample(which(data.test$subtype == 'LumA'), 70 - length(keep.sample.her2) - length(keep.sample.basal))
#keep.sample = sample(1:nrow(data.test$mirna), 70, replace = FALSE)
summary(data.test$subtype[c(keep.sample.basal, keep.sample.her2, keep.sample.lumA)])
keep.sample.test = c(keep.sample.basal, keep.sample.her2, keep.sample.lumA)


data.test.small = list(
  #clinical = data.test$clinical[keep.sample,],
  #methylation = data.test$methylation[keep.sample,],
  mirna = data.test$mirna[keep.sample.test,],
  mrna = data.test$mrna[keep.sample.test, ],
  #protein = data.test$protein[keep.sample,],
  subtype = data.test$subtype[keep.sample.test]
)
lapply(data.test.small, dim)

# check
summary(data.test.small$subtype)

data.train = data.train.small
data.test = data.test.small

breast.TCGA= list(data.train=data.train,data.test=data.test)
save(breast.TCGA, file="mixOmics/data/breast.TCGA.rda",compress="xz")


# -----------------------------------------
# erase then test
rm(list = ls())
load('mixOmics/data/breast.TCGA.rda')  # 
lapply(breast.TCGA$data.train, dim)
lapply(breast.TCGA$data.test, dim)



# -------------------------
# breast.tumours: for sPLS-DA with missing values
# ------------------------

data(breast.tumors)
X = breast.tumors$gene.exp
dim(X)

set.seed(8)
keep.gene = sample(1:ncol(X), 1000, replace = FALSE)
X2 = X[, keep.gene]
gene.name  = as.character(breast.tumors$genes$name[keep.gene])
description = as.character(breast.tumors$genes$description[keep.gene])


breast.tumors = list(gene.exp = X2, genes = list(name = gene.name, description), sample = breast.tumors$sample)
breast.tumors$genes

lapply(breast.tumors$genes, length)
lapply(breast.tumors, dim)

sum(is.na(breast.tumors$gene.exp))  # 2400

save(breast.tumors, file="package-mixomics/mixOmics/data/breast.tumors.rda",compress="xz")



# -------------------
# for mixMC
# -------------------

# ------- for multilevel
# read data in my directory
load('Data/mixMC/HMP-Diverse.RData')
dim(data.TSS)
dim(taxonomy)
dim(indiv)
bodysite = indiv$HMPbodysubsite
sample = indiv$RSID
indiv = indiv[, 1:5]
diverse.16S = list(data.TSS = data.TSS, data.raw = data.raw, taxonomy = taxonomy, indiv = indiv, bodysite=bodysite, sample=sample)
save(diverse.16S, file="package-mixomics/mixOmics/data/diverse.16S.rda",compress="xz")

# ----------- non multilevel
# read data in my directory
load('Data/mixMC/Koren.RData')
dim(data.TSS)
dim(taxonomy)
dim(indiv)
colnames(indiv[25:49])
indiv = indiv[, 25:46]
bodysite = indiv$sample_type

Koren.16S = list(data.TSS = data.TSS, data.raw = data.raw, taxonomy = taxonomy, indiv = indiv, bodysite=bodysite)
save(Koren.16S, file="package-mixomics/mixOmics/data/Koren.16S.rda",compress="xz")


# -------------------
# for liver.toxicity: rat ID name
# -------------------
data("liver.toxicity")
data = liver.toxicity$gene
rownames(data) = substr(rownames(liver.toxicity$gene), 2, 5)
rownames(data)  = paste('ID', rownames(data), sep = '')
data[1:5,1:5]

liver.toxicity$gene = data
lapply(liver.toxicity, dim)
save(liver.toxicity, file="package-mixomics/mixOmics/data/liver.toxicity.rda",compress="xz")

# -------------------
# for SRBCT
# -------------------
data(srbct)
srbct$gene[1:5,1:5]
srbct$gene.name[1:5,]

colnames(srbct$gene) = paste('g', 1:ncol(srbct$gene), sep ='')
rownames(srbct$gene.name) = colnames(srbct$gene)

save(srbct, file="package-mixomics/mixOmics/data/srbct.rda",compress="xz")
