# ----------------------------------------
# Compile for devel version
# date: 16/02/2015
# latest update
# ----------------------------------------


# ------ notes for me to compile the package (if need be)
R CMD build --resave-data mixOmics
R CMD INSTALL -l MyR/ mixOmics_5.0-4.tar.gz 
R CMD check mixOmics --as-cran --timings
# --------------------------------

# list of authors in DESCRIPTION file
Author: Sebastien Dejean, Ignacio Gonzalez, Kim-Anh Le Cao with
contributions from Pierre Monget, Jeff Coquery, FangZou Yao,
Benoit Liquet, Florian Rohart, Benoit Gautier

# woud need to put as a Author@R but did not work
# see http://r-pkgs.had.co.nz/description.html help!)
Authors@R: c(
  person("Kim-Anh", "Le Cao", email = "k.lecao@uq.edu.au", role = "cre"),
  person("Ignacio", "Gonzalez", email = "igonzalez@toulouse.inra.fr", role = "aut"),
  person("Sébastien", "Déjean", email = "sebastien.dejean@math.univ-toulouse.fr", role = "aut"),
  person("Florian", "Rohart", email = "k.lecao@uq.edu.au", role = "aut"),
  person("Benoit", "Gautier", email = "k.lecao@uq.edu.au", role = "aut"),
  person("Benoit", "Liquet", email = "b.liquet@uq.edu.au", role = "ctb", comment = 'contribution to the multilevel module'),
  person("FangZou", "Yao", role = "ctb", comment = 'contribution to the IPCA module')
  person("Pierre", "Monget", role = "ctb", comment = 'contribution to the sPLS-DA module')
  person("Jeff", "Coquery", role = "ctb", comment = 'contribution to the IPCA module')
)


sessionInfo()
# R version 3.1.0 (2014-04-10)
# Platform: x86_64-apple-darwin10.8.0 (64-bit)
# 
# locale:
#   [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] tools_3.1.0

