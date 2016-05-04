# ----------------------------------------
# compile for colors
# date: 16/03/2015
# latest update: 
# ----------------------------------------


# ------ notes for me to compile the package (if need be)
# need to be root
R CMD build --resave-data mixOmics
R CMD INSTALL -l MyR/ mixOmics_5.0-4.tar.gz 
R CMD check mixOmics --as-cran --timings
# --------------------------------

sessionInfo()