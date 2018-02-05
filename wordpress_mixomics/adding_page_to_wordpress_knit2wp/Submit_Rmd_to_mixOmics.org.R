############################
# submit Rmd files to the mixOmics.org website
# Kim-Anh Le Cao
# last updatedL 24 May 2017
# ##########################



# set up (see the BloggingfromR.Rmd file)
library(RWordPress)
library(knitr)
# Set figure dimensions
#opts_chunk$set(fig.width=5, fig.height=5)
# Set figures to upload to imgur.com
opts_knit$set(upload.fun = imgur_upload, base.url = NULL) 

# Set your WP username, password, and your site URL
options(WordpressLogin = c(your.username = 'your.password'), #update these
        WordpressURL = 'http://mixomics.org/xmlrpc.php')








# ======================================
#         Methods tab: first trash the existing files!
# ======================================

knit2wp(input = 'PCA.Rmd', title = 'PCA',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

knit2wp(input = 'sIPCA.Rmd', title = 'IPCA',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

knit2wp(input = 'rCCA.Rmd', title = 'rCCA',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

knit2wp(input = 'sPLS.Rmd', title = 'sPLS',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

knit2wp(input = 'PLS-DA.Rmd', title = 'PLS-DA',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

knit2wp(input = 'Multilevel.Rmd', title = 'Multilevel',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))


# ======================================
#         Graphics tab: first trash the existing files!
# ======================================

knit2wp(input = 'plotIndiv.Rmd', title = 'plotIndiv',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))



# ======================================
#         Case studies tab: first trash the existing files!
# ======================================

knit2wp(input = 'PLSDA_SRBCT.Rmd', title = 'sPLSDA: SRBCT',  shortcode = FALSE ,publish = FALSE, 
        action = c("newPage"))


# ======================================
#         DIABLO tab: first trash the existing files!
# ======================================
# intro
# somehow struggles to get the link diablo instead of diablo-2
# check in the appearance ->menu first
knit2wp(input = 'Diablo.Rmd', title = 'DIABLO',  shortcode = FALSE ,publish = FALSE, action = c("newPage"))

# analysis example with TCGA: 
# !!! set to 'tcga-example in the quick edit slug; and set no parents, this is important for our publications links
knit2wp(input = 'DIABLO_TCGA.Rmd', title = 'Case study: TCGA',  shortcode = FALSE ,publish = FALSE, 
        action = c("newPage"))

# ======================================
#         MINT tab: first trash the existing files!
# ======================================

# analysis example with stem cells
# !!!! set to 'stemcells-example in the quick edit slug; and set no parents, this is important for our publications links
knit2wp(input = 'MINT_stemcells.Rmd', title = 'Case study: stem cells',  shortcode = FALSE ,publish = FALSE, 
        action = c("newPage"))



# ======================================
#         mixMC tab: first trash the existing files!
# ======================================

knit2wp(input = 'mixMC_parent.Rmd', title = 'mixMC',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

knit2wp(input = 'mixMC_Normalisation.Rmd', title = 'Pre-processing',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

knit2wp(input = 'mixMC_nonmultilevel_Koren.Rmd', title = 'Case study: Koren diverse bodysites',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

knit2wp(input = 'mixMC_Multilevel_HMP_16S_Data.Rmd', title = 'Case study: HMP bodysites repeated measures',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

# make sure the link is mixOmics.org/mixkernel!
knit2wp(input = 'mixKernelUsersGuide.Rmd', title = 'Case study: Tara ocean with mixKernel',  shortcode = FALSE, publish = FALSE, 
        action = c("newPage"))

