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

options(WordpressLogin = c(klecao = 'mixOmics8'), #update these
        WordpressURL = 'http://mixomics.org/xmlrpc.php')





# ======================================
#         mixDIABLO tab: first trash the existing files!
# ======================================
# intro
knit2wp(input = 'DIABLO.Rmd', title = 'DIABLO',  shortcode = FALSE ,publish = FALSE, action = c("newPage"))

# analysis example with TCGA: make sure the parent is from mixDIABLO!
knit2wp(input = 'DIABLO_TCGA.Rmd', title = 'TCGA example',  shortcode = FALSE ,publish = FALSE, 
        action = c("newPage"))



# ======================================
#         MINT tab: first trash the existing files!
# ======================================

# analysis example with stem cells: make sure the parent is from MINT!
knit2wp(input = 'MINT_stemcells.Rmd', title = 'Stemcells example',  shortcode = FALSE ,publish = FALSE, 
        action = c("newPage"))



# ======================================
#         Case studies tab: first trash the existing files!
# ======================================

# analysis example with stem cells: make sure the parent is from MINT!
knit2wp(input = 'PLSDA_SRBCT.Rmd', title = 'SRBCT example',  shortcode = FALSE ,publish = FALSE, 
        action = c("newPage"))
