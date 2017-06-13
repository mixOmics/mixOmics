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
#         mixDIABLO tab: first trash the existing files!
# ======================================
# intro
knit2wp(input = 'mixDIABLO.Rmd', title = 'mixDIABLO',  shortcode = FALSE ,publish = FALSE, action = c("newPage"))

# analysis example with TCGA: make sure the parent is from mixDIABLO!
knit2wp(input = 'mixDIABLO_TCGA_example.Rmd', title = 'Example of mixDIABLO analysis',  shortcode = FALSE ,publish = FALSE, action = c("newPage"))
