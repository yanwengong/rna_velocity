## install pagoda 2

## mac dependencies

# brew update
# brew install cmake boost eigen gsl curl openssl wget



# Install Bioconductor dependencies
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "org.Hs.eg.db","org.Mm.eg.db", "pcaMethods"), suppressUpdates=TRUE)
library(devtools)
install_github("igraph/rigraph") # Don't install with install.packages()
install_github("jkrijthe/Rtsne",ref="openmp")
install.packages(c("Cairo","urltools"))

# Install pagoda
install_github("hms-dbmi/pagoda2")
library('pagoda2')
# Pagoda2 is now ready to use