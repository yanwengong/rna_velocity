# ## terminal
# brew install hdf5 --enable-cxx
# brew install cmake boost eigen gsl curl openssl wget
# download clang4-r.pkg
# download gfortran-6.1-ElCapitan.dmg
# ##

install.packages("dplyr")
install.packages("ggplot2")
install.packages("devtools")

source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
install_github("mannau/h5")


library(devtools)  
install_github("velocyto-team/velocyto.R")


library(velocyto.R)
#wget http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom)
ldat <- read.loom.matrices("/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/SCG71.loom")
str(ldat)


# Using spliced expression matrix as input to pagoda2

emat <- ldat$spliced
hist(log10(colSums(emat)),col='wheat',xlab='cell size')

# filter ( this dataset already been pre-filtered)

emat <- emat[,colSums(emat)>=1e3]









# reduce the cell names to the short well labels
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("_unique.bam","",gsub(".*:","",colnames(x)))
  x
})


# Spliced expression magnitube dist across gene
hist(log10(rowSums(ldat$spliced)+1), col = 'wheat', xlab = 'log10[ number of reads + 1]',main='number of reads per gene')

# exonic read (spliced) expression matrix
emat <- ldat$spliced;
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$ambiguous;
# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.5)
# look at the resulting gene set
length(intersect(rownames(emat),rownames(nmat)))



