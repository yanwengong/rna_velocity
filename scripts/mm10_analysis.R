## analyze mm10 data

## the loom file is prepared on HPC through velocyto

library(velocyto.R)
library(pagoda2)
library(igraph)

## to avoid out of memory, reset . Renviron in terminal
## cd ~
## touch .Renviron
## open .Renviron
## R_MAX_VSIZE=50Gb
mm10_loom <- read.loom.matrices("/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/possorted_genome_bam_E7XYF.loom")
str(mm10_loom)

### use spliced expression matrix as input to pagoda2
emat <- mm10_loom$spliced
#### cell size here reflect the number of expressed spliced gene of individual cells 
hist(log10(colSums(emat)),col='wheat',xlab='cell size')

### filter to omit cells with small number of RNAs
emat <- emat[,colSums(emat)>1e3] ## may change to 200
str(emat)

### filter to omit genes expressed in few cells
hist(log10(rowSums(emat)+1),col='wheat',xlab='cell size')
emat <- emat[rowSums(emat)>=10,] ##
str(emat)

## Pagoda2 processing

### make gene names unique by appending sequence numbers to duplicates
rownames(emat) <- make.unique(rownames((emat)))
str(emat)
### create pagoda2 object
r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
### adjust the variance, to normalize the extent to which genes with very
### different expression will contribute to the downstream analysos
#### NOT SURE HOW TO UNDERSTAND THIS
r$adjustVariance(plot=T,do.par=T,gam.k=10)

### reduce the dataset dimentsions by running PCA
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')


r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
#### object 'multilevel.community' not found
#####r$getKnnClusters(type='PCA')


###generate t-SNE embeddingr4
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
###plot
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$depth,main='depth') 


## RNA velocity estimation 

### prepare matrices and clustering data
emat <- mm10_loom$spliced
nmat <- mm10_loom$unspliced
#### restrict to cells that passed p2 filter
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]

### take cluster labels
#### NOTE, i did community here because i don't have multilevel community
cluster.label <- r$clusters$PCA$multilevel # take the cluster factor that was calculated by p2
cell.colors <- pagoda2:::fac2col(cluster.label)
### take embedding form previous
emb <- r$embeddings$PCA$tSNE
### also take a cell-cell distance 
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

### filter genes by minimum average expression magnitude, output total number of resulting valid genes
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(nmat)))

### estimate RNA velocity  (using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)
#saveRDS(rvel.cd, file = "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/rvel.cd.rds")

### visualize velocity on t-SNE embedding

show.velocity.on.embedding.cor(emb,rvel.cd,n=200,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
