library(velocyto.R)
library(pagoda2)
library(igraph)
library(tidyverse)

## to avoid out of memory, reset . Renviron in terminal
## cd ~
## touch .Renviron
## open .Renviron
## R_MAX_VSIZE=50Gb

## read in the loom file 
input_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/mm10_10x_1010.loom"
cell_list_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/cell_groups/identity_CTL_BS_Krt14_Krt1_Prolif_k6.txt"
cell_list_path0 <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/tsne_info/identity_CTL_BS_FullSeq.txt"
mm10_loom <- read.loom.matrices(input_path)

## read in the cell list and convert cell barcode to the correct format 
cell_list <- read.delim(cell_list_path, header = FALSE, skip = 1)
cell_list$V1 <- gsub("CTL_BS_FullSeq_", "mm10_outs:", gsub(".1", "x", cell_list$V1))

## filter the spliced loom file with cell_list 
emat <- mm10_loom$spliced %>% as.matrix()%>%as.data.frame()%>%select(one_of(cell_list$V1))

## chekc cell size, shown at least 1e3
## cell size here reflect the number of expressed spliced gene of individual cells 
hist(log10(colSums(emat)),col='wheat',xlab='cell size')

## check gene expression 
hist(log10(rowSums(emat)+1),col='wheat',xlab='cell size')
### THIS MAY BE MODIFIED
emat <- emat[rowSums(emat)>=10,] ## gene num changed from 27998 to 11059
## 11059 obs. of  3021 variables

## Pagoda2 processing

## make gene names unique by appending sequence numbers to duplicates
rownames(emat) <- make.unique(rownames((emat)))

## convert back to dgCMatrix
emat <- emat %>% as.matrix() %>% as("dgCMatrix")

## create pagoda2 object; NOT SURE what is plain, trim, log.scale
r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
## adjust the variance, to normalize the extent to which genes with very
## different expression will contribute to the downstream analysos
#### NOT SURE HOW TO UNDERSTAND THIS
r$adjustVariance(plot=T,do.par=T,gam.k=10)

## reduce the dataset dimentsions by running PCA
## nPcs number of PCs\n
## n.odgenes whether a certain number of top overdispersed genes should be used
## maxit maximum number of irlba iterations to use
set.seed(0)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)

## make a knn graph
## k use k neighbors
## use cosine distance A*B/ (|A|*|B|)
set.seed(1)
r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine') ## optional for PCA?



## calculate clusters based on the KNN graph for PCA
r$getKnnClusters(method=infomap.community,type='PCA') ## optional for PCA
M <- 30; 
r$getEmbedding(type='PCA',M=M,perplexity=30,gamma=1/M,alpha=1) ## optional for PCA

## plot PCA
set.seed(1)
r$plotEmbedding(type='PCA',show.legend=F,
                mark.clusters=T, min.group.size=10, shuffle.colors=F, 
                mark.cluster.cex=1, alpha=0.3, main='PCA cluster')

#### object 'multilevel.community' not found
#####r$getKnnClusters(type='PCA')

## calculate clusters based on the KNN graph for t-sne
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')

###generate t-SNE embedding
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
###plot
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,
                mark.clusters=T,min.group.size=10,shuffle.colors=F,
                mark.cluster.cex=1,alpha=0.3,main='cell clusters (tSNE)')
## depth here is an expression pattern of gene'
## I am think it is the amount of expressed genes in single cell
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$depth,main='depth') 
dev.off()


## Velocity estimation on PCA
## prepare matrices and clustering data
emat <- mm10_loom$spliced #27998x5372
nmat <- mm10_loom$unspliced #27998x5372
## restrict to cells that passed p2 filter
emat <- emat[,rownames(r$counts)]; # 27998 3021
nmat <- nmat[,rownames(r$counts)]
# take cluster labels # this can be done after makeKnnGraph and getKnnClusters
cluster.label <- r$clusters$PCA$community # take the cluster factor that was calculated by p2 PCA

cell.colors <- pagoda2:::fac2col(cluster.label)
# take reduction form p2
 emb_pca <- r$embeddings$PCA$largeVis
#emb_pca <- r$reductions$PCA ## Suoqin

## In addition to clustering and the t-SNE embedding, from the p2 processing we will also take a 
## cell-cell distance, which will be better than the default whole-transcriptome correlation 
##distance that velocyto.R would normally use.
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA))) ### calculate the cell dist use pagoda or Suoqin

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.2) ## dramatically decreased:4726x3021
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)## dramatically decreased:5933x3021
length(intersect(rownames(emat),rownames(nmat)))
fit.quantile_pca <- 0.02
rvel.cd_pca <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile_pca)
show.velocity.on.embedding.cor(emb_pca,rvel.cd_pca,n=500,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=0.7,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)

