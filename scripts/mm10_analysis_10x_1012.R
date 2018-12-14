## analyze mm10 data
## 10/12/2018
## analysis template on the mm10_data with 10x velocyto run 

## the loom file is prepared on HPC through velocyto

library(velocyto.R)
library(pagoda2)
library(igraph)

## to avoid out of memory, reset . Renviron in terminal
## cd ~
## touch .Renviron
## open .Renviron
## R_MAX_VSIZE=50Gb

## set input loom file path
input_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/mm10_10x_1010.loom"

mm10_loom <- read.loom.matrices(input_path)
str(mm10_loom)


### use spliced expression matrix as input to pagoda2
emat <- mm10_loom$spliced
#### cell size here reflect the number of expressed spliced gene of individual cells 
hist(log10(colSums(emat)),col='wheat',xlab='cell size')

### filter to omit cells with small number of RNAs
emat <- emat[,colSums(emat)>1e3] ## may change to 200 ##this data set changed from 5372 to 5337
str(emat)

### filter to omit genes expressed in few cells
hist(log10(rowSums(emat)+1),col='wheat',xlab='cell size')
emat <- emat[rowSums(emat)>=10,] ## gene num changed from 27998 to 13082
str(emat)


## Pagoda2 processing

## make gene names unique by appending sequence numbers to duplicates
rownames(emat) <- make.unique(rownames((emat)))
str(emat)

## create pagoda2 object; NOT SURE what is plain, trim, log.scale
## why in this object the gene num changed from 13082 to 11186?
r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
## adjust the variance, to normalize the extent to which genes with very
## different expression will contribute to the downstream analysos
#### NOT SURE HOW TO UNDERSTAND THIS
r$adjustVariance(plot=T,do.par=T,gam.k=10)

## reduce the dataset dimentsions by running PCA
## nPcs number of PCs\n
## n.odgenes whether a certain number of top overdispersed genes should be used
## maxit maximum number of irlba iterations to use
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)

## make a knn graph
## k use k neighbors
## use cosine distance A*B/ (|A|*|B|)
r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')

## calculate clusters based on the KNN graph
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
#### object 'multilevel.community' not found
#####r$getKnnClusters(type='PCA')

###generate t-SNE embedding
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
###plot
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters (tSNE)')
## depth here is an expression pattern of gene'
## I am think it is the amount of expressed genes in single cell
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$depth,main='depth') 
dev.off()

## check one gene
gene <-"Rgs20"
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,gene],shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main=gene)



## Velocity estimation

## prepare matrices and clustering data
emat <- mm10_loom$spliced
nmat <- mm10_loom$unspliced
## restrict to cells that passed p2 filter
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]
# take cluster labels
cluster.label <- r$clusters$PCA$multilevel # take the cluster factor that was calculated by p2
cell.colors <- pagoda2:::fac2col(cluster.label)
# take embedding form p2
emb <- r$embeddings$PCA$tSNE

## In addition to clustering and the t-SNE embedding, from the p2 processing we will also take a 
## cell-cell distance, which will be better than the default whole-transcriptome correlation 
##distance that velocyto.R would normally use.
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

## filter genes by minimum average expression magnitude, output total number of resulting 
## valid genes
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(nmat)))

### estimate RNA velocity  (using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)
#saveRDS(rvel.cd, file = "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/rvel.cd.rds")

### visualize velocity on t-SNE embedding

show.velocity.on.embedding.cor(emb,rvel.cd,n=200,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)




## change the color label to SQ label

cluster.label <- r$clusters$PCA$multilevel # take the cluster factor that was calculated by p2
cell.colors <- pagoda2:::fac2col(cluster.label)

## read in the cell type file
fullseq_label <- read.delim("/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/tsne_info/identity_CTL_BS_FullSeq.txt",
                            header = FALSE, skip = 1)
## convert cell barcode to the correct format 
fullseq_label$V1 <- gsub("CTL_BS_FullSeq_", "mm10_outs:", gsub("-1", "x", fullseq_label$V1))

## convert df to factor with name attribute 
#fullseq_label_sub <- fullseq_label[1:5,]
for (i in 1:nrow(fullseq_label)){
  cell_i <- fullseq_label[i,1];
  group_i <- fullseq_label[i,2];
  names(group_i)<-cell_i
  
  if (i == 1){
    group <- group_i;
  } else {
    group <- c(group, group_i);
  }
}
group <- as.factor(group)

## assign color 
cell.colors_fullseq <- pagoda2:::fac2col(group)
