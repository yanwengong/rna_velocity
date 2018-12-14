## analyze mm10 data with 3 basals, 2 diff and 1 prolifering groups
## 10/18/2018
## analysis template on the mm10_data with 10x velocyto run 

## the loom file is prepared on HPC through velocyto

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
emat <- mm10_loom$spliced
nmat <- mm10_loom$unspliced
## restrict to cells that passed p2 filter
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]
# take cluster labels
cluster.label <- r$clusters$PCA$community # take the cluster factor that was calculated by p2 PCA
cell.colors <- pagoda2:::fac2col(cluster.label)
# take reduction form p2
 
## code2
emb_pca <- r$embeddings$PCA$largeVis

## code1
#emb_pca <- r$reductions$PCA ## Suoqin

## In addition to clustering and the t-SNE embedding, from the p2 processing we will also take a 
## cell-cell distance, which will be better than the default whole-transcriptome correlation 
##distance that velocyto.R would normally use.
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA))) ### calculate the cell dist use pagoda or Suoqin

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(nmat)))
fit.quantile_pca <- 0.05
rvel.cd_pca <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile_pca)
show.velocity.on.embedding.cor(emb_pca,rvel.cd_pca,n=500,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=0.7,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)


## Velocity estimation on t-SNE

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
emb_pca <- r$embeddings$PCA
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

### visualize velocity on PCAs 
pca.velocity.plot(rvel.cd_pca,nPcs=5,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1), max.grid.arrow.length = 0.0001, arrow.scale = 0.1)
pca.velocity.plot(rvel.cd_pca,nPcs=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1, max.grid.arrow.length = 0.0001, arrow.scale = 0.1)

# show.velocity.on.embedding.cor(emb_pca,rvel.cd,n=200,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)

#################### make the cell filtering less rigid - no filtration################

## prepare matrices and clustering data
emat <- mm10_loom$spliced
nmat <- mm10_loom$unspliced
## restrict to cells that passed p2 filter
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]
length(intersect(rownames(emat),rownames(nmat)))
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)
# filtered out 6847 out of 12077 genes due to low nmat-emat correlation
# filtered out 863 out of 5230 genes due to low nmat-emat slope
#saveRDS(rvel.cd, file = "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/rvel.cd.rds")

### visualize velocity on t-SNE embedding

show.velocity.on.embedding.cor(emb,rvel.cd,n=200,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)






## change the color label to SQ label


## convert df to factor with name attribute 
#fullseq_label_sub <- fullseq_label[1:5,]
for (i in 1:nrow(cell_list)){
  cell_i <- cell_list[i,1];
  group_i <- cell_list[i,2];
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


### visualize velocity on t-SNE embedding

show.velocity.on.embedding.cor(emb,rvel.cd,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.5),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)

### visualize velocity on PCA
show.velocity.on.embedding.cor(emb_pca,rvel.cd_pca,n=1000,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.5),cex=0.8,arrow.scale=0.7,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)

#pca.velocity.plot(rvel.cd_pca,nPcs=5,plot.cols=2,cell.colors=ac(cell.colors_fullseq,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1), max.grid.arrow.length = 0.0001, arrow.scale = 0.1)
#pca.velocity.plot(rvel.cd_pca,nPcs=2,cell.colors=ac(cell.colors_fullseq,alpha=0.7),cex=1.2,pcount=0.1, max.grid.arrow.length = 0.0001, arrow.scale = 0.1)


## plot the cell name with the same color 

## construct a table with cell group, count and color 
cell_type <- cell_list %>% group_by(V2) %>% summarise(count = n()) 
cell.colors_fullseq_df <- as.data.frame(cell.colors_fullseq)
cell_color <- cell.colors_fullseq_df %>% group_by(cell.colors_fullseq) %>%
  summarise(count = n())

cell_type_color<- cell_type %>%as.data.frame()%>%left_join(cell_color, by="count")
colors <-  unlist(cell_type_color$cell.colors_fullseq)

cell_type_color %>%
  ggplot(aes(x = V2, y = count, fill=V2)) + geom_bar(stat="identity") +
  scale_fill_manual(values =  as.character(unlist(cell_type_color$cell.colors_fullseq))) + geom_text(aes(label=count))+
  xlab("cell type") + ylab("counts")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12),
        axis.title.x  = element_text(size=12),
        legend.position = "none")





## Suoqing's PCA

pca_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/pca_info/ydata_scEpath.txt"
pca_dimension <- read.delim(pca_path)
cell_list$V1 <- gsub("CTL_BS_FullSeq_", "mm10_outs:", gsub(".1", "x", cell_list$V1))
