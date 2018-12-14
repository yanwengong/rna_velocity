## analyze 2 unwound(BS) samples separately 
## 12/04/2018
library(velocyto.R)
library(pagoda2)
library(igraph)
library(tidyverse)
## FUN generate pagoda object; conduct dimension reduction

## this is function that build pagoda object from 
## input: the filtered spliced data
## output: the pagoda object r with pca reduction 

build_pagoda <- function (input){
  ## make gene names unique by appending sequence numbers to duplicates
  rownames(input) <- make.unique(rownames((input)))
  
  ## NOTE: i need to do this step here not in filter_cell_gene function to avoid error in velocity calculation 
  input <- input %>% as.matrix() %>% as("dgCMatrix")
  
  ## create pagoda2 object; NOT SURE what is plain, trim, log.scale
  r <- Pagoda2$new(input,modelType='plain',trim=10,log.scale=T)
  ## adjust the variance, to normalize the extent to which genes with very
  ## different expression will contribute to the downstream analysos
  #### NOT SURE HOW TO UNDERSTAND THIS
  r$adjustVariance(plot=T,do.par=T,gam.k=10)
  
  ## reduce the dataset dimentsions by running PCA
  ## nPcs number of PCs\n
  ## n.odgenes whether a certain number of top overdispersed genes should be used
  ## maxit maximum number of irlba iterations to use
  ## used r package irlba
  set.seed(0)
  r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
  
  ## make a knn graph
  ## k use k neighbors
  ## use cosine distance A*B/ (|A|*|B|)
  set.seed(1)
  r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine') ## optional for PCA?
  
  ## calculate clusters based on the KNN graph for PCA
  r$getKnnClusters(method=infomap.community,type='PCA') ## optional for PCA
  
  return(r)
  
}

## assign cell group and color
## FUN generate cluster label color from cell list 
## this function will take cell list as input and generate color assignment for each cell
## input: cell list
## output: cluster info for each cell; color assignment for each cell

assign_color <- function(cell_list){
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
  return (list(group, cell.colors_fullseq))
}

## read in cell group and color file
file_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/mm10_mouse_skin_loom_all_1201"
cellColors_bs <- read.delim(paste(file_path, "niedaigroupmeeting/colorCells_bs_epi_cca_basal_spinious_prolif.txt", sep="/"), header = FALSE)
identity_bs <- read.delim(paste(file_path, "niedaigroupmeeting/identity_bs_epi_cca_basal_spinious_prolif.txt", sep="/"), header = TRUE) 
identity_bs$Row <- gsub("..$", "x", identity_bs$Row)
identity_bs$Row <- gsub("bs_1_", "mm10_outs:", gsub("bs_2_", "dh_wt_ctrl_ms_bs_2:", identity_bs$Row))
identity_bs_1 <- identity_bs%>%filter(str_detect(Row, "mm10_outs:"))
identity_bs_2 <- identity_bs%>%filter(str_detect(Row, "dh_wt_ctrl_ms_bs_2:"))

## read in loom files and concatenate
## bs
dh_wt_ctrl_ms_bs_2_path <- paste(file_path, "dh_wt_ctrl_ms_bs_2.loom", sep="/")
dh_wt_ctrl_ms_bs_1_path <- paste(file_path, "dh_wt_ctrl_ms_bs_1.loom", sep="/")
dh_wt_ctrl_ms_bs_2 <- read.loom.matrices(dh_wt_ctrl_ms_bs_2_path)
dh_wt_ctrl_ms_bs_1 <- read.loom.matrices(dh_wt_ctrl_ms_bs_1_path)

## filter emat based on cell_list
emat_bs_1_filtered <- dh_wt_ctrl_ms_bs_1$spliced%>%as.matrix()%>%as.data.frame()%>%select(one_of(identity_bs$Row))
emat_bs_2_filtered <- dh_wt_ctrl_ms_bs_2$spliced%>%as.matrix()%>%as.data.frame()%>%select(one_of(identity_bs$Row))

r_bs_1 <- build_pagoda(emat_bs_1_filtered)
r_bs_2 <- build_pagoda(emat_bs_2_filtered)


## plot tSNE
r_bs_1$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
set.seed(1)
r_bs_1$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
r_bs_1$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,
                   mark.clusters=T,min.group.size=10,shuffle.colors=F,
                   mark.cluster.cex=1,alpha=0.3,main='cell clusters (tSNE) bs_1')

r_bs_2$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
set.seed(1)
r_bs_2$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
r_bs_2$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,
                     mark.clusters=T,min.group.size=10,shuffle.colors=F,
                     mark.cluster.cex=1,alpha=0.3,main='cell clusters (tSNE) bs_2')
## calculate velocity

## restrict to cells that passed p2 filter
emat_bs_filtered_1 <- dh_wt_ctrl_ms_bs_1$spliced[,rownames(r_bs_1$counts)]; 
nmat_bs_filtered_1 <- dh_wt_ctrl_ms_bs_1$unspliced[,rownames(r_bs_1$counts)]; 

emat_bs_filtered_2 <- dh_wt_ctrl_ms_bs_2$spliced[,rownames(r_bs_2$counts)]
nmat_bs_filtered_2 <- dh_wt_ctrl_ms_bs_2$unspliced[,rownames(r_bs_2$counts)]


#identity_bs$cluster<- as.character(identity_bs$cluster)
cell_assignment_bs1 <- assign_color(identity_bs_1)
cluster.label_fullseq_bs1 <- cell_assignment_bs1[1] %>% unlist()
cell.colors_fullseq_bs1 <- cell_assignment_bs1[2] %>% unlist()

cell_assignment_bs2 <- assign_color(identity_bs_2)
cluster.label_fullseq_bs2 <- cell_assignment_bs2[1] %>% unlist()
cell.colors_fullseq_bs2 <- cell_assignment_bs2[2] %>% unlist()

## calculate cell distance
cell.dist_bs_1 <- as.dist(1-armaCor(t(r_bs_1$reductions$PCA)))
cell.dist_bs_2 <- as.dist(1-armaCor(t(r_bs_2$reductions$PCA)))


## filter genes by minimum average expression magnitude, output total number of resulting 
## valid genes
emat_bs_1 <- filter.genes.by.cluster.expression(emat_bs_filtered_1,cluster.label_fullseq_bs1,min.max.cluster.average = 0.2)
nmat_bs_1 <- filter.genes.by.cluster.expression(nmat_bs_filtered_1,cluster.label_fullseq_bs1,min.max.cluster.average = 0.05)
length(intersect(rownames(emat_bs_1),rownames(nmat_bs_1))) ## 2151

emat_bs_2 <- filter.genes.by.cluster.expression(emat_bs_filtered_2,cluster.label_fullseq_bs2,min.max.cluster.average = 0.2)
nmat_bs_2 <- filter.genes.by.cluster.expression(nmat_bs_filtered_2,cluster.label_fullseq_bs2,min.max.cluster.average = 0.05)
length(intersect(rownames(emat_bs_2),rownames(nmat_bs_2))) ## 3279

### estimate RNA velocity  (using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
##bs
rvel.cd_bs_1 <- gene.relative.velocity.estimates(emat_bs_1,nmat_bs_1,kCells=25,deltaT=1,cell.dist=cell.dist_bs_1,fit.quantile=fit.quantile)
rvel.cd_bs_2 <- gene.relative.velocity.estimates(emat_bs_2,nmat_bs_2,kCells=25,deltaT=1,cell.dist=cell.dist_bs_2,fit.quantile=fit.quantile)

## set emb
emb_bs1_pagoda_tsne <- r_bs_1$embeddings$PCA$tSNE
emb_bs2_pagoda_tsne <- r_bs_2$embeddings$PCA$tSNE

show.velocity.on.embedding.cor(emb_bs1_pagoda_tsne,rvel.cd_bs_1,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq_bs1,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
show.velocity.on.embedding.cor(emb_bs2_pagoda_tsne,rvel.cd_bs_2,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq_bs2,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)



