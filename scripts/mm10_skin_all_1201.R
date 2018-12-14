## analyze 2 unwound(BS) and 3 wound (SW) data
## 12/01/2018
## First combine 2 unwounded data into 1; 3 wounded data into 1
## Then, projected separated on four types of low dimensional space (PCA, DC, monocle and umap)
## color code based on provided list 

## the loom file is prepared on HPC through velocyto

library(velocyto.R)
library(pagoda2)
library(igraph)
library(tidyverse)
library(Matrix)
#library(loomR)


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

## FUN generate group label and color label from identity and color list
## this function will take color list as input and add cell name as the row name
## input: identity list, color list
## output: group assignment, color assignment for each cell

assign_group_color <- function(cell_list, color_list){
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
  
  cell_color_list <- cbind(cell_list, color_list)
  for (j in 1:nrow(cell_color_list)){
    print(j)
    cell_i <- cell_color_list[j,1];
    color_i <- cell_color_list[j,3];
    names(color_i)<-cell_i
    
    if (j == 1){
      color <- color_i;
    } else {
      color <- c(color, color_i);
    }
  }
  
  return(list(group, color))
  
  
}

## read in cell group and color file
file_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/mm10_mouse_skin_loom_all_1201"
cellColors_bs <- read.delim(paste(file_path, "niedaigroupmeeting/colorCells_bs_epi_cca_basal_spinious_prolif.txt", sep="/"), header = FALSE)
identity_bs <- read.delim(paste(file_path, "niedaigroupmeeting/identity_bs_epi_cca_basal_spinious_prolif.txt", sep="/"), header = TRUE) 
identity_bs$Row <- gsub("..$", "x", identity_bs$Row)
identity_bs$Row <- gsub("bs_1_", "mm10_outs:", gsub("bs_2_", "dh_wt_ctrl_ms_bs_2:", identity_bs$Row))
## sw
cellColors_sw <- read.delim(paste(file_path, "niedaigroupmeeting/colorCells_sw_epi_cca_basal_spinious_prolif.txt", sep="/"), header = FALSE)
identity_sw <- read.delim(paste(file_path, "niedaigroupmeeting/identity_sw_epi_cca_basal_spinious_prolif.txt", sep="/"), header = TRUE) 
identity_sw$Row <- gsub("..$", "x", identity_sw$Row)
identity_sw$Row <- gsub("sw_1_", "wounds_data_20181106:", gsub("sw_2_", "dh_wt_ctrl_ms_sw_2:",gsub("sw_3_", "dh_wt_ctrl_ms_sw_3:", identity_sw$Row)))


## read in loom files and concatenate
## bs
dh_wt_ctrl_ms_bs_2_path <- paste(file_path, "dh_wt_ctrl_ms_bs_2.loom", sep="/")
dh_wt_ctrl_ms_bs_1_path <- paste(file_path, "dh_wt_ctrl_ms_bs_1.loom", sep="/")
dh_wt_ctrl_ms_bs_2 <- read.loom.matrices(dh_wt_ctrl_ms_bs_2_path)
dh_wt_ctrl_ms_bs_1 <- read.loom.matrices(dh_wt_ctrl_ms_bs_1_path)
##sw
dh_wt_ctrl_ms_sw_1_path <- paste(file_path, "dh_wt_ctrl_ms_sw_1.loom", sep="/")
dh_wt_ctrl_ms_sw_2_path <- paste(file_path, "dh_wt_ctrl_ms_sw_2.loom", sep="/")
dh_wt_ctrl_ms_sw_3_path <- paste(file_path, "dh_wt_ctrl_ms_sw_3.loom", sep="/")
dh_wt_ctrl_ms_sw_1 <- read.loom.matrices(dh_wt_ctrl_ms_sw_1_path)
dh_wt_ctrl_ms_sw_2 <- read.loom.matrices(dh_wt_ctrl_ms_sw_2_path)
dh_wt_ctrl_ms_sw_3 <- read.loom.matrices(dh_wt_ctrl_ms_sw_3_path)

##combining loom file
# emat_bs_1 <- dh_wt_ctrl_ms_bs_1$spliced %>%as.matrix()%>%as.data.frame()
# emat_bs_2 <- dh_wt_ctrl_ms_bs_2$spliced %>%as.matrix()%>%as.data.frame()
emat_bs_1 <- dh_wt_ctrl_ms_bs_1$spliced 
emat_bs_2 <- dh_wt_ctrl_ms_bs_2$spliced 
emat_combined_bs<- cbind(emat_bs_1, emat_bs_2)
# names(emat_bs_1) <- gsub("mm10_outs:", "bs_1_", colnames(emat_bs_1))
# names(emat_bs_2) <- gsub("dh_wt_ctrl_ms_bs_2:", "bs_2_", colnames(emat_bs_2))
# emat_bs_12 <- cbind(emat_bs_1[ order(row.names(emat_bs_1)), ], emat_bs_2[ order(row.names(emat_bs_2)), ])
# emat_bs_12_matrix <- as.matrix(emat_bs_12, sparse = TRUE)

nmat_bs_1 <- dh_wt_ctrl_ms_bs_1$unspliced
nmat_bs_2 <- dh_wt_ctrl_ms_bs_2$unspliced 
nmat_combined_bs<- cbind(nmat_bs_1, nmat_bs_2)
# nmat_bs_1 <- dh_wt_ctrl_ms_bs_1$unspliced %>%as.matrix()%>%as.data.frame()
# nmat_bs_2 <- dh_wt_ctrl_ms_bs_2$unspliced %>%as.matrix()%>%as.data.frame()
# names(nmat_bs_1) <- gsub("mm10_outs:", "bs_1_", colnames(nmat_bs_1))
# names(nmat_bs_2) <- gsub("dh_wt_ctrl_ms_bs_2:", "bs_2_", colnames(nmat_bs_2))
# nmat_bs_12 <- cbind(nmat_bs_1[ order(row.names(nmat_bs_1)), ], emat_bs_2[ order(row.names(nmat_bs_2)), ])
# nmat_bs_12_matrix <- as.matrix(nmat_bs_12, sparse = TRUE)

## sw
emat_sw_1 <- dh_wt_ctrl_ms_sw_1$spliced 
emat_sw_2 <- dh_wt_ctrl_ms_sw_2$spliced 
emat_sw_3 <- dh_wt_ctrl_ms_sw_3$spliced 
emat_combined_sw<- cbind(emat_sw_1, emat_sw_2, emat_sw_3)
nmat_sw_1 <- dh_wt_ctrl_ms_sw_1$unspliced 
nmat_sw_2 <- dh_wt_ctrl_ms_sw_2$unspliced 
nmat_sw_3 <- dh_wt_ctrl_ms_sw_3$unspliced 
nmat_combined_sw<- cbind(nmat_sw_1, nmat_sw_2, nmat_sw_3)


## filter emat based on cell_list
emat_bs_12_filtered <- emat_combined_bs%>%as.matrix()%>%as.data.frame()%>%select(one_of(identity_bs$Row))
r_bs <- build_pagoda(emat_bs_12_filtered)
## sw
emat_sw_12_filtered <- emat_combined_sw%>%as.matrix()%>%as.data.frame()%>%select(one_of(identity_sw$Row))
r_sw <- build_pagoda(emat_sw_12_filtered)


## plot tSNE
r_bs$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
set.seed(1)
r_bs$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
r_bs$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,
                mark.clusters=T,min.group.size=10,shuffle.colors=F,
                mark.cluster.cex=1,alpha=0.3,main='cell clusters (tSNE)')
## sw
r_sw$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
set.seed(1)
r_sw$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
r_sw$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,
                   mark.clusters=T,min.group.size=10,shuffle.colors=F,
                   mark.cluster.cex=1,alpha=0.3,main='cell clusters (tSNE)')

## calculate velocity

## restrict to cells that passed p2 filter
emat_bs <- emat_combined_bs[,rownames(r_bs$counts)]; 
nmat_bs <- nmat_combined_bs[,rownames(r_bs$counts)]
## sw
emat_sw <- emat_combined_sw[,rownames(r_sw$counts)]; 
nmat_sw <- nmat_combined_sw[,rownames(r_sw$counts)]


## assign cell group and color
#identity_bs$cluster<- as.character(identity_bs$cluster)
cellColors_bs$V1<- as.character(cellColors_bs$V1)
cell_assignment <- assign_group_color(identity_bs, cellColors_bs)
assignedClusters_bs <- cell_assignment[1] %>% unlist()
assignedColors_bs <- cell_assignment[2] %>% unlist()

##sw
cellColors_sw$V1<- as.character(cellColors_sw$V1)
cell_assignment_sw <- assign_group_color(identity_sw, cellColors_sw)
assignedClusters_sw <- cell_assignment_sw[1] %>% unlist()
assignedColors_sw <- cell_assignment_sw[2] %>% unlist()

## calculate cell distance
cell.dist_bs <- as.dist(1-armaCor(t(r_bs$reductions$PCA)))
cell.dist_sw <- as.dist(1-armaCor(t(r_sw$reductions$PCA)))

## filter genes by minimum average expression magnitude, output total number of resulting 
## valid genes
emat_bs <- filter.genes.by.cluster.expression(emat_bs,assignedClusters_bs,min.max.cluster.average = 0.2)
nmat_bs <- filter.genes.by.cluster.expression(nmat_bs,assignedClusters_bs,min.max.cluster.average = 0.05)
length(intersect(rownames(emat_bs),rownames(nmat_bs))) ## 2497
##sw
emat_sw <- filter.genes.by.cluster.expression(emat_sw,assignedClusters_sw,min.max.cluster.average = 0.2)
nmat_sw <- filter.genes.by.cluster.expression(nmat_sw,assignedClusters_sw,min.max.cluster.average = 0.05)
length(intersect(rownames(emat_sw),rownames(nmat_sw))) ## 2506

### estimate RNA velocity  (using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
##bs
rvel.cd_bs <- gene.relative.velocity.estimates(emat_bs,nmat_bs,kCells=25,deltaT=1,cell.dist=cell.dist_bs,fit.quantile=fit.quantile)
rvel.cd_bs_k50 <- gene.relative.velocity.estimates(emat_bs,nmat_bs,kCells=50,deltaT=1,cell.dist=cell.dist_bs,fit.quantile=fit.quantile)
rvel.cd_bs_k200 <- gene.relative.velocity.estimates(emat_bs,nmat_bs,kCells=20,deltaT=1,cell.dist=cell.dist_bs,fit.quantile=fit.quantile)

##sw
rvel.cd_sw <- gene.relative.velocity.estimates(emat_sw,nmat_sw,kCells=25,deltaT=1,cell.dist=cell.dist_sw,fit.quantile=fit.quantile)

## set emb
emb_bs_pagoda_tsne <- r_bs$embeddings$PCA$tSNE
emb_sw_pagoda_tsne <- r_sw$embeddings$PCA$tSNE

pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_bs_epi_cca_basal_spinious_prolif_tsne.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(emb_bs_pagoda_tsne,rvel.cd_bs,n=200,scale='sqrt',cell.colors=ac(assignedColors_bs,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)
dev.off()
##sw
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_sw_epi_cca_basal_spinious_prolif_tsne.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(emb_sw_pagoda_tsne,rvel.cd_sw,n=500,scale='sqrt',cell.colors=ac(assignedColors_sw,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=0.5,do.par=T,cell.border.alpha = 0.001, xlab = "Component 1", ylab = "Component 2")
dev.off()

## load position file

## FUN prepare low dimension file for plotting
## Input is the file low dimension file name; the function is to organize the Row column;
## and convert the file to correct format

generate_loc_file <- function(low_dim_file){
  low_dim_file$Row <- gsub("..$", "x", low_dim_file$Row)
  low_dim_file$Row <- gsub("bs_1_", "mm10_outs:", gsub("bs_2_", "dh_wt_ctrl_ms_bs_2:", low_dim_file$Row))
  loc <- as.matrix(apply(low_dim_file[,-1],2,as.numeric))
  row.names(loc) <- low_dim_file[,1]
  return(loc)
}

generate_loc_file_sw <- function(low_dim_file){
  low_dim_file$Row <- gsub("..$", "x", low_dim_file$Row)
  low_dim_file$Row <- gsub("sw_1_", "wounds_data_20181106:", gsub("sw_2_", "dh_wt_ctrl_ms_sw_2:",gsub("sw_3_", "dh_wt_ctrl_ms_sw_3:", low_dim_file$Row)))
  loc <- as.matrix(apply(low_dim_file[,-1],2,as.numeric))
  row.names(loc) <- low_dim_file[,1]
  return(loc)
}

## PCA
bs_pca<-read.delim(paste(file_path, "niedaigroupmeeting/dr_pca_bs_epi_cca_basal_spinious_prolif.txt", sep="/"))
pca_loc_bs_func<-generate_loc_file(bs_pca)
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_bs_epi_cca_basal_spinious_prolif_pca.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(pca_loc_bs_func,rvel.cd_bs,n=200,scale='sqrt',cell.colors=ac(assignedColors_bs,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)
dev.off()
show.velocity.on.embedding.cor(pca_loc_bs_func,rvel.cd_bs_k50,n=500,scale='sqrt',cell.colors=ac(assignedColors_bs,alpha=0.8),cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1.5,do.par=T,cell.border.alpha = 0.001, xlab='Component 1', ylab = 'Component 2')
show.velocity.on.embedding.cor(pca_loc_bs_func,rvel.cd_bs,n=100,scale='sqrt',cell.colors=ac(assignedColors_bs,alpha=0.8),cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1.5,do.par=T,cell.border.alpha = 0.001, xlab='Component 1', ylab = 'Component 2')


## DC
bs_DC<-read.delim(paste(file_path, "niedaigroupmeeting/dr_DC_bs_epi_cca_basal_spinious_prolif.txt", sep="/"))
DC_loc_bs<-generate_loc_file(bs_DC)
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_bs_epi_cca_basal_spinious_prolif_DC.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(DC_loc_bs,rvel.cd_bs,n=200,scale='sqrt',cell.colors=ac(assignedColors_bs,alpha=0.8),cex=0.8,arrow.scale=0.1,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)
dev.off()


## monocole
bs_mono<-read.delim(paste(file_path, "niedaigroupmeeting/dr_monocle_bs_epi_cca_basal_spinious_prolif.txt", sep="/"))
bs_mono$Row <- gsub("..$", "x", bs_mono$Row)
bs_mono$Row <- gsub("bs_1_", "mm10_outs:", gsub("bs_2_", "dh_wt_ctrl_ms_bs_2:", bs_mono$Row))

## construct DC_loc with pre-generated DC
#names(bs_DC) <- gsub("ydata", "PC", names(bs_pca))
mono_loc_bs <- as.matrix(apply(bs_mono[,-1],2,as.numeric))
row.names(mono_loc_bs) <- bs_mono[,1]
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_bs_epi_cca_basal_spinious_prolif_DC.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(mono_loc_bs,rvel.cd_bs,n=200,scale='sqrt',cell.colors=ac(assignedColors_bs,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)
dev.off()

show.velocity.on.embedding.cor(mono_loc_bs,rvel.cd_bs,n=50,scale='sqrt',cell.colors=ac(assignedColors_bs,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)
show.velocity.on.embedding.cor(mono_loc_bs,rvel.cd_bs_k200,n=200,scale='sqrt',cell.colors=ac(assignedColors_bs,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)


## ump
## manually added the colname of "Row"
bs_ump<-read.delim(paste(file_path, "niedaigroupmeeting/projectedData_ump_bs_epi_cca_basal_spinious_prolif.txt", sep="/"))
bs_ump$Row <- as.character(bs_ump$Row)
bs_ump$Row <- gsub("..$", "x", bs_ump$Row)
bs_ump$Row <- gsub("bs_1_", "mm10_outs:", gsub("bs_2_", "dh_wt_ctrl_ms_bs_2:", bs_ump$Row))

## construct DC_loc with pre-generated DC
#names(bs_DC) <- gsub("ydata", "PC", names(bs_pca))
ump_loc_bs <- as.matrix(apply(bs_ump[,-1],2,as.numeric))
row.names(ump_loc_bs) <- bs_ump[,1]
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_bs_epi_cca_basal_spinious_prolif_DC.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(ump_loc_bs,rvel.cd_bs,n=200,scale='sqrt',cell.colors=ac(assignedColors_bs,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1, xlab = "Component 1", ylab = "Component 2")
dev.off()


## sw
## PCA
sw_pca<-read.delim(paste(file_path, "niedaigroupmeeting/dr_pca_sw_epi_cca_basal_spinious_prolif.txt", sep="/"))
pca_loc_sw_func<-generate_loc_file_sw(sw_pca)
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_sw_epi_cca_basal_spinious_prolif_pca.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(pca_loc_sw_func,rvel.cd_sw,n=500,scale='sqrt',cell.colors=ac(assignedColors_sw,alpha=0.8),cex=0.8,arrow.scale=8,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1.5,do.par=T,cell.border.alpha = 0.001,  xlab = "Component 1", ylab = "Component 2")
dev.off()

## DC
sw_DC<-read.delim(paste(file_path, "niedaigroupmeeting/dr_DC_sw_epi_cca_basal_spinious_prolif.txt", sep="/"))
DC_loc_sw<-generate_loc_file_sw(sw_DC)
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_sw_epi_cca_basal_spinious_prolif_DC.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(DC_loc_sw,rvel.cd_sw,n=500,scale='sqrt',cell.colors=ac(assignedColors_sw,alpha=0.8),cex=0.8,arrow.scale=0.1,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=0.5,do.par=T,cell.border.alpha = 0.001,  xlab = "Component 1", ylab = "Component 2")
dev.off()


## monocole
sw_mono<-read.delim(paste(file_path, "niedaigroupmeeting/dr_monocle_sw_epi_cca_basal_spinious_prolif.txt", sep="/"))
mono_loc_sw<-generate_loc_file_sw(sw_mono)
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_sw_epi_cca_basal_spinious_prolif_mono.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(mono_loc_sw,rvel.cd_sw,n=500,scale='sqrt',cell.colors=ac(assignedColors_sw,alpha=0.8),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=0.5,do.par=T,cell.border.alpha = 0.001,  xlab = "Component 1", ylab = "Component 2")
dev.off()

## ump
## manually added the colname of "Row"
sw_ump<-read.delim(paste(file_path, "niedaigroupmeeting/projectedData_ump_sw_epi_cca_basal_spinious_prolif.txt", sep="/"))
ump_loc_sw<-generate_loc_file_sw(sw_ump)
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/velocity_sw_epi_cca_basal_spinious_prolif_ump.pdf', width = 5, height = 4.5)
show.velocity.on.embedding.cor(ump_loc_sw,rvel.cd_sw,n=500,scale='sqrt',cell.colors=ac(assignedColors_sw,alpha=0.8),cex=0.8,arrow.scale=2,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=0.5,do.par=T,cell.border.alpha = 0.001,  xlab = "Component 1", ylab = "Component 2")
dev.off()





### plot the name and color
## identity_bs, cellColors_bs
cluster_color <- cbind(identity_bs, cellColors_bs) 
cluster_color$cluster<- as.character(cluster_color$cluster)
cluster_color_summary <- cluster_color%>%group_by(cluster, V1)%>%summarise(count=n())%>%as.data.frame()
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/cell_count_summary.pdf', width = 7, height = 4.5)
cluster_color_summary %>% 
  ggplot(aes(x = cluster, y = count, fill=cluster)) + geom_bar(stat="identity") +
  scale_fill_manual(values =  as.character(unlist(cluster_color_summary$V1))) + geom_text(aes(label=count))+
  xlab("cell type") + ylab("counts")+
  theme(axis.text.x  = element_text(angle=30, vjust=0.5, size=12),
        axis.title.x  = element_text(size=12),
        legend.position = "none")
dev.off()
##sw
cluster_color_sw <- cbind(identity_sw, cellColors_sw) 
cluster_color_sw$cluster<- as.character(cluster_color_sw$cluster)
cluster_color_summary <- cluster_color_sw%>%group_by(cluster, V1)%>%summarise(count=n())%>%as.data.frame()
pdf('/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/mm10_skin_all_1203/cell_count_summary_sw.pdf', width = 7, height = 4.5)
cluster_color_summary %>% 
  ggplot(aes(x = cluster, y = count, fill=cluster)) + geom_bar(stat="identity") +
  scale_fill_manual(values =  as.character(unlist(cluster_color_summary$V1))) + geom_text(aes(label=count))+
  xlab("cell type") + ylab("counts")+
  theme(axis.text.x  = element_text(angle=30, vjust=0.5, size=12),
        axis.title.x  = element_text(size=12),
        legend.position = "none")
dev.off()

