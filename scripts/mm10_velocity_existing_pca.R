## Project velocity on existing PCA

## this script is a general script which can project RNA velocity on pangoda generated PCA/t-SNE or existing PCA

## input: loom file, cell group file for filtration and pca file

## outpot: velocity on pca plot 

################################### main script ##########################################

## the loom file is prepared on HPC through velocyto
library(velocyto.YG.R)
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
#input_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/wound_data/wounds_data_20181106.loom"
input_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/mm10_10x_1010.loom"

mm10_loom <- read.loom.matrices(input_path)

## read in the cell list and convert cell barcode to the correct format 
cell_list_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/cell_groups/identity_CTL_BS_Krt14_Krt1_Prolif_k6.txt"
#cell_list_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/wound_data/identity_CTL_SW_Krt14_Krt1_Prolif_Daniel.txt"
cell_list <- read.delim(cell_list_path, header = FALSE, skip = 1)
#cell_list$V1 <- gsub("CTL_SW_FullSeq_", "wounds_data_20181106:", gsub(".2", "x", cell_list$V1))
cell_list$V1 <- gsub("CTL_BS_FullSeq_", "mm10_outs:", gsub(".1", "x", cell_list$V1))

## read in the pca file
pca_path <- "/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/data/wound_data/ydata_scEpath_SW.txt"
pca_file <- read.delim(pca_path, header = TRUE)
pca_file$Row <- gsub(".2", "x", pca_file$Row)
pca_file$Row <- gsub("CTL_SW_FullSeq_", "wounds_data_20181106:", pca_file$Row)

## challenge - turn the pca_file which is df into numbers with 2 attrs as r$reductions$PCA
## what i need to do: change the format of the pca to the format of r$reductions$PCA and replce the r$reductions$PCA
## everything else in the r remains the same 

## however, in the original process, they use the pca output as in the input for knnGraph and knnCluster and embedding;
## and plotted the velocity on the embedding plot eventually .. so the shape may still be altered after embedding 


## it won't work, I should check this function: pca.velocity.plot to see how they calculate and project
## velocity on pca???
## https://github.com/velocyto-team/velocyto.R/blob/master/R/momentum_routines.R


## the one below, try the show.velocity.on.embedding.cor
## try take out prolifterating cells
#cell_list <- cell_list %>% filter(V2 != "IFE-P")
emat <- filter_cell_gene(mm10_loom, cell_list, 10, 1e3)
r <- build_pagoda(emat)


## plot PCA latgeVis
M <- 30; 
set.seed(2)
r$getEmbedding(type='PCA',M=M,perplexity=30,gamma=1/M,alpha=1) ## optional for PCA

r$plotEmbedding(type='PCA',show.legend=F,
                mark.clusters=T, min.group.size=10, shuffle.colors=F, 
                mark.cluster.cex=1, alpha=0.3, main='PCA cluster')

## plot tSNE
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
set.seed(1)
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,
                mark.clusters=T,min.group.size=10,shuffle.colors=F,
                mark.cluster.cex=1,alpha=0.3,main='cell clusters (tSNE)')

## calculate velocity

## prepare matrices and clustering data
emat <- mm10_loom$spliced
nmat <- mm10_loom$unspliced
## restrict to cells that passed p2 filter
emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]

# take cluster labels from cell list

cell_assignment <- assign_color(cell_list)


cluster.label_fullseq <- cell_assignment[1] %>% unlist()
cell.colors_fullseq <- cell_assignment[2] %>% unlist()


#cluster.label <- r$clusters$PCA$multilevel # take the cluster factor that was calculated by p2

###generate t-SNE embedding
emb <- r$embeddings$PCA$tSNE

### generate PCA largeVis embedding
emb_pca_l <- r$embeddings$PCA$largeVis
emb_pca_unwounded1<-emb_pca_l
emb_pca_wounded1<-emb_pca_l
#emb_pca_l<-emb_pca_unwounded1
# transform the pca file 
emb_pca <- r$reductions$PCA

# try input only the first and second pc
emb_pca12 <- emb_pca[, c(1,2)]

## construct pca_loc with pre-generated pca
pca_file <- pca_file[pca_file$Row %in% cell_list$V1,]
names(pca_file) <- gsub("ydata", "PC", names(pca_file))
pca_loc <- as.matrix(apply(pca_file[,-1],2,as.numeric))
row.names(pca_loc) <- pca_file[,1]


## In addition to clustering and the t-SNE embedding, from the p2 processing we will also take a 
## cell-cell distance, which will be better than the default whole-transcriptome correlation 
##distance that velocyto.R would normally use.
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))


##
cell.dist_fullseq155 <- as.dist(1-armaCor(t(pca_loc)))
cell.dist_fullseq50 <- as.dist(1-armaCor(t(pca_loc[,c(1:50)])))
## filter genes by minimum average expression magnitude, output total number of resulting 
## valid genes
emat <- filter.genes.by.cluster.expression(emat,cluster.label_fullseq,min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label_fullseq,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(nmat))) ## 2165 for wounded  ## 1829 for unwounded data



### estimate RNA velocity  (using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
fit.quantile2 <- 0.2

## test re-built package

rvel.cd_rebuilt <- velocyto.YG.R::gene.relative.velocity.estimates(K=0.5, N=3, emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)
rvel.cd <- velocyto.YG.R::gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)
rvel.cd <- velocyto.R::gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)


rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)
rvel.cd2 <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist_fullseq155,fit.quantile=fit.quantile)
rvel.cd3 <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist_fullseq50,fit.quantile=fit.quantile)
rvel.cd4 <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile2)

## show rna velocity 

## tsne
par(mfrow=c(1,2)) 

show.velocity.on.embedding.cor(emb,rvel.cd,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.6),cex=0.8,arrow.scale=4,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1.1,do.par=F,cell.border.alpha = 0.1, main='default')
show.velocity.on.embedding.cor(emb,rvel.cd_rebuilt,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.6),cex=0.8,arrow.scale=4,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1.1,do.par=F,cell.border.alpha = 0.1, main='hill function')
#dev.off()

## pagoda pca+largeVis
par(mfrow=c(1,2)) 
show.velocity.on.embedding.cor(emb_pca_l,rvel.cd,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.6),cex=0.8,arrow.scale=1,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1.1,do.par=F,cell.border.alpha = 0.1, main='default')
show.velocity.on.embedding.cor(emb_pca_l,rvel.cd_rebuilt,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.6),cex=0.8,arrow.scale=1,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1.1,do.par=F,cell.border.alpha = 0.1, main='hill function')
#dev.off()
## show rna velocity on the pagoda pca

## NOTE: if the cell group changed, the color has to be re-labled 
show.velocity.on.embedding.cor(emb_pca,rvel.cd,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.5),cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)
show.velocity.on.embedding.cor(emb_pca12,rvel.cd,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.5),cex=0.8,arrow.scale=15,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)
## no proliferating 
show.velocity.on.embedding.cor(emb,rvel.cd_no_p,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.5),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

## increase quantile
show.velocity.on.embedding.cor(emb_pca,rvel.cd4,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.5),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)


## show the velocity on existing pca
## Visualize RNA velocities on an existing embedding using correlation-based transition probability matrix within the kNN graph
show.velocity.on.embedding.cor(pca_loc,rvel.cd,n=200,scale='linear',cell.colors=ac(cell.colors_fullseq,alpha=0.5),
                               cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,
                               arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)



show.velocity.on.embedding.cor(pca_loc,rvel.cd2,n=200,scale='linear',cell.colors=ac(cell.colors_fullseq,alpha=0.5),
                               cex=0.8,arrow.scale=30,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,
                               arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

show.velocity.on.embedding.cor(pca_loc,rvel.cd3,n=200,scale='linear',cell.colors=ac(cell.colors_fullseq,alpha=0.5),
                               cex=0.8,arrow.scale=30,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,
                               arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

## try using eucleadian-distance based transition probability matrix
show.velocity.on.embedding.eu(pca_loc,rvel.cd,n=30,cell.colors=ac(cell.colors_fullseq,alpha=0.5), scale='log', 
                              control.for.neighborhood.density=TRUE, ntop.trajectories=1, do.par=T, 
                              show.cell.arrows=NULL, show.cell.trajectories=NULL, show.trajectories=FALSE, 
                              show.all.trajectories=FALSE, show.cell.diffusion.posterior=NULL, 
                              diffusion.steps=10, cell.dist=NULL, trajectory.spline.shape=1, 
                              cell.color.alpha=0.5, ,arrow.scale=30, cex=0.8, show.grid.flow=TRUE,
                              min.grid.cell.mass=0.5,grid.n=40,
                              arrow.lwd=1,cell.border.alpha = 0.1) 
show.velocity.on.embedding.eu(pca_loc,rvel.cd,n=200,scale='linear',cell.colors=ac(cell.colors_fullseq,alpha=0.5),
                              cex=0.8,arrow.scale=30,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,
                              arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)

## function to filter the spliced genes based on cell list; filter the gene if it shows up in fewer than ncell;
## filter the cell if it contains fewer than ngenes
## input: loom file, cell list, ncell, ngene
## output: filtered cell list

## FUN filters the spliced loom file with cell_list 
filter_cell_gene <- function(loom_file, cell_list_file, ncell, ngene){
  
  if (!is.null(cell_list_file)){
    emat <- loom_file$spliced %>% as.matrix()%>%as.data.frame()%>%select(one_of(cell_list_file$V1))
  } else {
    emat <- loom_file$spliced
  }
  
  ## convert back to dgCMatrix
  
  
  ## chekc cell size, shown at least 1e3
  
  ### filter to omit genes show up in fewer than 10 cells
  emat <- emat[rowSums(emat)>=ncell,]
  ### filter to omit cells with small number of RNAs
  emat <- emat[,colSums(emat)>ngene] 
  
  return (emat);
  
}

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

## FUN calculate velocity, plot PCA with the velocity 
## this is the function to calculate velocity with the spliced and unspliced genes
## show the velocity on the velocyto generated pca
## input: loom
## output: velocity on pca plot; a pca file generated by velocity 

velocity_self_pca <- function (loom_file, pagoda_object) {
  
  emat <- loom_file$spliced
  nmat <- loom_file$unspliced
  ## restrict to cells that passed p2 filter
  emat <- emat[,rownames(pagoda_object$counts)]; nmat <- nmat[,rownames(pagoda_object$counts)]
  
  ##distance that velocyto.R would normally use.
  cell.dist <- as.dist(1-armaCor(t(pagoda_object$reductions$PCA)))
  
  ## take t
  cluster.label <- pagoda_object$clusters$PCA$community
  
  ## filter genes by minimum average expression magnitude, output total number of resulting 
  ## valid genes
  emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.2)
  nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
  fit.quantile <- 0.02
  
  ### calculate velocity
  ### when running the code, get error: Error: $ operator is invalid for atomic vectors
  ### In addition: Warning message:
  ### In parallel::mclapply(sn(rownames(conv.emat.norm)), function(gn) { :
  ### all scheduled cores encountered errors in user code
  rvel.cd_pca <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile_pca)
  
  ### visualize velocity on PCAs 
  velocity_pca <- pca.velocity.plot(rvel.cd_pca,return.details=T, nPcs=5,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1), max.grid.arrow.length = 0.0001, arrow.scale = 0.1)
  return (velocity_pca)
}


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


## assign color 
cell.colors_fullseq <- pagoda2:::fac2col(group)

emat <- filter_cell_gene(mm10_loom, cell_list, 10, 1e3)

r <- build_pagoda(emat)

## okay, the returned list is empty 
velocyto_pca_file <- velocity_self_pca(mm10_loom, r)




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
