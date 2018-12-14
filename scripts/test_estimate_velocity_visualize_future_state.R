## study each step of estimate velocity script and visualize the data

################################ get the unwounded 1 data ready #############

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

############################ used functions #########################
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

# transform the pca file 
emb_pca <- r$reductions$PCA


## In addition to clustering and the t-SNE embedding, from the p2 processing we will also take a 
## cell-cell distance, which will be better than the default whole-transcriptome correlation 
##distance that velocyto.R would normally use.
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

## filter genes by minimum average expression magnitude, output total number of resulting 
## valid genes
emat <- filter.genes.by.cluster.expression(emat,cluster.label_fullseq,min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label_fullseq,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(nmat))) ## 2165



### estimate RNA velocity  (using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
fit.quantile2 <- 0.2

## test re-built package

rvel.cd_rebuilt <- velocyto.YG.R::gene.relative.velocity.estimates(K=0.5, N=3, emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)
rvel.cd <- velocyto.YG.R::gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=25,cell.dist=cell.dist,fit.quantile=fit.quantile)

######################################









#######################################  study projected velocity ###########################
# estimate projected delta given x'=(y-o) - gamma*x solution
# em - normalized expression matrix
# nm - normalized nascent matrix
# gamma - inferred degradation coefficients
# o - inferred offset (assumed to be zero by default)
# delta - time to project forward
t.get.projected.delta <- function(em,nm,gamma,offset=rep(0,length(gamma)),delta=0.5) {
  # adjust rownames
  gn <- intersect(names(gamma),rownames(em)); ## select intersect between gamma and spliced
  if(is.null(names(offset))) { names(offset) <- names(gamma); }
  em <- em[gn,]; nm <- nm[gn,]; gamma <- gamma[gn]; offset <- offset[gn]; ## filter by the pass-filter genes
  # time effect constant
  egt <- exp(-gamma*delta);
  ## y is 1116(genes) x 3021 cells for unwouded 1
  y <- nm-offset; y[y<0] <- 0; # zero out entries with a negative n levels after offset adjustment
  em*egt + (1-egt)*y/gamma  - em ## each gene, each cell have a projected.delta
}

### exam the output 
gn <- intersect(names(rvel.cd$gamma),rownames(rvel.cd$conv.emat.norm))
em <- rvel.cd$conv.emat.norm[gn,]; nm <- rvel.cd$conv.nmat.norm[gn,]; gamma <- rvel.cd$gamma[gn];
delta=1
egt <- exp(-gamma*delta);
y=nm
## below is the output of t.get.projected.delta function
## each gene, each cell have one projected.delta
projected.delta_unwounded1<-em*egt + (1-egt)*y/gamma  - em  

## to add the hill function on nm(conv.nmat.norm), I analyze the dist unspliced expression used in t.get.projected.delta
## random sample 20 cell with 100 genes, plot distribution 
gene_r <- sample(1:1116,100,replace=F)
cell_r <- sample(1:3021,20,replace=F)
nm.df<- as.data.frame(as.matrix(nm[c(gene_r), c(cell_r)]))
nm.df_long <- gather(nm.df, key = cells, value = expression)
ggplot_th <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12), 
                   legend.title = element_text(size=12),
                   title = element_text(size=12))
png("/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/hill_func1124/explore/normalized_unpliced_expression.png", width=600, height = 300)
ggplot(nm.df_long, aes(x=expression)) +
  geom_histogram(binwidth=.1, colour="black", fill="white")+xlim(0, 5) +  
  theme_bw()+
  theme_set(ggplot_th)
dev.off()
summary(nm.df_long$expression)

## plot the effect after hill function
## select k=1 and n=4
nm.df_long$k1_n4 <-nm.df_long$expression^4/(nm.df_long$expression^4+1^4)
nm.df_long$k05_n4 <-nm.df_long$expression^4/(nm.df_long$expression^4+0.5^4)
nm.df_long$k05_n3 <-nm.df_long$expression^3/(nm.df_long$expression^3+0.5^3)
nm.df_long$k1_n3 <-nm.df_long$expression^3/(nm.df_long$expression^3+1^3)
nm.df_long_long <- gather(nm.df_long, key =hill, value = expression_hill, -cells, -expression)


png("/Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/plots/hill_func1124/explore/hill_vs_noHill_deltaS.png", width=500, height = 300)
ggplot(aes(x=expression, y=expression_hill, group=hill, color=hill), data = nm.df_long_long)+
  geom_line()+xlim(0,3)+ggtitle("Apply hill function on unspliced expression")+
  theme_bw()+
  theme_set(ggplot_th)
dev.off()

######################## study function to generate rvel.cd$projected ###################
# calculates the difference in the number of counts based on the library size, renormalizes
# note: also introduces centripetal velocity
emn <- t.get.projected.cell2(emat.norm,emat.sz,as.matrix(deltaE),mult = mult,delta=deltaT2);
resl <- c(resl,list(projected=emn,current=emat.norm,deltaE=deltaE,deltaT=deltaT,ko=full.ko,mult=mult,kCells=kCells));
deltae <- as.matrix(em*egt + (1-egt)*y/gamma  - em)
## size: sum of gene expression of each cell
mult=1e3
emat.size <- Matrix::colSums(emat);nmat.size <- Matrix::colSums(nmat);
emat.cs <- emat.size[colnames(emat)]/mult;
nmat.cs <- nmat.size[colnames(nmat)]/mult;
emat.sz <- emat.cs; cellSize<-emat.sz
## reduced cell normalization on genes for which moentun was used
emat.norm <- emat[rownames(emat) %in% rownames(deltae),]
emat.norm <- t(t(emat.norm)/(emat.sz));

## t.get.projected.cell2 function
## construct 0 matrix with the same colname and rowname of emat.norm
rz <- matrix(0,nrow=nrow(emat.norm),ncol=ncol(emat.norm)); colnames(rz) <- colnames(emat.norm); rownames(rz) <- rownames(emat.norm)
gn <- intersect(rownames(deltae),rownames(rz))
## rz is identical to deltae for my test case
rz[match(gn,rownames(rz)),colnames(deltae)] <- deltae[gn,]; 
# translate fpm delta into the number of molecules based on the current cell size
rz <- t(t(rz)*emat.sz)
emm <- t(t(em)*emat.sz)
## this is the equation of S = S0 + vt
emn <- emm + rz*delta;
emn[emn<0] <- 0;
newCellSize <- (emat.sz+Matrix::colSums(emn-emm)/mult)
## at the end, still return the 1116(genes) x 3021(cells) array  -------------
emn <- t(t(emn)/newCellSize)


t.get.projected.cell2 <- function(em,cellSize,deltae,mult=1e3,delta=1) {
  rz <- matrix(0,nrow=nrow(em),ncol=ncol(em)); colnames(rz) <- colnames(em); rownames(rz) <- rownames(em)
  gn <- intersect(rownames(deltae),rownames(rz))
  rz[match(gn,rownames(rz)),colnames(deltae)] <- deltae[gn,]; 
  # translate fpm delta into the number of molecules based on the current cell size
  rz <- t(t(rz)*cellSize)
  emm <- t(t(em)*cellSize)
  emn <- emm + rz*delta; 
  emn[emn<0] <- 0;
  newCellSize <- (cellSize+Matrix::colSums(emn-emm)/mult)
  emn <- t(t(emn)/newCellSize)
  
  #emn <- t(t(emn)/Matrix::colSums(emn)*Matrix::colSums(em))
  emn
}

## this does not include unspliced reads, pause for now 

