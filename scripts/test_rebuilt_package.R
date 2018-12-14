## this is to test gene.relative.velocity.estimates_mod function write locally 

rvel.cd_mod <- gene.relative.velocity.estimates_mod(emat,nmat,deltaT=1,kCells=1,cell.dist=cell.dist,fit.quantile=fit.quantile)
show.velocity.on.embedding.cor(emb,rvel.cd_mod,n=200,scale='sqrt',cell.colors=ac(cell.colors_fullseq,alpha=0.5),cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1)


## this is to test the re-built package velocyto.YG.R

library(velocyto.YG.R)
rvel.cd_rerbuilt <- velocyto.YG.R::gene.relative.velocity.estimates(K=0.5, N=3,emat,nmat,deltaT=1,kCells=1,cell.dist=cell.dist,fit.quantile=fit.quantile)
