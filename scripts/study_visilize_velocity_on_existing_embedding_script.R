##' Visualize RNA velocities on an existing embedding using correlation-based transition probability matrix within the kNN graph
##'
##' @param emb embedding onto which to project the velocities; The dimensions of coordinates should be on the order of 10x10 for the default values to make sense.
##' @param vel velocity estimates (e.g. returned by gene.relative.velocity.estimates() )
##' @param n neighborhood size (default=100 cells)
##' @param cell.colors name vector of cell colors
##' @param corr.sigma sigma parameter used to translate velocity-(expression delta) correlation into a transition probability
##' @param show.grid.flow whether to show grid velocity summary
##' @param grid.n number of grid points along each axis
##' @param grid.sd standard deviation (in embedding coordinate space) used to determine the weighting of individual cells around each grid point
##' @param min.grid.cell.mass minimal cell "mass" (weighted number of cells) around each grid point required for it to show up
##' @param min.arrow.size minimal arrow size
##' @param arrow.scale arrow scale multiplier
##' @param max.grid.arrow.length minimal arrow size
##' @param fixed.arrow.length whether to use fixed arrow width (default=FALSE)
##' @param plot.grid.points whether to mark all grid points with dots (even if they don't have valid velocities)
##' @param scale velocity scale to use (default: 'log', other values: 'sqrt','rank','linear')
##' @param nPcs number of PCs to use for velocity regularization (default NA, turns off regularization)
##' @param arrow.lwd arrow width (under fixed.arrow.length=T)
##' @param xlab x axis label
##' @param ylab y axls label
##' @param n.cores number of cores to use
##' @param do.par whether to reset plotting parameters
##' @param show.cell whether to show detailed velocity estimates for a specified cell
##' @param cell.border.alpha transparency for the cell border
##' @param cc velocity-(exprssion delta) correlation matrix (can be passed back from previous results, as $cc) to save calculation time when replotting the same velocity estimates on the same embedding with different parameters
##' @param return.details whether to return detailed output (which can be passed to p1 app for visualization)
##' @param expression.scaling whether to scale the velocity length by the projection of velocity onto the expected expression change (based on the transition probability matrix)
##' @param ... extra parameters passed to plot() function
##' @return if return.details=F, returns invisible list containing transition probability matrix ($tp) and the velocity-(expression delta) correlation matrix ($cc). If return.details=T, returns a more extended list that can be passed as veloinfo to pagoda2::p2.make.pagoda1.app() for visualization
##' @export
show.velocity.on.embedding.cor <- function(emb,vel,n=100,cell.colors=NULL, corr.sigma=0.05, show.grid.flow=FALSE, grid.n=20, grid.sd=NULL, min.grid.cell.mass=1, min.arrow.size=NULL, arrow.scale=1, max.grid.arrow.length=NULL, fixed.arrow.length=FALSE, plot.grid.points=FALSE, scale='log', nPcs=NA,  arrow.lwd=1, xlab="", ylab="", n.cores=defaultNCores(), do.par=T, show.cell=NULL, cell.border.alpha=0.3,cc=NULL, return.details=FALSE, expression.scaling=FALSE,  ...) {
  randomize <- FALSE;
  if(do.par) par(mfrow=c(1,1), mar = c(3.5,3.5,2.5,1.5), mgp = c(2,0.65,0), cex = 0.85);
  celcol <- 'white'
  if(is.null(show.cell)) { celcol <- cell.colors[rownames(emb)] }
  ## this is to just plot the t-SNE or whatever embedding it is
  plot(emb,bg=celcol,pch=21,col=ac(1,alpha=cell.border.alpha), xlab=xlab, ylab=ylab, ...);
  
  #plot(emb,bg=cell.colors[rownames(emb)],pch=21,col=ac(1,alpha=0.3), xlab=xlab, ylab=ylab);
  ## the normalized matrix of spliced gene count
  em <- as.matrix(vel$current); 
  ccells <- intersect(rownames(emb),colnames(em));
  em <- em[,ccells]; emb <- emb[ccells,]
  nd <- as.matrix(vel$deltaE[,ccells])
  cgenes <- intersect(rownames(em),rownames(nd));
  nd <- nd[cgenes,]; em <- em[cgenes,]
  #browser()
  if(randomize) {
    # randomize cell and sign for each gene
    nd <- t(apply(nd,1,function(x) (rbinom(length(x),1,0.5)*2-1)*abs(sample(x))))
  }
  #vg <- rownames(em) %in% rownames(r)
  
  ## skip, cc is null
  if(is.null(cc)) {
    # cosine projections
    cat("delta projections ... ")
    
    if(scale=='log') {
      cat("log ")
      cc <- colDeltaCorLog10(em,(log10(abs(nd)+1)*sign(nd)),nthreads=n.cores);
    } else if(scale=='sqrt') {
      cat("sqrt ")
      cc <- colDeltaCorSqrt(em,(sqrt(abs(nd))*sign(nd)),nthreads=n.cores);
    } else if(scale=='rank') {
      cat("rank ")
      cc <- colDeltaCor((apply(em,2,rank)),(apply(nd,2,rank)),nthreads=n.cores);
    } else { # linear
      cat("linear ")
      cc <- colDeltaCor(em,nd,nthreads=n.cores);
    }
    colnames(cc) <- rownames(cc) <- colnames(em)
    diag(cc) <- 0;
  }
  
  ## start
  cat("knn ... ")
  if(n>nrow(cc)) { n <- nrow(cc) }
  # TODO: add kNN based on high-dimensional correlation or Euclidean distances
  # define kNNs based on the embedding (L2 distance)
  emb.knn <- balancedKNN(t(emb),k=n,maxl=nrow(emb),dist='euclidean',n.threads=n.cores)
  diag(emb.knn) <- 1
  # caluclate transition probabilities (from col to row)
  cat("transition probs ... ")
  tp <- exp(cc/corr.sigma)*emb.knn
  #diag(tp) <- 0; #  should we allow the self-corelation for scaling?
  tp <- t(t(tp)/Matrix::colSums(tp)); # tp shows transition from a given column cell to different row cells
  tp <- as(tp,'dgCMatrix')
  cat("done\n")
  if(!is.null(show.cell)) {
    i <- match(show.cell,rownames(emb));
    if(is.na(i)) stop(paste('specified cell',i,'is not in the embedding'))
    # plot transition prob for a given cell
    points(emb,pch=19,col=ac(val2col(tp[rownames(emb),show.cell],gradient.range.quantile=1),alpha=0.5))
    points(emb[show.cell,1],emb[show.cell,2],pch=3,cex=1,col=1)
    di <- t(t(emb)-emb[i,])
    di <- di/sqrt(Matrix::rowSums(di^2))*arrow.scale; di[i,] <- 0;
    dir <- Matrix::colSums(di*tp[,i]) 
    dic <- Matrix::colSums(di*(tp[,i]>0)/sum(tp[,i]>0)); # relative to expected kNN center
    dia <- dir-dic;
    #browser()
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dic[1],emb[colnames(em)[i],2]+dic[2],length=0.05,lwd=1,col='blue'))
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dir[1],emb[colnames(em)[i],2]+dir[2],length=0.05,lwd=1,col='red'))
    suppressWarnings(arrows(emb[colnames(em)[i],1]+dic[1],emb[colnames(em)[i],2]+dic[2],emb[colnames(em)[i],1]+dir[1],emb[colnames(em)[i],2]+dir[2],length=0.05,lwd=1,lty=1,col='grey50'))
    suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+dia[1],emb[colnames(em)[i],2]+dia[2],length=0.05,lwd=1,col='black'))
  } else {
    # arrow estimates for each cell
    cat("calculating arrows ... ")
    arsd <- data.frame(t(embArrows(emb,tp,arrow.scale,n.cores)))
    rownames(arsd) <- rownames(emb);
    
    if(expression.scaling) {
      tpb <- tp>0; tpb <- t(t(tpb)/colSums(tpb));
      es <- as.matrix(em %*% tp) -as.matrix(em %*% as.matrix(tpb));
      # project velocity onto expression shift
      #pm <- as.matrix(t(vel$deltaE)/sqrt(colSums(vel$deltaE*vel$deltaE)))[colnames(es),] * (t(es)/sqrt(colSums(es*es)))
      #pl <- pmax(0,apply(pm,1,sum))
      pl <- pmin(1,pmax(0,apply(as.matrix(vel$deltaE[,colnames(es)]) * es, 2, sum)/sqrt(colSums(es*es))))
      
      
      arsd <- arsd * pl;
    }
    
    
    ars <- data.frame(cbind(emb,emb+arsd));
    colnames(ars) <- c('x0','y0','x1','y1')
    colnames(arsd) <- c('xd','yd')
    rownames(ars) <- rownames(emb);
    cat("done\n")
    
    ## we normally set this to true
    if(show.grid.flow) { # show grid summary of the arrows
      
      # set up a grid
      cat("grid estimates ... ")
      rx <- range(c(range(ars$x0),range(ars$x1)))
      ry <- range(c(range(ars$y0),range(ars$y1)))
      gx <- seq(rx[1],rx[2],length.out=grid.n)
      gy <- seq(ry[1],ry[2],length.out=grid.n)
      
      # for each grid point calculate Gaussian-weighted delta average
      if(is.null(grid.sd)) {
        grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
        cat("grid.sd=",grid.sd," ")
      }
      if(is.null(min.arrow.size)) {
        min.arrow.size <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)*1e-2;
        cat("min.arrow.size=",min.arrow.size," ")
      }
      if(is.null(max.grid.arrow.length)) {
        max.grid.arrow.length <- sqrt(sum((par('pin')/c(length(gx),length(gy)))^2))*0.25
        cat("max.grid.arrow.length=",max.grid.arrow.length," ")
      }
      
      garrows <- do.call(rbind,lapply(gx,function(x) {
        # cell distances (rows:cells, columns: grid points)
        cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
        cw <- dnorm(cd,sd=grid.sd)
        # calculate x and y delta expectations
        gw <- Matrix::colSums(cw)
        cws <- pmax(1,Matrix::colSums(cw));
        gxd <- Matrix::colSums(cw*arsd$xd)/cws
        gyd <- Matrix::colSums(cw*arsd$yd)/cws
        
        al <- sqrt(gxd^2+gyd^2);
        vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
        
        cbind(rep(x,sum(vg)),gy[vg],x+gxd[vg],gy[vg]+gyd[vg])
      }))
      colnames(garrows) <- c('x0','y0','x1','y1')
      
      # plot
      if(fixed.arrow.length) {
        suppressWarnings(arrows(garrows[,1],garrows[,2],garrows[,3],garrows[,4],length=0.05,lwd=arrow.lwd))
      } else {
        alen <- pmin(max.grid.arrow.length,sqrt( ((garrows[,3]-garrows[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + ((garrows[,4]-garrows[,2])*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2))
        # can't specify different arrow lengths in one shot :(
        #suppressWarnings(arrows(garrows[,1],garrows[,2],garrows[,3],garrows[,4],length=alen,lwd=arrow.lwd))
        suppressWarnings(lapply(1:nrow(garrows),function(i) arrows(garrows[i,1],garrows[i,2],garrows[i,3],garrows[i,4],length=alen[i],lwd=arrow.lwd)))
      }
      if(plot.grid.points) points(rep(gx,each=length(gy)),rep(gy,length(gx)),pch='.',cex=1e-1,col=ac(1,alpha=0.4))
      
      cat("done\n")
      
      if(return.details) { # for the p1 app
        # calculate expression shift
        cat("expression shifts .")
        # for individual cells
        
        scale.int <- switch(scale,'log'=2,'sqrt'=3,1)
        #es <- expectedExpressionShift(e=as.matrix(em),tp=tp,scale=scale.int,nthreads=n.cores); colnames(es) <- colnames(em); rownames(es) <- rownames(em);
        if(!expression.scaling) { #otherwise it has already been calculated
          tpb <- tp>0; tpb <- t(t(tpb)/colSums(tpb));
          #es <- expectedExpressionShift(e=as.matrix(em %*% as.matrix(tpb)),tp=tp,scale=scale.int,nthreads=n.cores); colnames(es) <- colnames(em); rownames(es) <- rownames(em);
          es <- as.matrix(em %*% tp) -as.matrix(em %*% as.matrix(tpb));
        }
        cat(".");
        # for the grid
        gs <- do.call(cbind,parallel::mclapply(gx,function(x) {
          # cell distances (rows:cells, columns: grid points)
          cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
          cw <- dnorm(cd,sd=grid.sd)
          # calculate x and y delta expectations
          gw <- Matrix::colSums(cw)
          cws <- pmax(1,Matrix::colSums(cw));
          cw <- t(t(cw)/cws)
          gxd <- Matrix::colSums(cw*arsd$xd)
          gyd <- Matrix::colSums(cw*arsd$yd)
          al <- sqrt(gxd^2+gyd^2);
          vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
          if(any(vg)) {
            z <- es %*% cw[,vg]
          } else { NULL }
        },mc.cores=n.cores,mc.preschedule=T))
        
        if(scale=='log') {
          nd <- (log10(abs(nd)+1)*sign(nd))
        } else if(scale=='sqrt') {
          nd <- (sqrt(abs(nd))*sign(nd))
        }
        cat(".");
        # velocity for the grid
        gv <- do.call(cbind,parallel::mclapply(gx,function(x) {
          # cell distances (rows:cells, columns: grid points)
          cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
          cw <- dnorm(cd,sd=grid.sd)
          # calculate x and y delta expectations
          gw <- Matrix::colSums(cw)
          cws <- pmax(1,Matrix::colSums(cw));
          cw <- t(t(cw)/cws)
          gxd <- Matrix::colSums(cw*arsd$xd)
          gyd <- Matrix::colSums(cw*arsd$yd)
          al <- sqrt(gxd^2+gyd^2);
          vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
          if(any(vg)) {
            z <- nd %*% cw[,vg]
          } else { NULL }
        },mc.cores=n.cores,mc.preschedule=T))
        cat(". done\n")
        
        
        return(invisible(list(tp=tp,cc=cc,garrows=garrows,arrows=as.matrix(ars),vel=nd,eshifts=es,gvel=gv,geshifts=gs,scale=scale)))
      }
      
      
      
      
    } else { # draw individual arrows
      # calculate arrows, draw
      # lapply(1:nrow(emb),function(i) {
      #   # normalized directions to each point
      #   di <- t(t(emb)-emb[i,])
      #   di <- di/sqrt(Matrix::rowSums(di^2))*arrow.scale; di[i,] <- 0;
      #   di <- Matrix::colSums(di*tp[,i]) - Matrix::colSums(di*(tp[,i]>0)/sum(tp[,i]>0)); # relative to expected kNN center
      #   
      #   if(fixed.arrow.length) {
      #     suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+di[1],emb[colnames(em)[i],2]+di[2],length=0.05,lwd=arrow.lwd))
      #   } else {
      #     ali <- sqrt( (di[1] * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + (di[2]*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2)
      #     suppressWarnings(arrows(emb[colnames(em)[i],1],emb[colnames(em)[i],2],emb[colnames(em)[i],1]+di[1],emb[colnames(em)[i],2]+di[2],length=min(0.05,ali),lwd=arrow.lwd))
      #   }
      # })
      
      apply(ars,1,function(x) {
        if(fixed.arrow.length) {
          suppressWarnings(arrows(x[1],x[2],x[3],x[4],length=0.05,lwd=arrow.lwd))
        } else {
          ali <- sqrt( ((x[3]-x[1]) * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + ((x[4]-x[2])*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2)
          suppressWarnings(arrows(x[1],x[2],x[3],x[4],length=min(0.05,ali),lwd=arrow.lwd))
        }
      })
      
      
    }
  }
  return(invisible(list(tp=tp,cc=cc)))
}