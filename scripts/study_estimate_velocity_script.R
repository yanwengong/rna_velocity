## This is to study the gene.relative.velocity.estimates function in RNA velocity package

gene.relative.velocity.estimates_mod <- function(emat,nmat ,deltaT=1,smat=NULL,steady.state.cells=colnames(emat),kCells=10,
                                             cellKNN=NULL,kGenes=1,old.fit=NULL,mult=1e3,min.nmat.smat.correlation=0.05,
                                             min.nmat.emat.correlation=0.05, min.nmat.emat.slope=0.05, zero.offset=FALSE,deltaT2=1, 
                                             fit.quantile=NULL, diagonal.quantiles=FALSE, show.gene=NULL, do.par=TRUE, 
                                             cell.dist=NULL, emat.size=NULL, nmat.size=NULL, cell.emb=NULL, cell.colors=NULL, 
                                             expression.gradient=NULL,residual.gradient=NULL, n.cores=1, verbose=TRUE) {
  
  if(!all(colnames(emat)==colnames(nmat))) stop("emat and nmat must have the same columns (cells)");
  if(!is.null(smat)) { if(!all(colnames(emat)==colnames(smat))) stop("smat must have the same columns (cells) as emat") }
  ## create an empty list 
  resl <- list();
  # bring matrices to the same gene set (just in case)
  ## vg are the name of common genes between emat and nmat
  ## start1 - to make emat and nmat share the same gene
  vg <- intersect(rownames(emat),rownames(nmat));
  if(is.null(smat)) {
    ## now, the emat and nmat share the same genes
    emat <- emat[vg,]; nmat <- nmat[vg,]
  } else {
    vg <- intersect(vg,rownames(smat))
    emat <- emat[vg,]; nmat <- nmat[vg,]; smat <- smat[vg,]
  }
  ## end1
  if(!is.null(show.gene)) {
    if(!show.gene %in% rownames(emat)) { stop(paste("gene",show.gene,"is not present in the filtered expression matrices")) }
  }
  # TODO: add gene filtering options
  pcount <- 1;
  
  
  ## start2 - this is to check on cell dist which will be used to determine cell group with knnCell 
  if(!is.null(cell.dist)) {
    if(class(cell.dist)!='dist') { stop("cell.dist must be of a class dist") }
    if(!all(labels(cell.dist)==colnames(emat))) {
      cat("matching cells between cell.dist and emat/nmat ... ")
      cell.dist <- as.matrix(cell.dist)
      cn <- intersect(colnames(emat),colnames(cell.dist))
      cell.dist <- as.dist(cell.dist[cn,cn]);
      emat <- emat[,cn]; nmat <- nmat[,cn];
      if(!is.null(smat)) { smat <- smat[,cn] }
      cat("done\n")
    }
  }
  ## end2
  
  ## start3 - transfer emat and nmat, based on cell size
  ## mult is library scale factor, default is 10e3
  if(is.null(emat.size)) { emat.size <- Matrix::colSums(emat); } ## cell size is total number of expressed genes, named with cell
  if(is.null(nmat.size)) { nmat.size <- Matrix::colSums(nmat); }
  emat.cs <- emat.size[colnames(emat)]/mult; ## emat.cs is emat.size devided by mult(10e3)
  nmat.cs <- nmat.size[colnames(nmat)]/mult;
  
  
  ## emat.log.norm is used in knnCells; normalize by emat.cs then take log 
  ## default of pcount is 1
  emat.log.norm <- log(as.matrix(t(t(emat)/emat.cs))+pcount);
  ## end3
  
  if(!is.null(old.fit)) { cellKNN <- old.fit[['cellKNN']]}
  knn.maxl <- 1e2
  ## start4 - work on knnCells 
  ## called function balancedKNN, which further called _velocyto_R_balanced_knn
  ## https://github.com/velocyto-team/velocyto.R/blob/master/src/RcppExports.cpp
  ## here i assigned the knn I have to this variable resl$cellKNN <- rvel.cd$cellKNN
  ## cellKNN <- rvel.cd$cellKNN
  if(kCells>1) {
    ## here, to see whether there is a pre-calculated cellKNN matrix
    if(is.null(cellKNN)) {
      cat("calculating cell knn ... ")
      if(is.null(cell.dist)) {
        cellKNN <- balancedKNN(emat.log.norm,kCells,kCells*knn.maxl,n.threads=n.cores);
      } else {
        ## we will input cell distance matrix normally 
        cellKNN <- balancedKNN(emat.log.norm,kCells,kCells*knn.maxl,n.threads=n.cores,dist=cell.dist);
      }
      ## make diagnol to 1
      diag(cellKNN) <- 1;
      ## cellKNN is #cell * #cell matrix, diagnol is 1, some cell with the other cell is 1
      resl$cellKNN <- cellKNN;
      cat("done\n")
    }
    rm(emat.log.norm);
    # smoothed matrices
    cat("calculating convolved matrices ... ")
    ## %*% is matrix multiplication  
    ## not sure what exactly is going on here
    ## conv is after processed with KNN
    conv.emat <- emat %*% cellKNN[colnames(emat),colnames(emat)]
    conv.nmat <- nmat %*% cellKNN[colnames(nmat),colnames(nmat)]
    conv.emat.cs <- (emat.cs %*% cellKNN[colnames(emat),colnames(emat)])[1,]
    conv.nmat.cs <- (nmat.cs %*% cellKNN[colnames(nmat),colnames(nmat)])[1,]
    cat("done\n")
  } else { 
    conv.emat <- emat; conv.nmat <- nmat; cellKNN <- NULL;
    conv.emat.cs <- emat.cs; conv.nmat.cs <- nmat.cs;
  }
  
  #browser()
  
  # size-normalized counts
  ## normalize the KNNized matrix with cell size 
  conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
  conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)
  
  # size-normalized counts
  emat.norm <- t(t(emat)/emat.cs)
  nmat.norm <- t(t(nmat)/nmat.cs)
  ## end4
  
  ## start5 - can be skipped
  ## this is to work on knnGene, can skip because default is knnGene = 1;
  if(kGenes>1) {
    if(!is.null(old.fit) && !is.null(old.fit$geneKNN)) {
      geneKNN <- old.fit$geneKNN; 
    } else {
      cat("gene kNN ... ")
      geneKNN <- balancedKNN(t(log(as.matrix(conv.emat.norm)+pcount)),kGenes,kGenes*1.2e3,n.threads=n.cores); diag(geneKNN) <- 1;
    }
    resl$geneKNN <- geneKNN;
    
    
    # normalize contribution of different neighbor genes to match the median totals (to avoid distortions due to high-yielding genes)
    cat("scaling gene weights ... ")
    gt <- rowSums(conv.emat.norm)
    scaledGeneKNN <- t(apply(geneKNN,2,function(ii) pmin(1,median(gt[which(ii>0)])/gt) * ii))
    cat("convolving matrices ... ")
    conv.emat.norm <- scaledGeneKNN %*% conv.emat.norm;
    conv.nmat.norm <- scaledGeneKNN %*% conv.nmat.norm;
    
    cat("done\n")
  }
  
  
  ## can be skipped
  if(!is.null(smat)) {
    
    if(kCells>1) {
      conv.smat <- smat %*% cellKNN[colnames(smat),colnames(smat)]
    } else {
      conv.smat <- smat
    }
    conv.smat.cs <- Matrix::colSums(conv.smat)/mult;
    conv.smat.norm <- t(t(conv.smat)/conv.smat.cs)
    
    if(kGenes>1) {
      conv.smat.norm <- scaledGeneKNN %*% conv.smat.norm;
    } 
    
    # use spanning reads to fit offset for the intronic reads, test correlation
    if(is.null(old.fit)) {
      cat("fitting smat-based offsets ... ")
      sfit <- data.frame(do.call(rbind,parallel::mclapply(sn(rownames(conv.emat.norm)),function(gn) {
        df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]),s=conv.smat.norm[gn,steady.state.cells])
        sd <- lm(n~s,data=df)
        r <- with(df[df$s>0,],cor(n,s,method='spearman'),3)
        return(c(o=pmax(0,as.numeric(sd$coef[1])),s=as.numeric(sd$coef[2]),r=r))
      },mc.cores=n.cores,mc.preschedule=T)))
      cat("done\n")
      
    } else {
      sfit <- old.fit$sfit;
    }
  }
  
  ## end5
  
  ## FINISHED normalization (size norm, log, knn, etc.)
  ## store normalized matrix of spiced and unspliced gene 
  resl$conv.nmat.norm <- conv.nmat.norm;
  resl$conv.emat.norm <- conv.emat.norm;
  
  # fit gamma, using the offset above
  ## skip: default for show.gene is null; can skip to start6: cat("fitting gamma coefficients ... ")
  if(!is.null(show.gene)) {
    gn <- show.gene;
    if(!is.null(cell.emb)) {
      # show embedding heatmaps
      cc <- intersect(rownames(cell.emb),colnames(conv.emat.norm));
      if(do.par) { par(mfrow=c(1,4), mar = c(2.5,2.5,2.5,0.5), mgp = c(1.5,0.65,0), cex = 0.85); }
      plot(cell.emb[cc,],pch=21,col=ac(1,alpha=0.2),bg=val2col(conv.emat.norm[gn,cc],gradientPalette=expression.gradient),cex=0.8,xlab='',ylab='',main=paste(gn,'s'),axes=F); box();
      plot(cell.emb[cc,],pch=21,col=ac(1,alpha=0.2),bg=val2col(conv.nmat.norm[gn,cc],gradientPalette=expression.gradient),cex=0.8,xlab='',ylab='',main=paste(gn,'u'),axes=F); box();
    }
    do <- NULL;
    if(!is.null(smat)) { # use smat-based offsets
      df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]),o=sfit[gn,'o'])
      if(zero.offset) df$o <- 0;
    }  else { # calculate offset based on the nascent counts obsered for near-0 exonic levels
      df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]))
      o <- 0;
      df$o <- o;
      #zi <- emat[gn,steady.state.cells]==0;
      if(!zero.offset) { zi <- df$e<1/conv.emat.cs[steady.state.cells]; if(any(zi)) { o <- sum(df$n[zi])/(sum(zi)+1)} }
      df$o <- o;
      
      #table(zi)
      # if(any(zi)) {
      #   do <- lm(n~e,data=df[zi,])
      #   #summary(do)
      #   df$o <- max(0,do$coefficients[1])
      # }
      
      
    }
    
    
    #browser()
    d <- lm(n~e+offset(o)+0,data=df,weights=df$e^4+df$n^4);
    cell.col <- ac(rep(1,nrow(df)),alpha=0.1); names(cell.col) <- rownames(df)
    if(!is.null(cell.colors)) { 
      cc <- intersect(names(cell.colors),rownames(df)); 
      cell.col[cc] <- cell.colors[cc]
    }
    plot(df$e,df$n,pch=21,bg=ac(cell.col,alpha=0.3),col=ac(1,alpha=0.1),cex=0.8,xlab='s',ylab='u',main=paste(gn,'fit'))
    if(!is.null(do)) {
      abline(do,lty=2,col=8)
    }
    
    # min/max fit
    if(!is.null(fit.quantile)) {
      if(diagonal.quantiles) {
        # determine maximum ranges 
        emax <- quantile(df$e,p=0.99)
        nmax <- quantile(df$n,p=0.99)
        if(emax==0) emax <- max(max(df$e),1e-3)
        if(nmax==0) nmax <- max(max(df$n),1e-3)
        x <- df$e/emax + df$n/nmax;
        eq <- quantile(x,p=c(fit.quantile,1-fit.quantile))
        if(!is.null(smat)) { # will use smat offset, so disregard lower quantile
          pw <- as.numeric(x>=eq[2])
        } else {
          pw <- as.numeric(x>=eq[2] | x<=eq[1])
        }
      } else {
        eq <- quantile(df$e,p=c(fit.quantile,1-fit.quantile))
        if(!is.null(smat) || zero.offset) { # will use smat offset, so disregard lower quantile
          pw <- as.numeric(df$e>=eq[2])
        } else {
          pw <- as.numeric(df$e>=eq[2] | df$e<=eq[1])
        }
      }
      
      if(!is.null(smat) || zero.offset) { # use smat offset
        d <- lm(n~e+offset(o)+0,data=df,weights=pw);
      } else {
        d <- lm(n~e,data=df,weights=pw)
      }
      
      ## eq <- quantile(df$e,p=c(fit.quantile,1-fit.quantile))
      ## pw <- as.numeric(df$e>=eq[2] | df$e<=eq[1])
      ## if(!is.null(smat)) { # use smat offset
      ##   d <- lm(n~e+offset(o)+0,data=df,weights=pw);
      ## } else {
      ##   d <- lm(n~e,data=df,weights=pw)
      ## }
    } 
    
    
    df <- df[order(df$e,decreasing=T),]; 
    lines(df$e,predict(d,newdata=df),lty=2,col=2)
    
    
    if(!is.null(cell.emb)) {
      plot(cell.emb[cc,],pch=21,col=ac(1,alpha=0.2),bg=val2col(resid(d)[cc],gradientPalette=residual.gradient),cex=0.8,xlab='',ylab='',main=paste(gn,'resid'),axes=F); box();
    }
    if(kGenes>1) { return(invisible(geneKNN)) } else { return(1) }
  }
  
  ## start6
  cat("fitting gamma coefficients ... ")
  if(is.null(old.fit)) {
    print("K=0.5, n=3")
    ## for each gene (row), all cells (all column; steady.state.cells(IT IS ALL THE CELLS??)) apply the function (gn is input), 
    ## then combine all the rows with do.call(rbind)
    ## note: # quick self-naming vector routine: sn <- function(x) { names(x) <- x; x}
    ## "sn(rownames(conv.emat.norm))" get the list of genes
    ko <- data.frame(do.call(rbind,parallel::mclapply(sn(rownames(conv.emat.norm)),function(gn) {
      ## skip
      if(!is.null(smat)) { # use smat-based offsets
        # steady.state.cells=colnames(emat), default of steady.state.cells is all cells
        df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]),o=sfit[gn,'o'])
        if(zero.offset) df$o <- 0;
      }  
      ## start
      else { # calculate offset based on the nascent counts obsered for near-0 exonic levels
        ## df is the df of normalized counts of spliced and unspliced of one gene with all cells
        df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),e=(conv.emat.norm[gn,steady.state.cells]))
        o <- 0;
        ## default of zero.offset is FALSE
        ## zi should be a list of true or false
        #### I comment the line below temperorily, because I cannot get conv.emat.cs ###
        #if(!zero.offset) { zi <- df$e<1/conv.emat.cs[steady.state.cells]; if(any(zi)) { o <- sum(df$n[zi])/(sum(zi)+1)} }
        df$o <- o;
      }
      if(is.null(fit.quantile)) {
        #d <- lm(n~e+offset(o)+0,data=df,weights=df$e^4+df$n^4);
        d <- lm(n~e+offset(o)+0,data=df,weights=df$e^4+df$n^4);
        return(c(o=df$o[1],g=as.numeric(coef(d)[1]),r=cor(df$e,df$n,method='spearman')))
      } else {
        ## default diagonal.quantiles=FALSE, so skip
        ##diagonal.quantiles whether extreme quantiles should be computed diagonally
        if(diagonal.quantiles) {
          # determine maximum ranges 
          emax <- quantile(df$e,p=0.99)
          nmax <- quantile(df$n,p=0.99)
          if(emax==0) emax <- max(max(df$e),1e-3)
          if(nmax==0) nmax <- max(max(df$n),1e-3)
          x <- df$e/emax + df$n/nmax;
          eq <- quantile(x,p=c(fit.quantile,1-fit.quantile))
          if(!is.null(smat)) { # will use smat offset, so disregard lower quantile
            pw <- as.numeric(x>=eq[2])
          } else {
            pw <- as.numeric(x>=eq[2] | x<=eq[1])
          }
        } ## skip ends 
        else {
          ## e is the normalized spliced counts
          ## eq is the two numbers corresponds to quantile  
          eq <- quantile(df$e,p=c(fit.quantile,1-fit.quantile))
          ## default: smat = null, zero.offset = FALSE
          if(!is.null(smat) || zero.offset) { # will use smat offset, so disregard lower quantile
            ## pw is o or 1, indicates whether included in the quantiles
            pw <- as.numeric(df$e>=eq[2])
          } else {
            ## pw is 0 or 1, indicates whether included in the quantiles
            ## for us, will include both higher and lower quantile as below
            pw <- as.numeric(df$e>=eq[2] | df$e<=eq[1])
          }
        }
        
        if(!is.null(smat) || zero.offset) { # use smat offset
          d <- lm(n~e+offset(o)+0,data=df,weights=pw);
          return(c(o=df$o[1],g=as.numeric(coef(d)[1]),r=cor(df$e,df$n,method='spearman')))
        } else {
          ## THIS IS THE LINE OF FITTING A LINEAR MODEL!! 
          ## for each gene across all the cells, fit linear model of n~e
          #print("origin:")
          #print(df$n[1:5])
          
          replaced_n=(df$n)^3/(0.5^3+(df$n)^3)
          #print("replaced")
          #print(replaced_n[1:5])
          df$n=replaced_n
          #print("replaced in df")
          #print(df$n[1:5])
          d <- lm(n~e,data=df,weights=pw)
          # note: re-estimating offset here
          ## o is the estimated intercept (offset), g is the gamma , r is the spearman correlation between spliced and uncpliced for this gene
          return(c(o=as.numeric(coef(d)[1]),g=as.numeric(coef(d)[2]),r=cor(df$e,df$n,method='spearman')))
        }
      }
      
    },mc.cores=n.cores,mc.preschedule=T)))
    ## ko is a data from, column name is: o (intercept), g(gamma), r(correlation); 
    ## each row is one gene
    ko <- na.omit(ko)
    cat("done. succesfful fit for",nrow(ko),"genes\n")
  } else { full.ko <- ko <- na.omit(old.fit$ko); }
  
  ## skip
  if(!is.null(smat)) {
    sfit <- na.omit(sfit)
    ko <- ko[rownames(ko) %in% rownames(sfit),]; # omit genes for which sfit didn't work
    vi <- sfit$r > min.nmat.smat.correlation
    ko <- ko[vi,]
    if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to low nmat-smat correlation\n")
  }
  
  ## start
  full.ko <- ko;
  ## min.nmat.emat.correlation=0.05 defalt
  ## filter out genes is correlation between spliced and unspliced is low.
  vi <- ko$r>min.nmat.emat.correlation
  if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to low nmat-emat correlation\n")
  ko <- ko[vi,]
  
  ##  min.nmat.emat.slope=0.05 default
  ## filter genes out if gamma is low
  vi <- ko$g>min.nmat.emat.slope
  if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to low nmat-emat slope\n")
  ko <- ko[vi,]
  
  ## each gene has one gamma, one offset
  gamma <- ko$g; offset <- ko$o; names(gamma) <- names(offset) <- rownames(ko);
  ## end6
  
  
  ## start7
  cat("calculating RNA velocity shift ... ")
  ## default kGenes is 1, so skip
  if(kGenes>1) { # gene-convolved estimation
    # estimate M value
    ## for each cell, calculate predicted U, by u = gamma*s + offset
    npred <- gamma*conv.emat.norm[names(gamma),] + ko$o;
    npred[npred<0] <- 0;
    ## mvel is log2[(normalized unsplice +pcount)/(predicted normalized unsplice + pcount)
    mval <- log2(conv.nmat.norm[names(gamma),]+pcount) - log2(npred+pcount);
    resl$mval <- mval;
    
    #resl$conv.deltaE <- t.get.projected.delta(conv.emat.norm,conv.nmat.norm,gamma,offset=offset,delta=deltaT)
    #resl$conv.projected <- t.get.projected.cell2(conv.emat.norm,emat.size,as.matrix(resl$conv.deltaE),mult = mult,delta=deltaT2);
    #resl$conv.projected[resl$conv.projected<0] <- 0;
    
    # switch back to non-gene-kNN conv.* matrices
    conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
    conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)
    # estimate gamma
    cat("re-estimating gamma of individual genes ... ")
    
    am <- conv.nmat.norm[rownames(mval),]-offset; am[am<0] <- 0;
    fm <- log2(am) - mval - log2(conv.emat.norm[rownames(mval),])
    wm <- is.finite(fm)
    fm[!is.finite(fm)] <- 0;
    gammaA <- 2^(rowSums(fm * wm)/rowSums(wm))
    gammaA <- gammaA[is.finite(gammaA)];
    
    
    gamma <- gammaA;
    cat("done\n")
    
    
    # can estimate deltaE from the mval
    cat("calculating RNA velocity shift ... ")  
    # estimate delta from M value
    deltaE <- t.get.projected.delta.from.log2ratio(em=conv.emat.norm,gamma=gamma,r=mval,delta=deltaT)
    #deltaE <- t.get.projected.delta2(conv.emat.norm,conv.nmat,conv.nmat.cs,gamma,offset=offset,delta=deltaT)
  } 
  ## skip done
  
  ## start
  else { # regular estimation
    ## t.get.projected.delta2 is another projected method 
    #deltaE <- t.get.projected.delta2(conv.emat.norm,conv.nmat,conv.nmat.cs,gamma,offset=offset,delta=deltaT)
    # estimate projected delta given x'=(y-o) - gamma*x solution
    ## deltaE is the diff expected spliced vs observed splice, by assuming Ut=Uo:
    ## St = So*exp(-gamma*t) + Uo/gamma*(1-exp(-gamma*t))
    ## default deltaT=1
    deltaE <- t.get.projected.delta(conv.emat.norm,conv.nmat.norm,gamma,offset=offset,delta=deltaT)
  }
  
  ## end7
  
  ## save gamma
  resl$gamma <- gamma;
  
  cat("done\n")
  cat("calculating extrapolated cell state ... ")
  
  ## start8 - celculate projected spliced, normalize current spliced with size:
  ##projected=emn,current=emat.norm
  # reduced cell normalization (only genes for which momentum was estimated)
  emat.norm <- emat[rownames(emat) %in% rownames(deltaE),]
  #emat.sz <- Matrix::colSums(emat.norm)/mult;
  #browser()
  emat.sz <- emat.cs;
  emat.norm <- t(t(emat.norm)/(emat.sz));
  
  ## calculates the difference in the number of counts based on the library size, renormalizes
  emn <- t.get.projected.cell2(emat.norm,emat.sz,as.matrix(deltaE),mult = mult,delta=deltaT2);
  #emn <- t.get.projected.cell(emat.norm,as.matrix(deltaE),target.mult = mult,model.mult=mult,delta=deltaT2,size.normalize=FALSE);
  
  #table(emat.norm[,cn]==0)
  #table(emn[,cn]==0)
  
  cat("done\n")
  full.ko$valid <- rownames(full.ko) %in% rownames(ko)
  resl <- c(resl,list(projected=emn,current=emat.norm,deltaE=deltaE,deltaT=deltaT,ko=full.ko,mult=mult,kCells=kCells));
  if(!is.null(smat)) { resl$sfit <- sfit }
  return(resl)
}



## used functions

balancedKNN <- function(val,k,maxl=k,return.distance.values=FALSE,n.threads=1,dist='cor') {
  if(class(dist)=="dist") { # actual distance was passed
    if(!all(labels(dist)==colnames(val))) { stop("balancedKNN(): supplied distance doesn't match the columns of val") }
    cd <- as.matrix(dist);
  }  else {
    if(dist=='cor') {
      cd <- 1-cor(val);
    } else if(dist=='euclidean') {
      cd <- as.matrix(dist(t(val)))
    } else {
      stop(paste("unknown distance",dist,"specified"))
    }
  }
  z <-  balanced_knn(cd,k,maxl,return.distance.values,n.threads);
  rownames(z) <- colnames(z) <- colnames(val);
  z
}

# quick self-naming vector routine
sn <- function(x) { names(x) <- x; x}

#' adjust colors, while keeping the vector names
#' 
#' @param x color vector
#' @param alpha transparenscy value (passed to adjustcolors as alpha.f)
#' @param ... parameters passsed to adjustcolor
#' @export
ac <- function(x, alpha=1, ...) { y <- adjustcolor(x, alpha.f=alpha, ...); names(y) <- names(x); return(y)}

# quick function to map value vector to colors
val2col <- function(x,gradientPalette=NULL,zlim=NULL,gradient.range.quantile=0.95) {
  if(all(sign(x)>=0)) {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- as.numeric(quantile(na.omit(x),p=c(1-gradient.range.quantile,gradient.range.quantile)))
      if(diff(zlim)==0) {
        zlim <- as.numeric(range(na.omit(x)))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])
    
  } else {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- c(-1,1)*as.numeric(quantile(na.omit(abs(x)),p=gradient.range.quantile))
      if(diff(zlim)==0) {
        zlim <- c(-1,1)*as.numeric(na.omit(max(abs(x))))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])
    
  }
  gp <- gradientPalette[x*(length(gradientPalette)-1)+1]
  if(!is.null(names(x))) { names(gp) <- names(x) }
  gp
}

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
  ########################## APPLY HILL FUNCTION HILL#################################
  N=3; K=0.5;
  nm_hill <- nm^N/(nm^N+K^N)
  ## y is 1116(genes) x 3021 cells for unwouded 1
  y <- nm_hill-offset; y[y<0] <- 0; # zero out entries with a negative n levels after offset adjustment
  
 
  
  em*egt + (1-egt)*y/gamma  - em ## each gene, each cell have a projected.delta
}

# estimate projected delta given log2 fold observed/expected nascent ratio 
# em - normalized expression matrix
# nm - normalized nascent matrix
# gamma - inferred degradation coefficients
# r - log2 observed/expected nascent cont ratio
# delta - time to project forward
t.get.projected.delta.from.log2ratio <- function(em,gamma,r,delta=0.5,min.val=1e-4) {
  # adjust rownames
  gn <- intersect(intersect(names(gamma),rownames(em)),rownames(r));
  em <- em[gn,]; gamma <- gamma[gn]; r <- 2^r[gn,];
  # time effect constant
  egt <- exp(-gamma*delta);
  (em+min.val)*(egt*(1-r) +r) - em
}


# calculates the difference in the number of counts based on the library size, renormalizes
# note: also introduces centripetal velocity
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

