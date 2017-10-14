#' A script wrapper
#'
#' This function wraps toghether the script to perform the SRSF shape analysis of reads on specified genomic regions
#' @param filenames names of the UCSC bed files containing read coordinates
#' @param TSSregions name of the UCSC bed file specifying genomic regions for the shape analysis (e.g. +/- 1000bp intervals from TSS)
#' @param MARG Default is 0. Specifies by how many bp should the regions (specified in TSSregions) be expended by. Allows negative numbers
#' @param CUTOFF Default is 100. Specifies minimum number of reads per region of interest for the region to be included in the analysis. Low read counts per region could result in unreliable and noisy read density shapes.
#' @param ADJ Default is 0.5. Bandwidth adjustment for density estimation
#' @param p Default is 1000. Size of the margins of the left and right end of genomic region to account for unnecessary shifts in empty regions (where read density is zero). This is different then MARG parameter and is only used to clear possible algorithm errors.
#' @param weight.adj Default is 0. Weight of the read coordinates in estimating the read desnity over region of interest.
#' @param penalty Default is 0.0001. Controls the penalty for the amount of shift.
#' @param quant Default is 0.999 Top quantile of the regions with the highest shift statistic that are to be visualized
#' @param max.cores Default is given by detectCores() function.
#' @keywords main wrapper
#' @export
#' @examples
#' TSSregions.bed.name="data/chr1.TSSgencode.bed";
#' MARG = 0; CUTOFF=100;  ADJ=0.5;p=1000;weight.adj=0;penalty=0.0001; quant=0.999; max.cores=56;
#' group.id<-NULL #Need to adjust for variable group sizes
#' filenames=c("data/chr1.GH1_cap_R1_clip.fastq_q20.bed", "data/chr1.GH2_cap_R1_clip.fastq_q20.bed",
#'             "data/chr1.IH1_cap_R1_clip.fastq_q20.bed", "data/chr1.IH2_cap_R1_clip.fastq_q20.bed")
#' results<-SRSFN_wrapper(filenames,TSSregions)


SRSFN_wrapper<-function(filenames, TSSregions, MARG = 0, CUTOFF=100,  ADJ=0.5, p=1000, weight.adj=0, penalty=0.0001, quant=0.999, max.cores=detectCores())
{
  ##--------------------------------------------------------------------------------
  cat("1.===== Loading Datasets:\n TSSregions")
  TSSregions<-fread(TSSregions.bed.name, sep="\t", fill=TRUE)

  heavy<-mclapply(filenames, function(filename)
    fread(input=filename, sep="\t", fill=TRUE, sep2=""),
    mc.cores=min(length(filenames),max.cores), mc.preschedule=F)
  names(heavy)<-filenames

  cat(": DONE ",date(), "\n")

  ##--------------------------------------------------------------------------------
  cat("2. ===== Extracting reads per TSSregion, per sample. Margin is ", MARG ,"\n")
  cat(paste(rep("|",49), collapse=""), "    '|' =  2%  = ",round(dim(TSSregions)[1]/50),
      " TSS processed out of ",dim(TSSregions)[1],  " \n")

  TSS.reads<-mclapply(1:dim(TSSregions)[1], function(i)
    get.TSS.reads(heavy, TSSinterval=c(TSSregions$V2[i], TSSregions$V3[i]),
                  margin=MARG, i, dim(TSSregions)[1]),
    mc.cores=min(detectCores()-1, max.cores), mc.preschedule=F)
  cat(date(), " DONE \n ")

  ##--------------------------------------------------------------------------------
  ## CALCULATING THE OPTIMAL WARPING
  cat("3. ===== Calculating shifts, parameters are",ADJ,p,weight.adj,penalty,MARG,
      CUTOFF,"\n", sep=" ")
  cat(paste(rep("|",49), collapse=""), "    '|' =  2%  = ",round(dim(TSSregions)[1]/50),
      " TSS processed out of ",dim(TSSregions)[1],  " \n")

  all.gammas <- mclapply(1:length(TSS.reads), function(i)
    get.gamma(TSS.reads[[i]],i,CUTOFF, ADJ, p, weight.adj, penalty, N=dim(TSSregions)[1]),
    mc.preschedule = F, mc.cores = min(detectCores()-1, max.cores) )
  cat(date(), "DONE \n")

  ##--------------------------------------------------------------------------------
  cat("4. ===== Calculating Statistics ")
  all.gammas.trunc<-lapply(all.gammas, function(gam) Jorges.function(gam, "both", p))

  gamma.stats<-unlist(lapply(1:length(all.gammas.trunc),
                             function(TSS.id)
                               get.stats.gamma(gam=all.gammas.trunc[[TSS.id]][,1]) ))
  table.results<-cbind(TSSregions, gamma.stats)
  cat(date(), " DONE \n")


  ##--------------------------------------------------------------------------------
  ##PLOTTING
  cat("5. ===== Plotting top ", quant, " quantile")
  big.norm.id<-which(gamma.stats>quantile(gamma.stats, quant))
  dir.create(paste("Figures",ADJ,p,weight.adj,penalty,sep="_"))
  lapply(big.norm.id,
         function(i) plot.TSS.data(TSS.reads[[i]],i,gamma.results=all.gammas[[i]],p,
                                   prefix=paste(ADJ,p,weight.adj,penalty, sep="_"),
                                   name=table.results[i,4]))

  return(table.results)
}


#' Auxilary function for loading and saving files
#'
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' load.save()
load.save<- function(filename)

{

  print(paste("Loading ", filename, sep=""))
  temp<-fread(filename)
  print("Selecting relevant columns")

  if (filename %in% c("GH1_cap_R1_clip.fastq_q20.bed",
                      "GH2_cap_R1_clip.fastq_q20.bed",
                      "IH1_cap_R1_clip.fastq_q20.bed",
                      "IH2_cap_R1_clip.fastq_q20.bed"))
    temp<-cbind(temp$V1, as.numeric(temp$V2),as.numeric(temp$V3), temp$V6);

  print("Saving in binary format")
  save(temp, file=paste(filename,".Rdata",sep="" ))

  return(NULL)
}

##############################################################################
                                        #Input:
                                        #Output:
                                        #Description:

##reads.per.TSSfunction(TSSinterval, midpoints, margin)
                                        #{
                                        #  TSSID=which(midpoints<(TSSinterval[2]+margin) &
                                        #              midpoints>(TSSinterval[1]-margin))
                                        # return(midpoints[TSSID])
                                        #}

#############################################################################
                                        #Input:
                                        #Output:
                                        #Description:
load.dataset<-function(filename)
{
  load(paste(filename, ".Rdata", sep=""))
  cat(filename," DONE\n")
  return(temp)
}

##############################################################################
                                        #Input:
                                        #Output:
                                        #Description:
get.TSS.reads<-function(heavy, TSSinterval, margin, i, N)
{
  if ((i %% round(N/50))==0)
    cat("|")

  output=list()
  for (file.id in 1:length(heavy) )
  {
    midpoints=((heavy[[file.id]][,2]) +
               (heavy[[file.id]][,3]))/2
    ##    midpoints=(as.numeric(heavy[[file.id]][,2]) +
    ##              as.numeric(heavy[[file.id]][,3]))/2

    TSSID=which(midpoints<(TSSinterval[2]+margin) &
                midpoints>(TSSinterval[1]-margin))
    ##    output[[file.id]]=reads.per.TSS(TSSinterval, midpoints, margin)

    output[[file.id]]=as.numeric(unlist(midpoints[TSSID]))
  }
  names(output)<-names(heavy)

  return(output)
}

#############################################################################
                                        #Input:
                                        #Output:
                                        #Description:

get.bandwidth<-function(current.TSS, n=256, method="nrd")
{
  if (method=="nrd")
    all.bandwidth<-lapply(current.TSS, function(ct) bw.nrd(ct))
  if (method=="ucv")
    all.bandwidth<-lapply(current.TSS, function(ct) bw.ucv(ct))
  bandwidth<-mean(unlist(all.bandwidth))
  return(bandwidth)
}

get.curves<-function(current.TSS,n=256,ADJ,p, weight.adj)
{
  domain=c(min(unlist(current.TSS))-p, max(unlist(current.TSS))+p)
  bandwidth=get.bandwidth(current.TSS)

  weights=list() ##this assumes that observations come sorted
  weights.list<-lapply(1:length(current.TSS),
                       function(sample.id)
                         approx(density(current.TSS[[sample.id]],
                                        bw=bandwidth,adjust=ADJ, n=n,
                                        from = domain[1], to = domain[2])$x,
                                density(current.TSS[[sample.id]]
                                       ,bw=bandwidth,adjust=ADJ, n=n,
                                        from = domain[1], to = domain[2])$y,
                                xout = sort(current.TSS[[sample.id]]))$y)
  weights.list<-lapply(weights.list,
                       function(weights)
                       (weights+weight.adj)/sum(weights + weight.adj))


  curvearray=array(0,c(n,length(current.TSS)))
  curvearray<-sapply(1:length(current.TSS),
                     function(sample.id)
                       density(sort(current.TSS[[sample.id]]),
                               weights=weights.list[[sample.id]],
                               bw=bandwidth,adjust=ADJ, n=n,
                               from = domain[1], to = domain[2])$y)
  return(curvearray)
}

##############################################################################
                                        #Input:
                                        #Output:
                                        #Description:

get.gamma<-function(current.TSS,i, cutoff , ADJ, p, weight.adj,penalty,N,n=256 )
{
  if ((i %% round(N/50))==0)
    cat("|")

  domain=c(min(unlist(current.TSS))-p, max(unlist(current.TSS))+p)
  time=seq(0,1,length.out = n)
  results=array(0,c(1,11))
  gamma<-seq(0,1,length.out=n)
  if(min(unlist(lapply(current.TSS, function(ct) length(ct))))>cutoff)
  {

    ##------------------------------------------------------------------------
    ##Obtaining densities:
    curvearray4<-get.curves(current.TSS,n, ADJ=ADJ, p=p, weight.adj=weight.adj)

    ##------------------------------------------------------------------------
    ##Obtaining aligned means (is SRVF space):
    sink("nowhere")
    ##The sink function prevents time_warping form dumping text to the terminal,
    ##just to keep it clear
    n=length(curvearray4[,1]);
    time=seq(0,1, length.out = n)
    all1<-time_warping(curvearray4[,1:2],
                       time=time, MaxItr=10, showplot = FALSE, lambda=penalty)
    all2<-time_warping(curvearray4[,3:4],
                       time=time, MaxItr=10, showplot = FALSE, lambda=penalty)
    f1=all1$fmean
    f2=all2$fmean
    q1=all1$mqn
    q2=all2$mqn
    sink()
    ##Restores the regular output in the
    ##terminal, stopped by previous sink
    ##------------------------------------------------------------------------
    ##Obtaining the shift between q1, q2
    gam=optimum.reparam(q1,time,q2,time, lambda = penalty)
    f2.gaminv<-approx(time,f2,xout = gam)$y
    gam<-approx(gam,time,xout = time)$y

                                        #Obtaining aligned f2

    gamma<-gam*diff(domain)+domain[1]

    results<-cbind(gam, gamma, f1, f2, q1, q2, f2.gaminv, curvearray4)
  }
  ##------------------------------------------------------------------------
  return(results)
}

################################################################################
                                        #Input:
                                        #Output:
                                        #Description:
get.stats.gamma<-function(gam, f1=0, f2=0)
{
                                        #  f2.gaminv<-approx(time,f2,xout = gam)$y
  dev.norm=0
  len<-length(gam)
  if (len>1)
  {
    dev.norm<-sum((1-sqrt(diff(gam)*length(gam)))^2)
    dev.max<-max(sqrt(diff(gam)*256))
    gam.max<-max(abs(gam-seq(0,1,length.out = length(gam)) ) )
  }
  return(dev.norm)
}


################################################################################
##Input:
##Output:
##Description:
plot.TSS.data<-function(current.TSS,i,gamma.results,p, CEX=1, prefix, name)
{
  domain=c(min(unlist(current.TSS))-p, max(unlist(current.TSS))+p)
  n=dim(gamma.results)[1];
  genomic.region=seq(domain[1], domain[2], length.out = n)
  time=seq(0,1,length.out = n)

  curvearray4<-gamma.results[,8:11]
  f1=gamma.results[,3]; f2=gamma.results[,4]
  f2.gaminv=gamma.results[,7]
  gam=gamma.results[,2]
  gamma=gamma.results[,1]
  ##-----------------------------------------------------------------------------
  ##Reads (points over genome)
  setEPS()
  postscript(paste("Figures_",prefix,"/TSS_",i,".eps",sep=""),
             width=5, height=9)
  ## layout(matrix(c(1,2,2,3,3,4,4,5,5),9, byrow=T) )
  layout(matrix(c(1,1,2,2,2,3,3,3,4,4,4),11, byrow=T) )
  plot(current.TSS[[1]], rep(1,length(current.TSS[[1]])),
       ylim=c(-1-length(current.TSS),0),
       xlim=domain, cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,
       main=paste("TSS of gene ",name,": Mapped reads", sep=""),
       xlab="Reference genome position",
       ylab="Sample index")
  lapply(1:length(current.TSS), function(sample.id)
    points(current.TSS[[sample.id]],
           rep(-sample.id,length(current.TSS[[sample.id]]))+
           rnorm(length(current.TSS[[sample.id]]), 0, 0.05),
           col=ceiling(sample.id/2)))
  abline(v=c(domain[1]+p, domain[2]-p), lty=2, lwd=1)
  ##----------------------------------------------------------------------------
  ##Densities
                                        #  curvearray4<-get.curves(current.TSS,n, ADJ=ADJ, p=p, weight.adj=weight.adj)

  plot(genomic.region, curvearray4[,1],
       xlim=domain, ylim=c(0,max(curvearray4)),col=0,
       cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,
       main="Read densities",
       xlab="Reference genome position",
       ylab="Density")


  lapply(1:dim(curvearray4)[2], function(sample.id)
    lines(genomic.region, curvearray4[,sample.id],
          col=ceiling(sample.id/2), lwd=2))
  abline(v=c(domain[1]+p, domain[2]-p), lty=2, lwd=1)
  ##----------------------------------------------------------------------------
  ##Aligned means

  plot(genomic.region,f1, type="l", lwd=2,
       cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,
       xlim=domain, ylim=c(0,max(curvearray4)),
       main="Averaged vs Aligned read desnities",
       xlab="Reference genome position",
       ylab="Density")

  lines(seq(domain[1], domain[2],length.out = length(f2)),f2,
        type="l", col=2, lwd=2)
  lines(genomic.region, f2.gaminv, lwd=2, lty=2, col=2 )
  abline(v=c(domain[1]+p, domain[2]-p), lty=2, lwd=1)
  ##-----------------------------------------------------------------------------
  ## the shifting (gam)
  plot(seq(domain[1], domain[2], length.out = length(diff(gam))),
       rev(diff(gam)*n/diff(domain)), type="l", col=2, lwd=2,
       cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,
       xlim=domain,
       main="Shift necessary for alignment",
       xlab="Reference genome position",
       ylab="Local shift size")

  abline(h=1, lwd=2,cex.lab=CEX)
  abline(v=c(domain[1]+p, domain[2]-p), lty=2, lwd=1)

  ##  plot(seq(domain[1], domain[2], length.out = length(gam)),
  ##       gam-seq(domain[1], domain[2], length.out = length(gam)),
  ##       cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,
  ##       type="l", col=2, lwd=3, ylim=c(-diff(domain), diff(domain)));
  ##  abline(a=0, b=0, lwd=2)
  ##  abline(a=-domain[1], b=1, lty=2, lwd=0.5)
  ##  abline(a=domain[1], b=-1, lty=2, lwd=0.5)
  ##  abline(a=domain[2], b=-1, lty=2, lwd=2)
  ##  abline(a=-domain[2], b=1, lty=2, lwd=2)

  dev.off()
}


##############################################################################
##Input:
##Output:
##Description:
Jorges.function<-function(gam,method,p)
{
  if (dim(gam)[1]<10)
    truncated.gamma<-gam
  if (dim(gam)[1]>9)
  {
    gamma<-gam[,2]
    cutoff.left=min(gamma)+p
    cutoff.right=max(gamma)-p
    
    p.left=c(min(which(gamma > cutoff.left)))-1
    p.right=c(max(which(gamma < cutoff.right)))
    if (method=="both")
      truncated.gamma <- gam[p.left:p.right,]
    if (method=="left")
      truncated.gamma <- gam[p.left:dim(gam)[1],]
    if (method=="right")
      truncated.gamma <- gam[1:p.right,]
    if (method=="none")
      truncated.gamma <- gam
  }
  return(truncated.gamma)
}

