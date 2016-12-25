oathPlot<-function(filename=NULL,n=0,alpha=0.05,outfile="oathPlot_res.pdf"){
  par(ps=12)
  pdf(paste0(outfile), paper="a4")
  
  #rawDat=read.table(arg[1], as.is=T)
  
  n=156744
  rawDat=read.table(filename, as.is=T) 
  
  pThreshold1=-log10(alpha/as.numeric(n))
  
  dimC=1
  while(TRUE)
  {
    if((2^dimC-1)==nrow(rawDat)) break
    dimC=dimC+1
  }
  
  modCnt=0
  for(i in 1:dimC)
  {
    modCnt = modCnt+ncol(combn(dimC, i))
  }
  
  mat=matrix(0, modCnt, dimC)
  cnt=1
  for(i in 1:dimC)
  {
    cb=combn(dimC, i)
    for(j in 1:ncol(cb))
    {
      mat[cnt,cb[,j]] = 1
      cnt=cnt+1
    }
  }
  
  sortDat=rawDat[sort(rawDat$V8, index.return = TRUE)$ix,]
  mode=matrix(0, modCnt, dimC)
  Str1=unlist(strsplit(sortDat$V1, ".oath"))
  COL=c(1:dimC)
  for(i in 1:modCnt)
  {
    subStr1=unlist(strsplit(Str1[i*2-1], "-"))
    num=as.numeric(subStr1[length(subStr1)])
    mode[i,]=mat[num,]
  }
  
  ##########snp
  par(mai=c(0.2, 1, 0.2, 0.5))
  mat<-matrix(c(rep(1,dimC*3),rep(2,dimC*3),rep(3,dimC*3),rep(4:(dimC+3)),(dimC+4),(dimC+4),(dimC+4)),ncol=1)
  layout(mat)
  barplot(-log10(pnorm(-abs(sortDat$V8))*2), border = F,ylim =c(0,max(-log10(pnorm(-abs(sortDat$V8))*2),pThreshold1)+1), ylab=expression(paste(-log[10], "(p)")), col=ifelse(-log10(pnorm(-abs(sortDat$V8))*2) > pThreshold1, "grey30", "grey70"))
  abline(h=pThreshold1, lty=2)
  barplot(sortDat$V6, border=F, ylab=expression(beta), col=ifelse(-log10(pnorm(-abs(sortDat$V8))*2) > pThreshold1, "grey30", "grey70"))
  barplot(sortDat$V7, border=F, ylab=expression(sigma[beta]), col=ifelse(-log10(pnorm(-abs(sortDat$V8))*2) > pThreshold1, "grey30", "grey70"))
  mt=matrix(1, dimC, modCnt)
  for(i in 1:modCnt){
    mode[i,mode[i,]==1]=0.000001
  }
  par(mai=c(0, 1, 0, 0.5))
  for(i in 1:dimC){
    barplot(t(mode)[i,], border = FALSE, axes = F, col=i,bty="n")
  }
  axis(side=1,at=5,xpd = F,font=6,labels = c("Models"), tick = F)
  plot.new()
  legend(ncol=dimC, "topright",legend = paste("Cov", seq(1, dimC)), col=c(1:dimC), pch=15, bty='n')
  dev.off()
  
  return()
}


