#oathX for naCohort
cmd1=paste0("java -jar ../gear.jar oathX --nss-batch Na.chr1.nss.list.nss --cm Na.chr1.nss.m.nss --n 3191 --out Na_md")
print(cmd1)
system(cmd1)

#oathX for sleCohort
cmd2=paste0("java -jar ../gear.jar oathX --nss-batch sle.chr1.nss.list.nss --cm sle.chr1.nss.m.nss --n 2309 --out sle_md")
print(cmd2)
system(cmd2)

#meta-analysis
for (i in 1:7)
{
  cmdMeta=paste0("java -jar ../gear.jar gmeta --meta-batch meta", i, ".txt --qt 3191 2309 --out meta", i)
  print(cmdMeta)
  system(cmdMeta)
}

par(ps=9)
mn=c("s={1}", "s={2}", "s={3}", "s={1,2}", "s={1,3}", "s={2,3}", "s={1,2,3}")
layout(matrix(1:8, 2, 4, byrow = T))
for(i in 1:7)
{
  gm=read.table(paste0("meta", i,".gmeta"), as.is=T, header=T)
  plot(main=mn[i], -log10(gm$P), xlab="Chr 1", pch=16, cex=0.5, bty='n', ylab=expression(-log[10](italic(p))))
  abline(h=-log10(0.05/nrow(gm)), lty=2)
}
