homedir=getwd()
source("oathPlot.R")

#analysis for LowMgS

if(!dir.exists("LowMgS"))
{
  dir.create("LowMgS")
}
setwd("LowMgS")

cmdLowMgSnns=paste("java -jar ../gear.jar  nss --bfile ../ArabData/Arab295Line --pheno ../ArabData/Arab295Line.phe --mpheno 11 --covar ../ArabData/Arab295Line.eigenvec --covar-number 1 2 3 4 5 --out Arab295Line_LowMgS")
system(cmdLowMgSnns)
cmdLowMgSoathX=paste("java -jar ../gear.jar  oathX --nss-batch Arab295Line_LowMgS.list.nss --cm Arab295Line_LowMgS.m.nss --n 295 --out Arab295Line_LowMgS")
system(cmdLowMgSoathX)
system("grep 3_8965883 *LowMgS*.oath > 3_8965883LowMgS.txt")
oathPlot("3_8965883LowMgS.txt", 156744, 0.05, "3_8965883LowMgS.pdf")


#analysis for HighMgRGT
setwd(homedir)
if(!dir.exists("HighMgRGT"))
{
  dir.create("HighMgRGT")
}
setwd("HighMgRGT")
cmdHighMgRGTnns=paste("java -jar ../gear.jar  nss --bfile ../ArabData/Arab295Line --pheno ../ArabData/Arab295Line.phe --mpheno 29 --covar ../ArabData/Arab295Line.eigenvec --covar-number 1 2 3 4 5 --out Arab295Line_HighMgRGT")
system(cmdHighMgRGTnns)
cmdHighMgRGToathX=paste("java -jar ../gear.jar  oathX --nss-batch Arab295Line_HighMgRGT.list.nss --cm Arab295Line_HighMgRGT.m.nss --n 295 --out Arab295Line_HighMgRGT")
system(cmdHighMgRGToathX)
system("grep 4_6353940 *HighMgRGT*.oath > 4_6353940HighMgRGT.txt")
oathPlot("4_6353940HighMgRGT.txt", 156744, 0.05, "4_6353940HighMgRGT.pdf")


#analysis for HighMgK
setwd(homedir)
if(!dir.exists("HighMGK"))
{
  dir.create("HighMgK")
}
setwd("HighMgK")
cmdHighMgKnns=paste("java -jar ../gear.jar  nss --bfile ../ArabData/Arab295Line --pheno ../ArabData/Arab295Line.phe --mpheno 36 --covar ../ArabData/Arab295Line.eigenvec --covar-number 1 2 3 4 5 --out Arab295Line_HighMgK")
system(cmdHighMgKnns)
cmdHighMgKoathX=paste("java -jar ../gear.jar  oathX --nss-batch Arab295Line_HighMgK.list.nss --cm Arab295Line_HighMgK.m.nss --n 295 --out Arab295Line_HighMgK")
system(cmdHighMgKoathX)
system("grep 5_20010406 *HighMgK*.oath > 5_20010406HighMgK.txt")
oathPlot("5_20010406HighMgK.txt", 156744, 0.05, "5_20010406HighMgK.pdf")
