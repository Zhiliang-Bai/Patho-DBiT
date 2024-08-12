library(ggplot2)
library(ggsci)

chromlens <- function(chromlengthfile,h=10,ratio=TRUE,segmentlen=10000){
  h1 = h*-0.5
  h2 = h*-1
  if (ratio) {
    h1=100*h1/segmentlen
    h2=100*h2/segmentlen
  }
  chromlength <- read.table(chromlengthfile, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  chromlength = chromlength[!chromlength[,1] %in% c('M'),]
  chromlength = data.frame(chromlength,cumsum(as.numeric(chromlength$V2))-chromlength$V2)
  colnames(chromlength) = c('chrName','Len','CumLen')
  rownames(chromlength) = chromlength[,1]
  midpos <- chromlength$CumLen +  chromlength$Len/ 2
  chromlength$midpoint = midpos
  chromlength$height = c(rep(c(h1,h2),12))
  chromlength4vline = rbind(chromlength,c('end',0,sum(chromlength[,2]),0,2.5))
  chromlength4vline[,3] = as.numeric(chromlength4vline[,3])
  chromlength2 = chromlength
  chromlength2[,1] = gsub('chr','',chromlength2[,1])
  chromlenslist = list()
  chromlenslist[[1]] = chromlength
  chromlenslist[[2]] = chromlength4vline
  chromlenslist[[3]] = chromlength2
  return(chromlenslist)
}
getdrawmat <- function(vcfinfile,chromlength,segmentlen,ratio = TRUE){
  vcfin <- read.table(vcfinfile, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  vcfin = vcfin[!vcfin[,1] %in% c('chrM'),]
  allsites = list()
  allcounts = list()
  for (ii in 1:dim(chromlength)[1]) {
    chrname = chromlength[ii,1]
    len = floor(chromlength[ii,2]/segmentlen)
    sites = 0:(len-1)*segmentlen+5000
    sites = c(sites,min(len*segmentlen+5000,chromlength[ii,2]))
    allsites[[chrname]] = sites
    allcounts[[chrname]] = rep(0,length(sites))
  }
  
  positions4all = floor(vcfin[,2]/segmentlen)
  for (ii in 1:dim(vcfin)[1]) {
    chrname = vcfin[ii,1]
    allcounts[[chrname]][positions4all[ii]] = allcounts[[chrname]][positions4all[ii]] + 1
  }
  
  fordraw = list()
  for (ii in 1:dim(chromlength)[1]) {
    chrname = chromlength[ii,1]
    fordraw[[ii]] = cbind(chrname,allsites[[chrname]],allcounts[[chrname]])
  }
  
  fordraw = do.call(rbind,fordraw)
  fordraw = data.frame(fordraw)
  fordraw[,c(2,3)] = apply(fordraw[,c(2,3)],2,as.numeric)
  drawmat = data.frame(cumpos = chromlength$CumLen[match(fordraw[,1],chromlength$chrName)]+fordraw[,2],value = fordraw[,3])
  if (ratio) {
    drawmat$value = 100*drawmat$value/segmentlen
  }
  drawmat$chr = fordraw[,1]
  return(drawmat)
}
getdrawmat0 <- function(vcfinfile,chromlength,segmentlen,ratio = TRUE){
  vcfin <- read.table(vcfinfile, header = FALSE, sep = ':', stringsAsFactors = FALSE)
  vcfin = vcfin[!vcfin[,1] %in% c('chrM'),]
  allsites = list()
  allcounts = list()
  for (ii in 1:dim(chromlength)[1]) {
    chrname = chromlength[ii,1]
    len = floor(chromlength[ii,2]/segmentlen)
    sites = 0:(len-1)*segmentlen+5000
    sites = c(sites,min(len*segmentlen+5000,chromlength[ii,2]))
    allsites[[chrname]] = sites
    allcounts[[chrname]] = rep(0,length(sites))
  }
  
  positions4all = floor(vcfin[,2]/segmentlen)
  for (ii in 1:dim(vcfin)[1]) {
    chrname = vcfin[ii,1]
    allcounts[[chrname]][positions4all[ii]] = allcounts[[chrname]][positions4all[ii]] + 1
  }
  
  fordraw = list()
  for (ii in 1:dim(chromlength)[1]) {
    chrname = chromlength[ii,1]
    fordraw[[ii]] = cbind(chrname,allsites[[chrname]],allcounts[[chrname]])
  }
  
  fordraw = do.call(rbind,fordraw)
  fordraw = data.frame(fordraw)
  fordraw[,c(2,3)] = apply(fordraw[,c(2,3)],2,as.numeric)
  drawmat = data.frame(cumpos = chromlength$CumLen[match(fordraw[,1],chromlength$chrName)]+fordraw[,2],value = fordraw[,3])
  if (ratio) {
    drawmat$value = 100*drawmat$value/segmentlen
  }
  drawmat$chr = fordraw[,1]
  return(drawmat)
}

alllist = list()
h = 60
segmentlen=10000
chromlenslist = chromlens('scRNAmir/data/Genomes/GRCh38/StarIndex/chrNameLength.txt',h=h*0.2,segmentlen=segmentlen)
drawmat = getdrawmat0('result/VarMAT/Pixel/FiltGermlineRNAvar/PASS.site.pos.txt',chromlenslist[[1]],segmentlen)
alllist[['all']] = drawmat

for (variable in 0:8) {
  h = 60
  segmentlen=10000
  chromlenslist = chromlens('scRNAmir/data/Genomes/GRCh38/StarIndex/chrNameLength.txt',h=h*0.2,segmentlen=segmentlen)
  drawmat = getdrawmat(paste0('result/Figures/VCFs4Figs/filteredSNVclusters',variable,'.vcf'),chromlenslist[[1]],segmentlen)
  alllist[[as.character(variable)]] = drawmat
}

chrorder = c("1" ,"2" ,"3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9" ,"10" ,"11" ,"12" ,"13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19" ,"20" ,"21" ,"22" ,"X" ,"Y")

outputfreqcount = list()
for (ii in 1:20) {
  outa = sapply(alllist, function(xx){
    xx[,2] = xx[,2] >= ii/100
    out = aggregate(xx$value,list(xx$chr),sum)
    out = out[match(chrorder,out[,1]),]
    return(out[,2])
  })
  outa = cbind(ii/100,chrorder,outa)
  outputfreqcount[[ii]] = outa
}
outputfreqcount = do.call(rbind,outputfreqcount)
outputfreqcount = rbind(colnames(outputfreqcount),outputfreqcount)
outputfreqcount[1,] = c('Cutoff of SNV Frequency(%) >=','chr','all_sample','c0','c1','c2','c3','c4','c5','c6','c7','c8')
write.table(data.frame(outputfreqcount), sep = '\t',file = 'result/Figures/output_chr_freqcount.txt',col.names = F,row.names = F, quote = FALSE)