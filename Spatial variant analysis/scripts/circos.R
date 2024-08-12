library(Cairo)
library(circlize)


germline4circos <- function(filein){
  germlinebed = read.table(gzfile(filein), sep = "\t")
  deeps = mapply(function(aa,bb){
    selectpos = aa %in% c('DP','DPF')
    values = bb
    values = values[selectpos]
    return(sum(as.numeric(values)))
  }, aa = strsplit(germlinebed$V9,':'),bb = strsplit(germlinebed$V10,':'))
  germlinebed = germlinebed[germlinebed[,1] != 'M' & grepl("PASS",germlinebed$V7) & deeps >= 20,]
  
  breaks = "Sturges"
  include.lowest = TRUE
  right = TRUE
  yy = NULL
  for (i in unique(germlinebed$V1)) {
    nx = germlinebed$V2[germlinebed$V1 == i]
    h = hist(nx, plot = FALSE, breaks = breaks, include.lowest = include.lowest, 
             right = right)
    yy = c(yy, 0, h$counts)
  }
  summary(yy)
  breaks = "Sturges"
  include.lowest = TRUE
  right = TRUE
  yy = NULL
  for (i in unique(germlinebed$V1)) {
    nx = germlinebed$V2[germlinebed$V1 == i]
    h = hist(nx, plot = FALSE, breaks = breaks, include.lowest = include.lowest, 
             right = right)
    yy = c(yy, 0, h$counts)
  }
  print(summary(yy))
  
  return(germlinebed)
}
RNAtumorbed = germline4circos('result/Strelka/RNA/tumor_germline/results/variants/variants.vcf.gz')
RNAnormalbed = germline4circos("result/Strelka/RNA/normal_germline//results/variants/variants.vcf.gz")
DNAnormalbed = germline4circos("result/Strelka/DNA/normal_germline/results/variants/variants.vcf.gz")
DNAtumorbed = germline4circos("result/Strelka/DNA/tumor_germline/results/variants/variants.vcf.gz")


cytoband.df = read.table("/home/dz287/gibbs/scRNAmir/data/Genomes/T2T_CHM13/chm13v2.0_cytobands_allchrs.bed", colClasses = c("character", "numeric","numeric", "character", "character"), sep = "\t")
cytoband.df[,1] = gsub('chr','',cytoband.df[,1])

CairoSVG("result/circos.svg",pointsize = 8,height = 7,width = 7)
circos.clear()
par(cex = 1.25)
circos.par(gap.degree = c(rep(1,13),12,rep(1,10)),track.height = 0.11)
circos.initializeWithIdeogram(cytoband.df,plotType = c("ideogram", "labels"))
circos.trackHist(DNAtumorbed$V1, x = DNAtumorbed$V2, ylim = c(0,30000),col = '#72A2C0', border = FALSE)
circos.yaxis(side = "left", sector.index = '15',labels.cex = 0.7,at = c(0,15000,30000))
circos.trackHist(DNAnormalbed$V1, x = DNAnormalbed$V2, ylim = c(0,30000),col = 'red', border = FALSE)
circos.yaxis(side = "left", sector.index = '15',labels.cex = 0.7,at = c(0,15000,30000))
circos.trackHist(RNAtumorbed$V1, x = RNAtumorbed$V2, ylim = c(0,560),col = '#FDD262', border = FALSE)
circos.yaxis(side = "left", sector.index = '15',labels.cex = 0.7,at = c(0,280,560))
circos.trackHist(RNAnormalbed$V1, x = RNAnormalbed$V2, ylim = c(0,560),col = '#595959', border = FALSE)
circos.yaxis(side = "left", sector.index = '15',labels.cex = 0.7,at = c(0,280,560))
circos.clear()
dev.off()


CairoSVG("result/circos2.svg",pointsize = 8,height = 7,width = 7)
circos.clear()
par(cex = 1.25)
circos.par(gap.degree = c(rep(1,13),15,rep(1,10)),track.height = 0.30)
circos.initializeWithIdeogram(cytoband.df,plotType = c("ideogram", "labels"))
circos.trackHist(RNAtumorbed$V1, x = RNAtumorbed$V2, ylim = c(0,560),col = '#FDD262', border = FALSE)
circos.yaxis(side = "left", sector.index = '15',labels.cex = 0.7,at = c(0,280,560))
circos.trackHist(RNAnormalbed$V1, x = RNAnormalbed$V2, ylim = c(0,560),col = '#595959', border = FALSE)
circos.yaxis(side = "left", sector.index = '15',labels.cex = 0.7,at = c(0,280,560))
circos.clear()
dev.off()
