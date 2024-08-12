library(Cairo)
library(circlize)
Normalin = read.table("result/Figures/VCFs4Figs/filteredSNVcluster_Normal.vcf", sep = "\t")
CLu1in = read.table("result/Figures/VCFs4Figs/filteredSNVclusters1.vcf", sep = "\t")
CLu3in = read.table("result/Figures/VCFs4Figs/filteredSNVclusters3.vcf", sep = "\t")

Normalin = Normalin[Normalin$V1 != 'M',]
CLu1in = CLu1in[CLu1in$V1 != 'M',]
CLu3in = CLu3in[CLu3in$V1 != 'M',]

getscale <- function(vcfin) {
  breaks = "Sturges"
  include.lowest = TRUE
  right = TRUE
  yy = NULL
  for (i in unique(vcfin$V1)) {
    nx = vcfin$V2[vcfin$V1 == i]
    h = hist(nx, plot = FALSE, breaks = breaks, include.lowest = include.lowest, 
             right = right)
    yy = c(yy, 0, h$counts)
  }
  return(summary(yy))
}

getscale(Normalin)
getscale(CLu1in)
getscale(CLu3in)

cytoband.df = read.table("karyotype.human.hg38.bed", colClasses = c("character", "numeric","numeric", "character", "character"), sep = "\t")
cytoband.df[,1] = gsub('chr','',cytoband.df[,1])
CairoPDF("result/Figures/circos/circos_germline_1_3.pdf",pointsize = 8,height = 7,width = 7,family = "Arial")
circos.clear()
par(cex = 1.25,family = "Arial")
circos.par(gap.degree = c(rep(1,13),10,rep(1,10)),track.height = 0.12)#,track.margin = c(0.02,0.02)
circos.initializeWithIdeogram(cytoband.df,plotType = c("ideogram", "labels"))
circos.trackHist(Normalin$V1, x = Normalin$V2, ylim = c(0,750),col = '#424242', border = FALSE)
circos.yaxis(side = "left", sector.index = '15',labels.cex = 0.5,at = c(0,375,750))
circos.trackHist(CLu1in$V1, x = CLu1in$V2, ylim = c(0,1200),col = '#EB545C', border = FALSE)
circos.yaxis(side = "left", sector.index = '15',labels.cex = 0.5,at = c(0,600,1200))
circos.trackHist(CLu3in$V1, x = CLu3in$V2, ylim = c(0,1000),col = '#1D65A6', border = FALSE)
circos.yaxis(side = "left", sector.index = '15',labels.cex = 0.5,at = c(0,500,1000))
circos.clear()
dev.off()

