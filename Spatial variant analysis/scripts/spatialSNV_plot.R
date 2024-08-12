library(patchwork)
library(ggplot2)

cluster_file = "/home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/cluster_and_position_LM0623_Large_RNAsomatic.csv"
mutfile = 'result/VarMAT/Pixel/FiltDNAvar/mut.mat'
outputfolder = 'result/Figures/Spatial/FiltDNAvar/'
tumorcluster = '1:3'
mutfile = 'result/VarMAT/Pixel/FiltRNAvar/mut.mat'
outputfolder = 'result/Figures/Spatial/FiltRNAvar/'
mutfile = 'result/VarMAT/Pixel/FiltRNAvar_SNVcluster/mut.mat'
outputfolder = 'result/Figures/Spatial/FiltRNAvar_SNVcluster/'
mutfile = 'result/VarMAT/Pixel/StrictFiltDNAvar/mut.mat'
outputfolder = 'result/Figures/Spatial/StrictFiltDNAvar/'
mutfile = 'result/VarMAT/Pixel/ChisqFiltDNAvar/mut.mat'
outputfolder = 'result/Figures/Spatial/ChisqFiltDNAvar/'


arg = commandArgs(T)
if (length(arg) == 4) {
  cluster_file = arg[1]
  mutfile = arg[2]
  outputfolder = arg[3]
  tumorcluster = arg[4]
}
tumorcluster = as.numeric(strsplit(tumorcluster,':')[[1]])

clusters = read.table(cluster_file, header = TRUE, sep = ',', stringsAsFactors = FALSE)
clusterkinds = unique(clusters[,2])


mutin = fread(mutfile)
mutname  =  mutin[,1]
mutin  =  mutin[,-1]
colnames(mutin) = gsub('.*/','',gsub('\\..*$','',colnames(mutin)))
selcted_cols = colnames(mutin)[colnames(mutin) %in% clusters[,1]]
mutin = mutin[,..selcted_cols]
mutin[is.na(mutin)] = 0

count1 <- apply(mutin, 1, function(x) sum(x == 1))
count2 <- apply(mutin, 1, function(x) sum(x == 2))
count3 <- apply(mutin, 1, function(x) sum(x == 3))

selcted_cols_tumor = colnames(mutin)[colnames(mutin) %in% clusters[clusters[,2] %in% tumorcluster,1]]
selcted_cols_normal = colnames(mutin)[colnames(mutin) %in% clusters[!clusters[,2] %in% tumorcluster,1]]
mut_tumor = mutin[,..selcted_cols_tumor]
mut_normal = mutin[,..selcted_cols_normal]




xyframe_tumor = apply(do.call(rbind,strsplit(colnames(mut_tumor),'x')),2,as.numeric)
xyframe_normal = apply(do.call(rbind,strsplit(colnames(mut_normal),'x')),2,as.numeric)


significant_is = which(count1 + count2 + count3 >= 5)
for (i in significant_is) {
  thevar = mutname[i]
  output_tumor = data.frame(xyframe_tumor,unlist(mut_tumor[i,]))
  output_normal = data.frame(xyframe_normal,unlist(mut_normal[i,]))
  colnames(output_normal) = c('x','y','value')
  colnames(output_tumor) = c('x','y','value')
  output = rbind(output_tumor,output_normal)
  output[,2] = 100-output[,2]
  output_tumor[,2] = 100-output_tumor[,2]
  output_normal[,2] = 100-output_normal[,2]
  Cairo::CairoPNG(filename = paste0(outputfolder,'/',thevar,'.png'),width = 2500,height = 2500)
  p <- ggplot(output, aes(x = x, y = y, color = as.factor(value))) + geom_point(size = 3.5) + scale_color_manual(values = c('0' = "lightyellow", '1' = "blue", '2' = "green", "3" = "red")) + theme_minimal() + labs(color = 'Value', title = thevar)
  p1 <- ggplot(output_tumor, aes(x = x, y = y, color = as.factor(value))) + geom_point(size = 3.5) + scale_color_manual(values = c('0' = "lightyellow", '1' = "blue", '2' = "green", "3" = "red")) + theme_minimal() + labs(color = 'Value', title = 'Tumor')
  p2 <- ggplot(output_normal, aes(x = x, y = y, color = as.factor(value))) + geom_point(size = 3.5) + scale_color_manual(values = c('0' = "lightyellow", '1' = "blue", '2' = "green", "3" = "red")) + theme_minimal() + labs(color = 'Value', title = 'Normal')
  p = (p + plot_spacer()) / (p1 + p2)
  print(p)
  dev.off()
}



