adepth = sapply(strsplit(parsed_data[,6],","), FUN = function(xx){as.numeric(xx)[2]})
dep=10;table(datain[adepth >= dep,]$exist)/sum(adepth >= dep)

clusters = read.table("/home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/Cluster_and_spatial.id_LM0623_Large.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
clusterkinds = unique(clusters[,2])


mutin = fread('result/PixelGenotypeMat.txt')
mutname  =  mutin[,1:2]
mutin  =  mutin[,-(1:2)]
selcted_cols = colnames(mutin)[colnames(mutin) %in% clusters[,1]]
mutin = mutin[,..selcted_cols]

locationstr = apply(mutname,1,FUN=function(xx){paste0(xx,collapse ='')})

count1 <- apply(mutin, 1, function(x) sum(x == 1))
count2 <- apply(mutin, 1, function(x) sum(x == 2))
count3 <- apply(mutin, 1, function(x) sum(x == 3))
xyframe = apply(do.call(rbind,strsplit(colnames(mutin),'x')),2,as.numeric)
mutin2 = mutin[count1 >= 20 & count3 >= 20,]
for (i in 1:dim(mutin2)[1]) {
  thevar = locationstr[i]
  output = data.frame(xyframe,unlist(mutin2[i,]))
  colnames(output) = c('x','y','value')
  output[,2] = 100-output[,2]
  Cairo::CairoPNG(filename = paste0('4567/',i,'.png'),width = 1200,height = 1200)
  p <- ggplot(output, aes(x = x, y = y, color = as.factor(value))) + geom_point(size = 3.5) + scale_color_manual(values = c('0' = "lightyellow", '1' = "blue", '2' = "green", "3" = "red")) + theme_minimal() + labs(color = 'Value')
  print(p)
  dev.off()
}



