library(dendextend)
library(FNN)
library(igraph)
library(RColorBrewer)
library(pheatmap)
library(Cairo)
library(proxy)
library(data.table)
library(digest)

cluster_file = "/home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/cluster_and_position_LM0623_Large_RNAsomatic.csv"
mutinfile = 'result/VarMAT/Pixel/StrictFiltDNAvar/mutratiostr.mat'
cluster_file = "/home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/Cluster_and_spatial.id_LM0623_Large.csv"

datain <- fread(mutinfile, header = TRUE)
colnames(datain) = gsub('^.*/','',gsub('\\..*$','',colnames(datain)))
clusters = read.table(cluster_file, header = TRUE, sep = ',', stringsAsFactors = FALSE)
clusterkinds = unique(clusters[,2])
therownames = datain[,1]
datain = datain[,-1]

selcted_cols = colnames(datain)[colnames(datain) %in% clusters[,1]]
datain = datain[,..selcted_cols]


mutcount = datain[, lapply(.SD, function(x) as.numeric(sub(",.*$", "", x))), .SDcols = colnames(datain)]
allcount = datain[, lapply(.SD, function(x) as.numeric(sub("^.*,", "", x))), .SDcols = colnames(datain)]


mutcount[is.na(mutcount)] <- 0
allcount[is.na(allcount)] <- 0

#identical(colnames(mutcount),clusters[,1])

mutcountlist = list()
allcountlist = list()

for (clusteri in unique(clusters[,2])) {
  selectpixels = clusters[clusters[,2] == clusteri,1]
  mutcountlist[[as.character(clusteri)]] = apply(mutcount[,..selectpixels],1,sum)
  allcountlist[[as.character(clusteri)]] = apply(allcount[,..selectpixels],1,sum)
}

mutcount_clu = do.call(cbind,mutcountlist)
allcount_clu = do.call(cbind,allcountlist)

allcount_clu = allcount_clu[apply(mutcount_clu,1,sum) != 0,]
mutcount_clu = mutcount_clu[apply(mutcount_clu,1,sum) != 0,]

threshold = 3
ratiocount_clu = mutcount_clu >= threshold
ratiocount_clu[allcount_clu <= threshold] = NA
ratiocount_clu[mutcount_clu < 0.25*allcount_clu] = 0
ratiocount_clu = cbind(0,ratiocount_clu)
colnames(ratiocount_clu)[1] = 'bk'

custom_distance <- function(x, y) {
  values = na.omit(x + y)
  values = sum(values == 1)/length(values)
  return(values)
}
distance_matrix <- as.matrix(dist(t(ratiocount_clu), method = custom_distance))
hc = hclust(as.dist(distance_matrix),method = 'average')

phylo_tree <- as.phylo(hc)
rooted_tree <- root(phylo_tree, outgroup = "bk", resolve.root = TRUE)
plot(rooted_tree, main = "Hierarchical Clustering Dendrogram with 'bk")


















ratiocount_clu <- matrix(mapply(FUN = function(xx,yy){paste0(xx,'/',yy)}, mutcount_clu, allcount_clu), ncol = 20)
colnames(ratiocount_clu) = colnames(mutcount_clu)
#ratiocount_clu = cbind('0/3',ratiocount_clu)
#colnames(ratiocount_clu)[1] = 'bk'




custom_distance1 <- function(x, y) {
  mutvalues1 = as.numeric(strsplit(x,'/')[[1]])
  mutvalues2 = as.numeric(strsplit(y,'/')[[1]])
  x1 <- mutvalues1[1]
  n1 <- mutvalues1[2]
  x2 <- mutvalues2[1]
  n2 <- mutvalues2[2]
  data <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2)
  result <- fisher.test(data)
  if (result$p.value <= 0.5) {
    return(1)
  }else{
    return(-1)
  }
}

custom_distance <- function(x,y) {
  outv = 0
  for (i in 1:length(x)) {
    outv = outv + custom_distance1(x[i],y[i])
  }
  return(outv)
}

custom_distance <- function(x,y) {
  xin = do.call(rbind,strsplit(x,'/'))
  xin = apply(xin,2,as.numeric)
  yin = do.call(rbind,strsplit(y,'/'))
  yin = apply(yin,2,as.numeric)
  ratiox = xin[,1]/xin[,2]
  ratiox[is.na(ratiox)] = 0
  ratioy = yin[,1]/yin[,2]
  ratioy[is.na(ratioy)] = 0
  outv = 1-cor(ratiox,ratioy)
  return(outv)
}


custom_distance <- function(x,y) {
  xin = do.call(rbind,strsplit(x,'/'))
  xin = apply(xin,2,as.numeric)
  yin = do.call(rbind,strsplit(y,'/'))
  yin = apply(yin,2,as.numeric)
  ratiox = xin[,1]/xin[,2]
  select = xin[,2] != 0 & yin[,2] != 0
  ratioy = yin[,1]/yin[,2]
  outv = 1-cor(ratiox[select],ratioy[select])
  return(outv)
}
ratiocount_clu <- as.matrix(ratiocount_clu)
n <- ncol(ratiocount_clu)
dist_matrix <- matrix(0,nrow = n, ncol = n)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    thisvalue = custom_distance(ratiocount_clu[, i], ratiocount_clu[, j])
    dist_matrix[i,j] <- thisvalue
    dist_matrix[j,i] <- thisvalue
  }
}
colnames(dist_matrix) = colnames(ratiocount_clu)
rownames(dist_matrix) = colnames(ratiocount_clu)
dist_object <- as.dist(dist_matrix)
hc = hclust(as.dist(dist_object),method = 'complete')

phylo_tree <- as.phylo(hc)
rooted_tree <- root(phylo_tree, outgroup = "11", resolve.root = TRUE)
plot(rooted_tree, main = "'bk' as Root")

