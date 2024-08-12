library(dendextend)
library(FNN)
library(igraph)
library(RColorBrewer)
library(pheatmap)
library(Cairo)
library(proxy)
library(data.table)
library(digest)


datain <- fread('result/VarMAT/Pixel/FiltDNAvar/mut.mat', header = TRUE)
colnames(datain) = gsub('^.*/','',gsub('\\..*$','',colnames(datain)))
cluster_file = "/home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/cluster_and_position_LM0623_Large_RNAsomatic.csv"
clusters = read.table(cluster_file, header = TRUE, sep = ',', stringsAsFactors = FALSE)
clusterkinds = unique(clusters[,2])
therownames = datain[,1]
datain = datain[,-1]

tumorpixels  = clusters[clusters[,2] %in% c(1,3),1]
normalpixels = clusters[!clusters[,2] %in% c(1,3),1]

selcted_cols = colnames(datain)[colnames(datain) %in% clusters[,1]]
datain = datain[,..selcted_cols]
datain[datain == 0] <- NA
datain[datain == 1] <- 0
datain[datain > 1] <- 1
datain[,'bk'] = 0
selcted_cols = colnames(datain)[which(apply(datain,2,FUN=function(xx){sum(!is.na(xx))}) != 0)]
datain = datain[,..selcted_cols]



hash_values <- sapply(datain, function(col) digest(col))
duplicated_columns <- which(duplicated(hash_values) | duplicated(hash_values, fromLast = TRUE))

todeletecols = c()
pixelsextend = list()
if (length(duplicated_columns) > 0) {
  for (col in unique(hash_values[duplicated_columns])) {
    thedupcols = which(hash_values == col)
    pixelsextend[[as.character(thedupcols[1])]] = colnames(datain)[thedupcols]
    todeletecols = c(todeletecols,thedupcols[-1])
  }
}
datain[, (todeletecols) := NULL]
datain = datain[apply(datain,1,FUN = function(xx){length(na.omit(xx))}) >= 2,]

#mutnumber = apply(datain != 0,1,FUN = function(xx){sum(xx,na.rm = T)})
#nornumber = apply(datain == 0,1,FUN = function(xx){sum(xx,na.rm = T)})+1
#logpossnumber = -log(mutnumber/(mutnumber+nornumber))

custom_distance <- function(x, y) {
  values = na.omit(x + y)
  values = -0.2*sum(values == 2) + sum(values == 1) - 0.1*sum(values == 0)
  return(values)
}
distance_matrix <- as.matrix(dist(t(datain), method = custom_distance))
distance_matrix <- readRDS('distance.rds')

hc = hclust(as.dist(distance_matrix),method = 'complete')
finalCluster = cutree(hc, k = 20)
table(finalCluster[normalpixels])
table(finalCluster[normalpixels])/length(normalpixels)
table(finalCluster[tumorpixels])
table(finalCluster[tumorpixels])/length(tumorpixels)


dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 8)
plot(dend, main = "Hierarchical Clustering Dendrogram with 8 Clusters", 
     ylab = "Height")


finalCluster = finalCluster[1:6726]
adds = c()
for (ii in pixelsextend) {
  select = ii[1]
  add = rep(finalCluster[select],length(ii)-1) 
  names(add) = ii[-1]
  adds = c(adds,add)
}
finalCluster = c(finalCluster,adds)

position = do.call(rbind, strsplit(names(finalCluster),'x'))
position = cbind(position,finalCluster)
position = apply(position,2,as.numeric)
outputmat = matrix(rep(NA,10000),ncol =100, nrow = 100)
position[,3] = as.factor(position[,3])
for (ii in 1:dim(position)[1]) {
  xvalue = position[ii,1]
  yvalue = position[ii,2]
  outputmat[yvalue,xvalue] = position[ii,3]
}
pheatmap(outputmat, cluster_cols=FALSE,cluster_rows=FALSE,na_col = "white")

#colnames(position) = c('x','y','cluster')
#write.table(position, file = paste0(outfile,".clu.txt"), row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE, append = FALSE)

CairoPNG(outfile,width = 1000,height = 1000)
pheatmap(outputmat, cluster_cols=FALSE,cluster_rows=FALSE,na_col = "white")
dev.off()

CairoPNG(paste0(outfile,".tree.png"),width = 1100,height = 1000)
dend <- as.dendrogram(hc)
labels_cex(dend) = NA
valuecolor = c("skyblue", "grey", "orange","purple",'cyan' ,'green', 'pink','black','darkred','chocolate','coral','drakblue')
dend %>%  set("branches_k_color", value = valuecolor[1:k], k=k) %>% plot(axes = FALSE,horiz =TRUE)
rect.dendrogram(dend, k=4, lty = 5,  lwd = 0.1, horiz = TRUE, x=1, col=rgb(0.1, 0.2, 0.4, 0.1) ) 
dev.off()







HCbyCor <- function(filein,cluin,outfile,k) {
  library(dendextend)
  library(FNN)
  library(igraph)
  library(RColorBrewer)
  library(pheatmap)
  library(Cairo)
  library(proxy)
  
  colors = unique(c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"), brewer.pal(8, "Set1")))
  datain <- read.table(filein, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  clusters = read.table(cluin, header = TRUE, sep = ',', stringsAsFactors = FALSE)
  clusters = paste0('X',clusters[,1])
  datain = datain[,colnames(datain) %in% clusters]
  datain[datain == 0] = NA
  datain = as.matrix(apply(datain,2,as.numeric))
  
  custom_similar <- function(x, y) {
    selectcols = (! is.na(x)) & x != 0 & y!= 0  & (! is.na(y))
    if (sum(selectcols) <= 2) {
      return(-2)
    }
    if (var(x[selectcols]) == 0 || var(y[selectcols]) == 0 ) {
      return(-2)
    }
    output = cor(x[selectcols],y[selectcols])
    return(output)
  }
  similar_matrix <- as.matrix(dist(t(datain), method = custom_similar))
  similar_matrix[similar_matrix == -2] = mean(similar_matrix[similar_matrix != -2])
  write.table(similar_matrix, file = paste0(outfile,".txt"), row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE, append = FALSE)
  #similar_matrix <- read.table(paste0(outfile,".txt"), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  distmat = as.dist(1-similar_matrix)
  hc = hclust(distmat,method = 'ward.D2')
  finalCluster = cutree(hc, k = k)
  position = do.call(rbind, strsplit(gsub(pattern = 'X', replacement = '', x = colnames(datain)),'x'))
  position = cbind(position,finalCluster)
  position = apply(position,2,as.numeric)
  outputmat = matrix(rep(NA,2500),ncol =50, nrow = 50)
  
  for (ii in 1:dim(position)[1]) {
    xvalue = position[ii,1]
    yvalue = position[ii,2]
    outputmat[yvalue,xvalue] = position[ii,3]
  }
  #colnames(position) = c('x','y','cluster')
  #write.table(position, file = paste0(outfile,".clu.txt"), row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE, append = FALSE)
  
  CairoPNG(outfile,width = 1000,height = 1000)
  pheatmap(outputmat, color = colors, breaks = 0.5:(length(unique(finalCluster))+0.5), cluster_cols=FALSE,cluster_rows=FALSE,na_col = "white")
  dev.off()
  
  CairoPNG(paste0(outfile,".tree.png"),width = 1100,height = 1000)
  dend <- as.dendrogram(hc)
  labels_cex(dend) = NA
  valuecolor = c("skyblue", "grey", "orange","purple",'cyan' ,'green', 'pink','black','darkred','chocolate','coral','drakblue')
  dend %>%  set("branches_k_color", value = valuecolor[1:k], k=k) %>% plot(axes = FALSE,horiz =TRUE)
  rect.dendrogram(dend, k=k, lty = 5,  lwd = 0.1, horiz = TRUE, x=1, col=rgb(0.1, 0.2, 0.4, 0.1) ) 
  dev.off()
  
}


filein = gzfile("LM0623ratio.tsv.gz")
cluin  = "cluster_and_position_LM0623.csv"
outfile = "LM0623ratio.png"
k = 4
HCbyCor(filein,cluin,outfile,k)

filein = "LM1ratio.tsv"
cluin  = "cluster_and_spatial.id_LM1.csv"
outfile = "LM1ratio.png"
k = 3
HCbyCor(filein,cluin,outfile,k)

filein = "HD837ratio.tsv"
cluin  = "cluster_and_spatial.id_HD837.csv"
outfile = "HD837ratio.png"
k = 6
HCbyCor(filein,cluin,outfile,k)








KNNbyOwnDist <- function(filein,cluin,outfile) {
  
  library(FNN)
  library(igraph)
  library(RColorBrewer)
  library(pheatmap)
  library(Cairo)
  library(proxy)
  
  colors = unique(c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"), brewer.pal(8, "Set1")))
  datain <- read.table(filein, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  clusters = read.table(cluin, header = TRUE, sep = ',', stringsAsFactors = FALSE)
  clusters = paste0('X',clusters[,1])
  datain = datain[,colnames(datain) %in% clusters]
  datain[datain == 0] = NA
  datain = as.matrix(apply(datain,2,as.numeric))
  
  custom_distance <- function(x, y) {
    selectcols = (! is.na(x)) & (! is.na(y))
    if (sum(selectcols) <= 20) {
      return(1)
    }
    output = sqrt(mean((x[selectcols]-y[selectcols])^2))
    return(output)
  }
  
  custom_distance0 <- function(x, y) {
    selectcols = (! is.na(x)) & (! is.na(y))
    if (sum(selectcols) <= 10) {
      return(1)
    }
    x <- as.factor(x[selectcols])
    y <- as.factor(y[selectcols])
    chisq_test <- tryCatch({
      chisq.test(x, y)
    }, error = function(e) NA)
    if (is.na(chisq_test)) {
      return(1)
    }
    n <- sum(chisq_test$observed)  
    k <- min(length(unique(x)), length(unique(y))) 
    cramers_v <- sqrt(chisq_test$statistic / (n * (k - 1)))
    return(cramers_v)
  }
  custom_distance <- function(x, y) {
    selectcols = (! is.na(x)) & x != 0 & y!= 0  & (! is.na(y))
    if (sum(selectcols) <= 20) {
      return(-1)
    }
    output = cor(x[selectcols],y[selectcols])
    return(output)
  }
  distance_matrix <- as.matrix(dist(t(datain), method = custom_distance))
  similar_matrix = 1 - distance_matrix
  
  k <- 30
  knn_indices <- t(apply(similar_matrix, 1, order, decreasing = FALSE))[, 2:(k+1)]
  calculate_shared_neighbors <- function(knn_indices) {
    shared_neighbors <- matrix(0, nrow = nrow(knn_indices), ncol = nrow(knn_indices))
    for (i in 1:nrow(knn_indices)) {
      for (j in 1:nrow(knn_indices)) {
        shared_neighbors[i, j] <- sum(knn_indices[i, ] %in% knn_indices[j, ])
      }
    }
    return(shared_neighbors)
  }
  shared_neighbors <- calculate_shared_neighbors(knn_indices)
  threshold <-15 #change
  adj_matrix <- shared_neighbors >= threshold
  snn_graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  clusters <- cluster_louvain(snn_graph)
  finalCluster <- factor(membership(clusters))
  write.table(distance_matrix, file = "dist_.txt", row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE, append = FALSE)
  
  position = do.call(rbind, strsplit(gsub(pattern = 'X', replacement = '', x = colnames(datain)),'x'))
  position = cbind(position,finalCluster)
  position = apply(position,2,as.numeric)
  outputmat = matrix(rep(NA,2500),ncol =50, nrow = 50)
  
  for (ii in 1:dim(position)[1]) {
    xvalue = position[ii,1]
    yvalue = position[ii,2]
    outputmat[yvalue,xvalue] = position[ii,3]
  }
  
  CairoPNG(outfile,width = 1000,height = 1000)
  pheatmap(outputmat, color = colors, breaks = -0.5:16.5, cluster_cols=FALSE,cluster_rows=FALSE,na_col = "white")
  dev.off()
  return(distance_matrix)
}


HCbyCor <- function(filein,cluin,outfile,k) {
  library(dendextend)
  library(FNN)
  library(igraph)
  library(RColorBrewer)
  library(pheatmap)
  library(Cairo)
  library(proxy)
  
  colors = unique(c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"), brewer.pal(8, "Set1")))
  datain <- read.table(filein, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  clusters = read.table(cluin, header = TRUE, sep = ',', stringsAsFactors = FALSE)
  clusters = paste0('X',clusters[,1])
  datain = datain[,colnames(datain) %in% clusters]
  datain[datain == 0] = NA
  datain = as.matrix(apply(datain,2,as.numeric))
  
  similar_matrix <- read.table(paste0(outfile,".txt"), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  distmat = as.dist(1-similar_matrix)
  hc = hclust(distmat,method = 'ward.D2')
  finalCluster = cutree(hc, k = k)
  position = do.call(rbind, strsplit(gsub(pattern = 'X', replacement = '', x = colnames(datain)),'x'))
  position = cbind(position,finalCluster)
  position = apply(position,2,as.numeric)
  outputmat = matrix(rep(NA,2500),ncol =50, nrow = 50)
  
  for (ii in 1:dim(position)[1]) {
    xvalue = position[ii,1]
    yvalue = position[ii,2]
    outputmat[yvalue,xvalue] = position[ii,3]
  }
  #colnames(position) = c('x','y','cluster')
  #write.table(position, file = paste0(outfile,".clu.txt"), row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE, append = FALSE)
  
  CairoPNG(outfile,width = 1000,height = 1000)
  pheatmap(outputmat, color = colors, breaks = 0.5:(length(unique(finalCluster))+0.5), cluster_cols=FALSE,cluster_rows=FALSE,na_col = "white")
  dev.off()
  
  CairoPNG(paste0(outfile,".tree.png"),width = 1100,height = 1000)
  dend <- as.dendrogram(hc)
  labels_cex(dend) = NA
  valuecolor = c("skyblue", "grey", "orange","purple",'cyan' ,'green', 'pink','black','darkred','chocolate','coral','drakblue')
  dend %>%  set("branches_k_color", value = valuecolor[1:k], k=k) %>% plot(axes = FALSE,horiz =TRUE)
  rect.dendrogram(dend, k=k, lty = 5,  lwd = 0.1, horiz = TRUE, x=1, col=rgb(0.1, 0.2, 0.4, 0.1) ) 
  dev.off()
  
}


filein = gzfile("LM0623ratio.tsv.gz")
cluin  = "cluster_and_position_LM0623.csv"
outfile = "LM0623ratio.png"
k = 4
HCbyCor(filein,cluin,outfile,k)

filein = "LM1ratio.tsv"
cluin  = "cluster_and_spatial.id_LM1.csv"
outfile = "LM1ratio.png"
k = 3
HCbyCor(filein,cluin,outfile,k)

filein = "HD837ratio.tsv"
cluin  = "cluster_and_spatial.id_HD837.csv"
outfile = "HD837ratio.png"
k = 6
HCbyCor(filein,cluin,outfile,k)




library(dplyr)
library(data.table)
library(parallel)
library(stringr)

datain <- fread('result/VarMAT/Pixel/FiltDNAvar/mutratio.mat', header = TRUE)
colnames(datain) = gsub('^.*/','',gsub('\\..*$','',colnames(datain)))
cluster_file = "/home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/cluster_and_position_LM0623_Large_RNAsomatic.csv"
clusters = read.table(cluster_file, header = TRUE, sep = ',', stringsAsFactors = FALSE)
clusterkinds = unique(clusters[,2])

selcted_cols = colnames(datain)[colnames(datain) %in% clusters[,1]]
datain = datain[,..selcted_cols]
selcted_cols = colnames(datain)[which(apply(datain,2,FUN=function(xx){sum(!is.na(xx))}) != 0)]
datain = datain[,..selcted_cols]
datain[is.na(datain)] = 0

sum(apply(datain,2,FUN=function(xx){sum(!is.na(xx))}) <= 5)




library(data.table)
library(digest) # 用于生成列的哈希值

# 创建一个示例data.table
dt <- data.table(
  A = c(1, 2, 3, 4),
  B = c(1, 2, 3, 4),
  C = c(1, 1, 1, 1),
  D = c(4, 3, 2, 1)
)

# 生成每一列的哈希值
hash_values <- sapply(dt, function(col) digest(col))

# 找出哈希值相同的列
duplicated_columns <- which(duplicated(hash_values) | duplicated(hash_values, fromLast = TRUE))

# 打印结果
if (length(duplicated_columns) > 0) {
  cat("完全相同的列：\n")
  for (col in unique(hash_values[duplicated_columns])) {
    cat(names(hash_values[hash_values == col]), "\n")
  }
} else {
  cat("没有完全相同的列。\n")
}