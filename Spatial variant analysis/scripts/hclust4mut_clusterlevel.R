library(dendextend)
library(FNN)
library(igraph)
library(RColorBrewer)
library(pheatmap)
library(Cairo)
library(proxy)
library(data.table)
library(digest)


# An example of Cluster_and_spatial.id.csv
#"","SCT_snn_res.1.2"
#"100x1","8"
#"100x2","4"
#"100x3","15"
#"100x4","15"
#"100x5","4"
#"100x6","3"
#"100x7","18"
#"100x8","8"
#"100x9","3"

cluster_file = "Cluster_and_spatial.id.csv"
mutinfile = 'VarMAT_ChisqFiltDNAvar/mutratiostr.mat'



datain <- fread(mutinfile, header = TRUE)
locationin = datain$location
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

locationin = locationin[apply(mutcount_clu,1,sum) != 0]
allcount_clu = allcount_clu[apply(mutcount_clu,1,sum) != 0,]
mutcount_clu = mutcount_clu[apply(mutcount_clu,1,sum) != 0,]

threshold = 3
ratiocount_clu = mutcount_clu >= threshold
ratiocount_clu[allcount_clu < threshold] = NA
ratiocount_clu[mutcount_clu < 0.25*allcount_clu] = 0
ratiocount_clu = cbind(0,ratiocount_clu)
colnames(ratiocount_clu)[1] = 'Reference'

custom_distance <- function(x, y) {
  values = na.omit(x + y)
  values = sum(values == 1)/length(values)
  return(values)
}
distance_matrix <- as.matrix(dist(t(ratiocount_clu), method = custom_distance))
hc = hclust(as.dist(distance_matrix),method = 'average')

phylo_tree <- as.phylo(hc)
rooted_tree <- root(phylo_tree, outgroup = "Reference", resolve.root = TRUE)
plot(rooted_tree, main = "Hierarchical Clustering Dendrogram with 'bk")
ggtree(rooted_tree) + geom_tiplab(aes(label = label), size = 3) #+ geom_tree(aes(color = label))


#we can merge clusters which are too similar to reference to reference in the level of leaf
deletenums = which(colnames(distance_matrix) != 'Reference' & distance_matrix[,1] == 0)
distance_matrix2 = distance_matrix[-deletenums,-deletenums]
hc = hclust(as.dist(distance_matrix2),method = 'average')
phylo_tree <- as.phylo(hc)
rooted_tree <- root(phylo_tree, outgroup = "Reference", resolve.root = TRUE)
plot(rooted_tree)

ggtree(rooted_tree) + geom_tiplab(aes(label = label), size = 3) #+ geom_tree(aes(color = label))


# we can look at which mutations are important mutation for each cluster
bb = ratiocount_clu
bb[is.na(bb)] = 0
aa = ratiocount_clu[,! colnames(ratiocount_clu) %in% c(3,6,7,14,18)]
aa[is.na(aa)] = 0
tosave = which(apply(aa,1,sum) > 0)

locationin[which(bb[, colnames(bb) == 6] == 1 & tosave)]
locationin[which(bb[, colnames(bb) == 7] == 1 & tosave)]
locationin[which(bb[, colnames(bb) == 14] == 1 & tosave)]
locationin[which(bb[, colnames(bb) == 18] == 1 & tosave)]








