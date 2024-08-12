library(infercnv)

matin = "result/CNV/inferCNV/result/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"
cluster_file = "scRNAmir/data/FFPE/Clusters/Cluster_and_spatial.id_LM0623_Large.csv"


datain = read.table(matin, header = TRUE, sep = ' ', stringsAsFactors = FALSE)
colnames(datain) = gsub('^X','',colnames(datain))
clusters = read.table(cluster_file, header = TRUE, sep = ',', stringsAsFactors = FALSE)
transloc = read.table('cache/inferCNV/transloc.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
transloc = transloc[match(rownames(datain),transloc[,1]),]


clu7 = clusters[clusters[,2] == 7,1]
matclu7 = datain[,clu7] > 1
ratioclu7 = apply(matclu7,1,FUN=function(xx){sum(xx)/length(xx)})


write.table(x = transloc[ratioclu7 >= 0.5,], file = "result/CNV/inferCNV/result/select.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)