library(infercnv)

matin = "infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"
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


datain = read.table(matin, header = TRUE, sep = ' ', stringsAsFactors = FALSE)
colnames(datain) = gsub('^X','',colnames(datain))
clusters = read.table(cluster_file, header = TRUE, sep = ',', stringsAsFactors = FALSE)
transloc = read.table('cache/inferCNV/transloc.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
transloc = transloc[match(rownames(datain),transloc[,1]),]


clu7 = clusters[clusters[,2] == 7,1]
matclu7 = datain[,clu7] > 1
ratioclu7 = apply(matclu7,1,FUN=function(xx){sum(xx)/length(xx)})


write.table(x = transloc[ratioclu7 >= 0.5,], file = "select.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)
