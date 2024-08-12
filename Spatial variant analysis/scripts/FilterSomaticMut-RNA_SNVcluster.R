library(dplyr)
library(data.table)
library(parallel)
library(stringr)


data = fread('VarMAT_RNAvar_SNVcluster/mutratiostr.mat')
colnames(data) = gsub('\\..*$','',gsub('^.*/','',colnames(data)))

data[, c("dedup_nor_1", "dedup_nor_2") := tstrsplit(dedup_nor, ",", fixed = TRUE)]
data[, c("dedup_tmr_1", "dedup_tmr_2") := tstrsplit(dedup_tmr, ",", fixed = TRUE)]
data[, c("normal_1", "normal_2") := tstrsplit(normal, ",", fixed = TRUE)]
data[, c("tumor_1", "tumor_2") := tstrsplit(tumor, ",", fixed = TRUE)]

data = data[,c(1,6:13)]
data[is.na(data)] <- 0
data[, (2:ncol(data)) := lapply(.SD, as.numeric), .SDcols = 2:ncol(data)]


#aa[, 'ratio' := (tumor_1 * normal_2)  / (tumor_2 * normal_1+1)]
aa = data[tumor_1 >= 0.5* tumor_2 & tumor_1  > 20 & tumor_1 * normal_2  > 1.25* tumor_2 * normal_1 ]
sum(aa[,4]/aa[,5] > aa[,2]/aa[,3],na.rm = TRUE)/(dim(aa)[1]);print(dim(aa)[1])

aa = data[tumor_1 >= 0.25* tumor_2 & tumor_1  > 10 & tumor_1 * normal_2  >  tumor_2 * normal_1 ]
sum(aa[,4]/aa[,5] > aa[,2]/aa[,3],na.rm = TRUE)/(dim(aa)[1]);print(dim(aa)[1])

write.table(aa[,1], file = "VarMAT_FiltRNAvar_SNVcluster/PASS.site.pos.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)
