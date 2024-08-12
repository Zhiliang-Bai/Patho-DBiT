library(dplyr)
library(data.table)
library(parallel)
library(stringr)


data = fread('result/VarMAT/Sample/DNAvar/mutratiostr.mat')
colnames(data) = gsub('\\..*$','',gsub('^.*/','',colnames(data)))

data[, c("dedup_nor_1", "dedup_nor_2") := tstrsplit(dedup_nor, ",", fixed = TRUE)]
data[, c("dedup_tmr_1", "dedup_tmr_2") := tstrsplit(dedup_tmr, ",", fixed = TRUE)]
data[, c("normal_1", "normal_2") := tstrsplit(normal, ",", fixed = TRUE)]
data[, c("tumor_1", "tumor_2") := tstrsplit(tumor, ",", fixed = TRUE)]

data = data[,c(1,6:13)]
data[is.na(data)] <- 0
data[, (2:ncol(data)) := lapply(.SD, as.numeric), .SDcols = 2:ncol(data)]


aa = data[dedup_tmr_1 >= 10 & dedup_tmr_1 * dedup_nor_2  > 2* dedup_nor_1 * dedup_tmr_2 ]
col_sums <- unlist(aa[, lapply(.SD, sum), .SDcols = 2:ncol(aa)])
(col_sums[7]/col_sums[8])/(col_sums[5]/col_sums[6])
#1.568279

write.table(aa[,1], file = "result/VarMAT/Pixel/FiltDNAvar/PASS.site.pos.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)



aa = data[dedup_nor_1 <= 1 & dedup_tmr_1  > 0.45* dedup_tmr_2 ]
col_sums <- unlist(aa[, lapply(.SD, sum), .SDcols = 2:ncol(aa)])
(col_sums[7]/col_sums[8])/(col_sums[5]/col_sums[6])