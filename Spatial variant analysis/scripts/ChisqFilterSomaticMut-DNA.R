library(dplyr)
library(data.table)
library(parallel)
library(stringr)


data = fread('VarMAT_DNAvar/mutratiostr.mat')
colnames(data) = gsub('\\..*$','',gsub('^.*/','',colnames(data)))

data[, c("dedup_nor_1", "dedup_nor_2") := tstrsplit(dedup_nor, ",", fixed = TRUE)]
data[, c("dedup_tmr_1", "dedup_tmr_2") := tstrsplit(dedup_tmr, ",", fixed = TRUE)]
data[, c("normal_1", "normal_2") := tstrsplit(normal, ",", fixed = TRUE)]
data[, c("tumor_1", "tumor_2") := tstrsplit(tumor, ",", fixed = TRUE)]

data = data[,c(1,6:13)]
data[is.na(data)] <- 0
data[, (2:ncol(data)) := lapply(.SD, as.numeric), .SDcols = 2:ncol(data)]

data2 = data[data$dedup_nor_1 * data$dedup_tmr_2 < 0.1* data$dedup_tmr_1 * data$dedup_nor_2 ]
data3 = data2[,2:5]
data3 = unique(data3)
fisherfun <- function(xx){
  x1 = xx[1]
  n1 = xx[2]
  x2 = xx[3]
  n2 = xx[4]
  data <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2)
  result <- fisher.test(data,alternative = 'less')
  result$p.value
}

pvs = apply(data3, 1, FUN = function(xx){
  fisherfun(xx)
})

#data.table::set(data3, j = "pvalue", value = pvs)
#data4 <- merge(data2, data3[, .(dedup_nor_1, dedup_nor_2, dedup_tmr_1, dedup_tmr_2, pvalue)], by = c("dedup_nor_1", "dedup_nor_2", "dedup_tmr_1", "dedup_tmr_2"), all.x = TRUE)
#padjustvs = p.adjust(data4$pvalue,method = 'BH')

data3 = data3[pvs <= 0.001,]

aa <- data2[data3, on = .(dedup_nor_1, dedup_nor_2, dedup_tmr_1, dedup_tmr_2)]
col_sums <- unlist(aa[, lapply(.SD, sum), .SDcols = 2:ncol(aa)])
(col_sums[7]/col_sums[8])/(col_sums[5]/col_sums[6])


write.table(aa[,1], file = "VarMAT_ChisqFiltDNAvar/PASS.site.pos.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)


