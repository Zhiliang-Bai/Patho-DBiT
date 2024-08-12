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

data = data[data$dedup_nor_2 > 0 & data$dedup_tmr_2 > 0 ]
data = data[data$normal_2 + data$tumor_2 != 0]
#p normal <= 0.5 or = 0
#p normal = 0.5 p tumor = 1

aa = data[dedup_tmr_1  > 0.45* dedup_tmr_2 & dedup_nor_1 <= 0.1*dedup_nor_2  ]
p_normal_no = apply(aa,1,FUN=function(xx){
  result <- binom.test(as.numeric(xx[2]), as.numeric(xx[3]),p=0.5,alternative =  'less')
  return(result$p.value)
})
aa1 = aa[(dedup_nor_1 == 0 | p_normal_no <= 0.0001)]



fisherfun <- function(x1,n1,x2,n2){
  data <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2)
  result <- fisher.test(data,alternative = 'greater')
  result$p.value
}
aa = data[dedup_tmr_1 * dedup_nor_2  > 1.5* dedup_nor_1 * dedup_tmr_2 & dedup_tmr_1  > 0.5* dedup_tmr_2 ]
p_normal_half = apply(aa,1,FUN=function(xx){
  resultp <- fisherfun(as.numeric(xx[4]), as.numeric(xx[5]), as.numeric(xx[2]),  as.numeric(xx[3]))
  return(resultp)
})
p_normal_no = apply(aa,1,FUN=function(xx){
  result <- binom.test(as.numeric(xx[2]), as.numeric(xx[3]),p=0.5)
  return(result$p.value)
})
aa = aa[p_normal_half <= 0.0001 & p_normal_no >= 0.95]
aa = rbind(aa1,aa)
aa = aa[!duplicated(aa$location)]

write.table(aa[,1], file = "result/VarMAT/Pixel/StrictFiltDNAvar/PASS.site.pos.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)



aa0 = aa[aa$normal_2 != 0 & aa$tumor_2 != 0]
col_sums <- unlist(aa0[, lapply(.SD, sum), .SDcols = 2:ncol(aa0)])
(col_sums[7]/col_sums[8])/(col_sums[5]/col_sums[6])