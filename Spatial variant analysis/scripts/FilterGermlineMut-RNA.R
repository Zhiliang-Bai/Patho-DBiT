library(dplyr)
library(data.table)
library(parallel)
library(stringr)

datain <- fread('VarMAT_RNAgermline/mutratiostr.mat', header = TRUE)
mutcount = apply(datain[,-1],2,FUN=function(xx){as.numeric(gsub(',.*$','',xx))})
mutcount[is.na(mutcount)] = 0
allcount = apply(datain[,-1],2,FUN=function(xx){as.numeric(gsub('^.*,','',xx))})
allcount[is.na(allcount)] = 0
norcount = allcount-mutcount

datain = data.frame(datain[,1],mutcount[,1]+mutcount[,2],mutcount[,3]+mutcount[,4])
colnames(datain) = c('location','DNA','RNA')

write.table(datain[datain[,3] >= 10,1], file = "VarMAT_FiltGermlineRNAvar/PASS.site.pos.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)

richwhether = mutcount > 0.3*allcount
datain$DNAexists = richwhether[,1] | richwhether[,2]
print(sum(datain[datain[,3] >= 10,]$DNAexists > 0 & datain[datain[,3] >= 10,2] > 10)/sum(datain[,3] >= 10))
filterdatain = datain[datain[,3] >= 10,]
filterdatain$location = gsub(':.*', '', filterdatain$location)
result <- filterdatain %>% group_by(location) %>% summarize(DNA_verified_SNP = sum(DNAexists), All_SNP_Count = n())
result = as.data.frame(result)
result$ratio = result[,2]/result[,3]
write.table(result, file = "result/Figures/RNAsomatic_DNAsource_by_chr.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE, append = FALSE)
