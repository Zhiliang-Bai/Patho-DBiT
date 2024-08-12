library(ggplot2)
library(Cairo)

clusters = read.table("/home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/Cluster_and_spatial.id_LM0623_Large.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
clusterkinds = unique(clusters[,2])

mutin = fread('result/VarMAT/Pixel/FiltDNAvar/mutratiostr.mat')
mutname  =  mutin[,1:2]
mutin  =  mutin[,-(1:2)]
colnames(mutin) = gsub('\\..*$', '', gsub("^.*/",'',colnames(mutin)))
selcted_cols = colnames(mutin)[colnames(mutin) %in% clusters[,1]]
mutin = mutin[,..selcted_cols]

countmut <- mutin[, lapply(.SD, function(xx){ as.numeric(gsub(",.*$", "", xx))}), .SDcols = names(mutin)]
countread <- mutin[, lapply(.SD, function(xx){ as.numeric(gsub("^.*,", "", xx))}), .SDcols = names(mutin)]

mutnumber <- apply(countmut, 2, function(x) sum(x,na.rm = TRUE))
readnumber <- apply(countread, 2, function(x) sum(x,na.rm = TRUE))
mutnumber <- mutnumber[readnumber != 0]
readnumber <- readnumber[readnumber != 0]
ratio = mutnumber/readnumber

xyframe = apply(do.call(rbind,strsplit(names(ratio),'x')),2,as.numeric)
output = data.frame(xyframe,ratio)
colnames(output) = c('x','y','value')
output[,2] = 100-output[,2]

p <- ggplot(output, aes(x = x, y = y, color = value)) + geom_point(size = 2)
print(p)
dev.off()


SeeCluMut <- clusters
colnames(SeeCluMut) = c('pixel','cluster')
countmat = data.frame('pixel' = names(ratio),ratio)
SeeCluMut = merge(SeeCluMut,countmat,by = 'pixel')
SeeCluMut = na.omit(SeeCluMut)
cluster_means <- SeeCluMut %>% group_by(cluster) %>% summarise(mean_ratio = median(ratio))

mutdepthframe <- clusters
colnames(mutdepthframe) = c('pixel','cluster')
countmat = data.frame('pixel' = names(ratio),mutnumber,readnumber)
mutdepthframe = merge(mutdepthframe,countmat,by = 'pixel')
mutdepthframe = na.omit(mutdepthframe)
mutdepthframe <- mutdepthframe %>% group_by(cluster) %>% summarise(mutnumber = sum(mutnumber),readnumber = sum(readnumber),)


fordraw = data.frame(mutdepthframe[,1],t(apply(mutdepthframe,1,FUN=function(xx){
  binomout <- binom.test(xx[2],xx[3])
  return(c(binomout$estimate,binomout$conf.int[1],binomout$conf.int[2]))
})))
colnames(fordraw) = c('Cluster','Estimate', 'Q2.5', 'Q97.5')
fordraw = arrange(fordraw,Estimate)
fordraw$Cluster = factor(fordraw$Cluster, levels = fordraw$Cluster)
CairoPNG(filename = 'result/VarMAT/Pixel/FiltDNAvar/cluster_CI.png',width = 500,height = 600)
p = ggplot() + geom_pointrange(data=fordraw, mapping=aes(x=Cluster, y=Estimate, ymin=Q2.5, ymax=Q97.5, color = Estimate), fatten = 2, size=0.5)  + scale_color_gradient(low=rgb(192,188,221,max = 255), high=rgb(111,95,164,max = 255)) + xlab('') + ylab('') + coord_flip() + theme(panel.grid.major.y = element_blank(),panel.grid.major.x = element_line(color = 'black', size = 0.03), panel.border = element_rect(fill = NA, color = 'black', size = 0.03), panel.background = element_blank()) + labs(color='Mutation rate') + xlab('Cluster')
print(p)
dev.off()
write.table(fordraw, file = "result/VarMAT/Pixel/FiltDNAvar/cluster_CI.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE, append = FALSE)