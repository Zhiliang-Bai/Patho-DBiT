library(dplyr)
library(data.table)
library(parallel)
library(stringr)

arg = commandArgs(T)
VCFs_path = arg[1]
outputMat = arg[2]
mutationlist = arg[3]
num_cores <- arg[4]


mutin <- read.table(mutationlist, header = FALSE, sep = ':', stringsAsFactors = FALSE)
mutin$mut = mutin[,4]
string = mutin[nchar(unlist(mutin[,4])) != 1,4]
mutin$mut[nchar(unlist(mutin[,4])) != 1] = paste0('+',nchar(string)-1,substring(string, 2))
string = mutin[nchar(unlist(mutin[,3])) != 1,3]
mutin$mut[nchar(unlist(mutin[,3])) != 1] = paste0('-',nchar(string)-1,substring(string, 2))
mutin$location = paste0(mutin[,1],':',mutin[,2])
mutin$mut = paste0(':',mutin$mut)
mutin = mutin[,c(3,4,5,6)]
mutin = as.data.table(mutin)
#mutin[, V3 := substr(V3, 1, 1)]
setnames(mutin, c("WT","VAR",'mut','location'))

file_list <- list.files(path = VCFs_path, pattern = "*.vcf", full.names = TRUE)
dedupmut <- function(xx){
  sapply(xx, FUN=function(xx){
    positions = str_locate_all(xx, '(\\+|-)(\\d+)')[[1]]
    startpos = positions[,1]
    endpos = positions[,2]
    endlength = as.numeric(str_match_all(xx, "(\\+|-)(\\d+)")[[1]][,3])
    endpos = endpos + endlength
    indels <- substr(rep(xx, length(startpos)), startpos, endpos)
    deletpos <- unlist(mapply(function(a1,b1){a1:b1},startpos,endpos,SIMPLIFY = F))
    if (length(deletpos) >= 1) {
      noindels <- strsplit(xx,'')[[1]][-deletpos]
      pileups = c(indels,noindels)
    }else{
      pileups <- strsplit(xx,'')[[1]]
    }
    pileups = gsub('.', ',', pileups,fixed = TRUE)
    pileups = paste0(':',paste0((pileups),collapse = ':'),':')
    return(pileups)
  })
}
dedupmut_noindel <- function(xx) {
  xx = gsub('.', ',', xx,fixed = TRUE)
  paste0(':',sapply(xx,function(xx){paste0((strsplit(xx,'')[[1]]),collapse = ':')}),':')
}

read_and_process_file <- function(file_path) {
  data <- fread(file_path, header = FALSE)
  if (ncol(data) == 0) {
    data <- data.table(location = character(), file_path = numeric())
    setnames(data, c("location", file_path))
    return(data)
  }

  setnames(data, c("V1", "V2", "V3", "V4", "V5", "V6"))
  data[, V5 := gsub('\\^.','',V5)]
  data[, V5 := toupper(gsub('(>|<|\\$)','',V5))]
  data = data[V5 != '',]
  data <- data[, .(V1, V2, V5)]
  data[grepl('(\\+|\\-)', V5), V5 := dedupmut(V5)]
  data[!grepl('(\\+|\\-)', V5), V5 := dedupmut_noindel(V5)]
  data$location = paste0(data$V1,':',data$V2)
  data = merge(mutin,data,by = 'location')
  data[, Mutation := 2*as.numeric(str_count(pattern = fixed(mut),V5) > 0)+as.numeric(grepl(pattern = ':,',x = V5,fixed=TRUE))]
  data[, location := paste(location, WT, VAR, sep = ":")]
  data <- data[, .(location, Mutation)]
  setnames(data, c("location", file_path))
  data <- data[!duplicated(location)]
  return(data)
}


data_list <- mclapply(file_list, read_and_process_file, mc.cores = num_cores)
n = floor(sqrt(length(data_list)))
data_list = split(data_list, cut(seq_along(data_list), breaks = n, labels = FALSE))
data_list <- mclapply(data_list, function(ii){Reduce(function(x, y) merge(x, y, by = 'location', all = TRUE), ii) },mc.cores = num_cores)
merged_data <- Reduce(function(x, y) merge(x, y, by = 'location', all = TRUE), data_list)
fwrite(merged_data,outputMat,sep = "\t")





read_and_process_file0 <- function(file_path) {
    
  data <- fread(file_path, header = FALSE)
  if (ncol(data) == 0) {
    data <- data.table(location = character(), file_path = numeric())
    setnames(data, c("location", file_path))
    return(data)
  }
  
  setnames(data, c("V1", "V2", "V3", "V4", "V5", "V6"))
  data = data[V4 != 0,]
  data[, V5 := toupper(gsub('(>|<|\\^|\\$)','',V5))]
  data = data[V5 != '',]
  data <- data[, .(V1, V2, V5)]
  data[grepl('(\\+|\\-)', V5), V5 := dedupmut(V5)]
  data[!grepl('(\\+|\\-)', V5), V5 := dedupmut_noindel(V5)]
  data[, Normal := str_count(string = V5,pattern = fixed(':,'))]
  data$location = paste0(data$V1,':',data$V2)
  data = merge(mutin,data,by = 'location')
  data[, Mutation := str_count(string = V5,pattern = fixed(mut))]
  data[, Mut := str_count(string = V5,pattern = fixed(':,'))]
  data = data.table::melt(data, measure.vars = c("Normal", "Mutation"), variable.name = "source_col")
  data[source_col == "Normal", location := paste(location, WT, sep = "")]
  data[source_col == "Mutation", location := paste(location, mut, sep = "")]
  data <- data[, .(location, value)]
  setnames(data, c("location", file_path))
  data <- data[!duplicated(location)]
  return(data)
}
