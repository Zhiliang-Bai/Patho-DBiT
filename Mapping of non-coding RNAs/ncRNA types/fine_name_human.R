transin <- read.table('name2type.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
transin = transin[!duplicated(transin[,2]),]

transin[,1] = gsub("^(transcribed_unprocessed_pseudogene|unprocessed_pseudogene|transcribed_processed_pseudogene|processed_pseudogene|transcribed_unitary_pseudogene|translated_processed_pseudogene|unitary_pseudogene|rRNA_pseudogene|IG_C_pseudogene|TR_J_pseudogene|TR_V_pseudogene||IG_C_pseudogene|TR_J_pseudogene|TR_V_pseudogene|IG_V_pseudogene|IG_pseudogene|IG_V_pseudogene|IG_J_pseudogene|IG_pseudogene)$",'pseudogene',transin[,1])
transin[,1] = gsub("^Mt_tRNA$",'tRNA',transin[,1])
transin[,1] = gsub("^(IG_V_gene|IG_C_gene|IG_J_gene|TR_C_gene|TR_J_gene|TR_V_gene|TR_D_gene|IG_D_gene)$",'protein_coding',transin[,1])
transin[,1] = gsub("^Mt_rRNA$",'rRNA',transin[,1])
transin[,1] = gsub("^TEC$",'misc_RNA',transin[,1])
transin[,1] = gsub("^ribozyme|artifact|sRNA$",'misc_RNA',transin[,1])

fine_nameit <- function(bedin,bedout){
  datain = fread(bedin)
  datain = as.data.frame(datain)
  #datain <- read.table(bedin, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
  gene_names0 = unique(datain[,3])
  gene_names = gene_names0
  double_dash = grep('--',gene_names)
  gene_names_start <- gsub('(.*--).*', '\\1', gene_names[double_dash], perl=TRUE)
  aa = gsub('.*','',gene_names)
  aa[double_dash] = gene_names_start
  gene_names_start = aa
  
  gene_names = gsub('.*--','',gene_names)
  types = gsub('^.*__','',gene_names)
  gene_names = gsub('__.*$','',gene_names)
  
  intronloc = types == "transcript"
  exonloc   = types == 'exon'
  modified_gene_types = transin[match(gene_names[intronloc | exonloc],transin[,2]),1]
  types[intronloc | exonloc] = modified_gene_types
  #types = gsub("^rRNA_pseudogene$",'pseudogene',types)
  
  intron_names = unique(gene_names[intronloc])
  exon_names = unique(gene_names[exonloc])
  short_genes = union(setdiff(exon_names,intron_names),setdiff(intron_names,exon_names))
  
  exoninron = ifelse(exonloc,'__exon','')
  exoninron[intronloc] = '__intron'
  exoninron[gene_names %in% short_genes] = ''
  
  outgene_name = paste(gene_names_start,gene_names,'__',types,exoninron,sep = '')
  outgene_name = cbind(gene_names0,outgene_name)
  outgene_name = outgene_name[match(datain[,3],outgene_name[,1]),2]
  datain[,3] = outgene_name
  write.table(datain, file = bedout, row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE)
}
fine_nameit('expmat.bed','fined_expmat.bed')


