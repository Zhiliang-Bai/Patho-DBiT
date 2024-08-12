cat gencode.v45.annotation.gtf | grep $'\t'gene$'\t' | egrep -v '^#' | sed 's/^.*gene_type "//' | sed 's/";.*gene_name "/\t/' | sed 's/";.*$//' > name2type.txt
Rscript fine_name_human.R
#for each fined_expmat.bed 
perl  genemat2tsv.pl -i fined_expmat.bed  -o fined_expmat.tsv -d spatial_barcodes100.txt
cat Sample_name/fined_expmat.bed | cut -f 3 | sed -E 's/(__exon|__intron)$//g' | sed 's/.*__//g' | sort | uniq -c

cat gencode.vM34.annotation.gtf | grep $'\t'gene$'\t' | egrep -v '^#' | sed 's/^.*gene_type "//' | sed 's/";.*gene_name "/\t/' | sed 's/";.*$//' > name2type_mus.txt
Rscript fine_name_mouse.R
#for each fined_expmat.bed
perl  genemat2tsv.pl -i fined_expmat.bed  -o S_PAP_fine_expmat.tsv -d spatial_barcodes100.txt

cat Sample_name/fined_expmat.bed | cut -f 3 | sed -E 's/(__exon|__intron)$//g' | sed 's/.*__//g' | sort | uniq -c

