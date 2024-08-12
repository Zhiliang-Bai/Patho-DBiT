#processing RNAcentral gtf:
########################################
########################################
# make sure bedtools is available for your envirnoment, it is needed for some scripts, if not, there will be error.
# If there is no error without bedtools, it means the scripts you uesd do not need bedtools


### this download maybe change with RNAcentral version, you can see the version at RNAcentral/release_notes.txt
wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz
###
mv homo_sapiens.GRCh38.gff3.gz RNAcentral/human/homo_sapiens.GRCh38.gff3
cat RNAcentral/human/homo_sapiens.GRCh38.gff3 | grep $'\t'transcript$'\t' | egrep -v '^[^0-9XY]+[\t]' | sed 's/^/chr/' > RNAcentral/human/homo_sapiens.GRCh38-addchr_findtrans.gff3
python3 div_clpsGTF.py -f RNAcentral/human/ -i RNAcentral/human/homo_sapiens.GRCh38-addchr_findtrans.gff3
########################################
########################################







#processing gencode gtf:
########################################
########################################
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
mv gencode.v45.annotation.gtf gencode/
cat gencode/gencode.v45.annotation.gtf | egrep -v '^#' | awk '/gene_name "Y_RNA"/ {gsub(/gene_type/,  "old_type"); gsub(/gene_name/,  "gene_type"); gsub(/gene_id/,  "gene_name"); print; next} 1' > gencode/gencode.v45.formod.gtf
########################################
########################################











#processing GtRNAdb gtf:
# download GtRNAdb/hg38-tRNAs.fa as GtRNAdb/Screenshot 2024-05-01 at 9.54.45 PM.png ,that is,
#wget https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa
# mv hg38-tRNAs.fa GtRNAdb/hg38-tRNAs.fa.txt
# but sometimes the wget does not work so just download it form browser manually
########################################
########################################
cat GtRNAdb/hg38-tRNAs.fa.txt | grep '>' | awk 'BEGIN{OFS="\t"} {split($(NF-1), parts, "[:-]"); gsub(/[()]/, "", $NF); gsub(/>/, "", $1); print parts[1], "ENSEMBL", "exon", parts[2], parts[3], ".", $NF, ".", $1"__tRNA"  }' > GtRNAdb/hg38-tRNAs.mod.gtf
########################################
########################################







#modifying gtfs:
########################################
########################################
python3 modGTF.py -i gencode/gencode.v45.formod.gtf   -T gene_type -G gene_name -f gencode/gencode.v45.mod.gtf
python3 modGTF.py -i RNAcentral/human/vault_RNA.gtf   -T gene_type -G gene_name -f RNAcentral/human/vault_RNA.mod.gtf
python3 modGTF.py -i RNAcentral/human/Y_RNA.gtf   -T gene_type -G gene_name -f RNAcentral/human/Y_RNA.mod.gtf
########################################
########################################







#final gtfs
########################################
########################################
python3 clpsGTF.py -i GtRNAdb/hg38-tRNAs.mod.gtf:RNAcentral/human/vault_RNA.mod.gtf:RNAcentral/human/Y_RNA.mod.gtf:gencode/gencode.v45.mod.gtf -o hsa.all.gtf
########################################
########################################













