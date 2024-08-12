#processing RNAcentral gtf:
########################################
########################################
### this download maybe change with RNAcentral version, you can see the version at RNAcentral/release_notes.txt
wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3/mus_musculus.GRCm39.gff3.gz
####
cat RNAcentral/mouse/mus_musculus.GRCm39.gff3 | grep $'\t'transcript$'\t' | egrep -v '^[^0-9XY]+[\t]' | sed 's/^/chr/' > RNAcentral/mouse/mus_musculus.GRCm39-addchr_findtrans.gff3
python3 div_clpsGTF.py -f RNAcentral/mouse/ -i RNAcentral/mouse/mus_musculus.GRCm39-addchr_findtrans.gff3
########################################
########################################







#processing gencode gtf:
########################################
########################################
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.annotation.gtf.gz
mv gencode.vM34.annotation.gtf.gz gencode
cat gencode/gencode.vM34.annotation.gtf | egrep -v '^#' | awk '/Gm55767|Gm56322|Gm56480|Gm56181|Gm55795|Gm56246|Gm56393|Gm55481|Gm54376|Gm54851/ {gsub(/misc_RNA/,  "Y_RNA"); print; next} 1' > gencode/gencode.vM34.formod.gtf

########################################
########################################











#processing GtRNAdb gtf:
# download GtRNAdb/mm39-tRNAs.fa as GtRNAdb/Screenshot 2024-08-11 at 7.01.37 PM.png ,that is,
#wget https://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/mm39-tRNAs.fa
# mv mm39-tRNAs.fa GtRNAdb/mm39-tRNAs.fa.txt
# but sometimes the wget does not work so just download it form browser manually

########################################
########################################
cat GtRNAdb/mm39-tRNAs.fa.txt | grep '>' | awk 'BEGIN{OFS="\t"} {split($(NF-1), parts, "[:-]"); gsub(/[()]/, "", $NF); gsub(/>/, "", $1); print parts[1], "ENSEMBL", "exon", parts[2], parts[3], ".", $NF, ".", $1"__tRNA"  }' > GtRNAdb/mm39-tRNAs.mod.gtf
########################################
########################################







#modifying gtfs:
########################################
########################################
python3 modGTF.py -i gencode/gencode.vM34.formod.gtf  -T gene_type -G gene_name -f gencode/gencode.vM34.mod.gtf
python3 modGTF.py -i RNAcentral/mouse/vault_RNA.gtf   -T gene_type -G gene_name -f RNAcentral/mouse/vault_RNA.mod.gtf
python3 modGTF.py -i RNAcentral/mouse/Y_RNA.gtf   -T gene_type -G gene_name -f RNAcentral/mouse/Y_RNA.mod.gtf
########################################
########################################







#final gtfs
########################################
########################################
python3 clpsGTF.py -i GtRNAdb/mm39-tRNAs.mod.gtf:gencode/gencode.vM34.mod.gtf:RNAcentral/mouse/vault_RNA.mod.gtf:RNAcentral/mouse/Y_RNA.mod.gtf -o mmu.all.gtf -C
########################################
########################################




