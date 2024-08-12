
make sure those software is avaiable in your environment: cnvkit
cnvkit.py batch result/LM0623_map/DNA/Tumor/dedup_tmr.bam --normal result/LM0623_map/DNA/Normal/dedup_nor.bam --fasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --output-reference result/CNV/cnvkit/my_reference.cnn --output-dir result/CNV/cnvkit/ --method wgs --diagram --scatter -p 16
cnvkit.py scatter -s result/CNV/cnvkit/dedup_tmr.cn{s,r} -c 2:47160083-61051990  -c 2:47160083-61051990  -o result/CNV/cnvkit/2_47160083_61051990.pdf
cnvkit.py scatter -s result/CNV/cnvkit/dedup_tmr.cn{s,r} -c 6:10492222-21232404  -c 6:10492222-21232404  -o result/CNV/cnvkit/6_10492222_21232404.pdf
cnvkit.py scatter -s result/CNV/cnvkit/dedup_tmr.cn{s,r} -c 16:10760918-28063714 -c 16:10760918-28063714 -o result/CNV/cnvkit/16_10760918_28063714.pdf


perl scripts/CNV/genemat2tsv4cnv.pl -o expmat4cnv.tsv -d spatial_barcodes100.txt -i expmat.bed
cat gencode.v45.annotation.gtf | egrep -v '^#' | awk '$3 == "gene"' | sed -E 's/\t[^\t]*gene_name "/\t/' | sed -E 's/";.+$//'  | cut -f 1,4,5,9 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $4, $1, $2, $3}'| awk '!seen[$1]++' "$input_file" >  result/CNV/inferCNV/transloc.txt
sed -n '2,$p' Cluster_and_spatial.id_LM0623_Large.csv | sed 's/"//g' | sed 's/,/\t/g' > result/CNV/inferCNV/cluster_id.csv

# those software were used in the environment: R(4.3.0) JAGS
Rscript scripts/CNV/inferCNV.R
Rscript scripts/CNV/deal_after_inferCNV.R