
make sure those software is avaiable in your environment: cnvkit
cnvkit.py batch dedup_tmr.bam --normal dedup_nor.bam --fasta Homo_sapiens.GRCh38.dna.fa --output-reference cnvkit/my_reference.cnn --output-dir cnvkit/ --method wgs --diagram --scatter -p 16
cnvkit.py scatter -s cnvkit/dedup_tmr.cn{s,r} -c 2:47160083-61051990  -c 2:47160083-61051990  -o cnvkit/2_47160083_61051990.pdf
cnvkit.py scatter -s cnvkit/dedup_tmr.cn{s,r} -c 6:10492222-21232404  -c 6:10492222-21232404  -o cnvkit/6_10492222_21232404.pdf
cnvkit.py scatter -s CNV/cnvkit/dedup_tmr.cn{s,r} -c 16:10760918-28063714 -c 16:10760918-28063714 -o cnvkit/16_10760918_28063714.pdf


perl scripts/CNV/genemat2tsv4cnv.pl -o expmat4cnv.tsv -d spatial_barcodes100.txt -i expmat.bed
cat gencode.v45.annotation.gtf | egrep -v '^#' | awk '$3 == "gene"' | sed -E 's/\t[^\t]*gene_name "/\t/' | sed -E 's/";.+$//'  | cut -f 1,4,5,9 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $4, $1, $2, $3}'| awk '!seen[$1]++' "$input_file" >  inferCNV/transloc.txt
# An example of Cluster_and_spatial.id.csv
#"","SCT_snn_res.1.2"
#"100x1","8"
#"100x2","4"
#"100x3","15"
#"100x4","15"
#"100x5","4"
#"100x6","3"
#"100x7","18"
#"100x8","8"
#"100x9","3"
sed -n '2,$p' Cluster_and_spatial.id.csv | sed 's/"//g' | sed 's/,/\t/g' > inferCNV/cluster_id.csv


# those software were used in the environment: R(4.3.0) JAGS
Rscript scripts/CNV/inferCNV.R
Rscript scripts/CNV/deal_after_inferCNV.R
