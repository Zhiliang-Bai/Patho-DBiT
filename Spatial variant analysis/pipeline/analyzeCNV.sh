#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=analyzeCNV
#SBATCH -o cache/analyzeCNV.%j.out
#SBATCH --ntasks=16 --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=23:00:00
#SBATCH --mail-type=END



#module load miniconda
#conda activate cnvkit
#cnvkit.py batch result/LM0623_map/DNA/Tumor/dedup_tmr.bam --normal result/LM0623_map/DNA/Normal/dedup_nor.bam --fasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --output-reference result/CNV/cnvkit/my_reference.cnn --output-dir result/CNV/cnvkit/ --method wgs --diagram --scatter -p 16
#cnvkit.py scatter -s result/CNV/cnvkit/dedup_tmr.cn{s,r} -c 2:47160083-61051990  -c 2:47160083-61051990  -o result/CNV/cnvkit/2_47160083_61051990.pdf
#cnvkit.py scatter -s result/CNV/cnvkit/dedup_tmr.cn{s,r} -c 6:10492222-21232404  -c 6:10492222-21232404  -o result/CNV/cnvkit/6_10492222_21232404.pdf
#cnvkit.py scatter -s result/CNV/cnvkit/dedup_tmr.cn{s,r} -c 16:10760918-28063714 -c 16:10760918-28063714 -o result/CNV/cnvkit/16_10760918_28063714.pdf


#perl scripts/CNV/genemat2tsv4cnv.pl -o result/CNV/inferCNV/expmat4cnv.tsv -d /home/dz287/gibbs/scRNAmir/data/FFPE/spatial_barcodes100.txt -i ~/Mark/WHOLEres/LM0623P100/expmat.bed
#cat /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/ncGTFs/gencode/gencode.v45.annotation.gtf | egrep -v '^#' | awk '$3 == "gene"' | sed -E 's/\t[^\t]*gene_name "/\t/' | sed -E 's/";.+$//'  | cut -f 1,4,5,9 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $4, $1, $2, $3}'| awk '!seen[$1]++' "$input_file" >  result/CNV/inferCNV/transloc.txt
#sed -n '2,$p' /home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/Cluster_and_spatial.id_LM0623_Large.csv | sed 's/"//g' | sed 's/,/\t/g' > result/CNV/inferCNV/cluster_id.csv
module purge
module load R/4.3.0-foss-2020b JAGS
Rscript scripts/CNV/inferCNV.R
Rscript scripts/CNV/deal_after_inferCNV.R