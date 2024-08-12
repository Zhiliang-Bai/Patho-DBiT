# those software were used in the environment:  Strelka vep

# to get variant-sample/pixel count matrix
zcat result/Strelka/DNA/somatic/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq > result/VarMAT/Sample/DNAvar/PASS.site.pos.txt
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Sample/DNAvar/ result/Cancer_map/SampleBam.list result/VarMAT/Sample/DNAvar/PASS.site.pos.txt scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/ChisqFilterSomaticMut-DNA.R
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/ChisqFiltDNAvar/ result/Cancer_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
conda activate perlr
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2,".",$3,$4}' result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt | sort -k1,1 -k2,2n > result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt.t.vcf
vep --dir data/vepcache/ --dir_cache data/vepcache/ --cache --offline --use_given_ref --fasta scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --symbol --force_overwrite --dont_skip --use_given_ref --cache_version 112 -i result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt.t.vcf -o result/VarMAT/Pixel/ChisqFiltDNAvar/SiteAnno.txt





bash pipeline/VCF2Mut/VCF2mut.sh bcd result/VarMAT/Sample/RNAgermline/ result/Cancer_map/SampleBam.list result/Strelka/RNA/germline/results/variants/variants.vcf.gz scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterGermlineMut-RNA.R
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/FiltGermlineRNAvar/ result/Cancer_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/FiltGermlineRNAvar/PASS.site.pos.txt scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
# unsupervised clustering using matrix in result/VarMAT/Pixel/FiltGermlineRNAvar/
zcat result/Strelka/RNA/somatic_SNVcluster/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq > result/VarMAT/Sample/RNAvar_SNVcluster/PASS.site.pos.txt
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Sample/RNAvar_SNVcluster/ result/Cancer_map/SampleBam.list result/VarMAT/Sample/RNAvar_SNVcluster/PASS.site.pos.txt scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterSomaticMut-RNA_SNVcluster.R
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/FiltRNAvar_SNVcluster/ result/Cancer_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/FiltRNAvar_SNVcluster/PASS.site.pos.txt scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/



# to get health variant-sample/pixel count matrix
STAR --genomeDir data/Genomes/GRCh38/StarIndex/ --readFilesIn health_combine.fq --outFileNamePrefix WHOLEres/health/STAR/temp --sjdbGTFfile data/Genomes/GRCh38/ncGTFs/hsa.all.gtf --runThreadN 16 --outSAMattributes NH HI AS nM NM --genomeLoad NoSharedMemory --limitOutSAMoneReadBytes 200000000 --outFilterMultimapNmax -1 --outFilterMultimapScoreRange 0 --readMatesLengthsIn NotEqual --limitBAMsortRAM 0 --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMorder Paired --outSAMprimaryFlag AllBestScore --outSAMmultNmax -1 --outFilterType Normal --outFilterScoreMinOverLread 0 --alignSJDBoverhangMin 30 --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0 --outFilterMismatchNoverLmax 0.1 --outFilterMismatchNoverReadLmax 0.15 --alignIntronMin 20 --alignIntronMax 1000000 --alignEndsType Local --twopassMode Basic
configureStrelkaGermlineWorkflow.py --bam WHOLEres/health/STAR/tempAligned.sortedByCoord.out.bam --referenceFasta scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --rna --runDir result/Strelka/health_RNA/germline
result/Strelka/health_RNA/germline/runWorkflow.py -m local -j 16
bcftools filter -i 'max(AD[0:1-]) > 10' result/Strelka/health_RNA/germline/results/variants/variants.vcf.gz > result/VarMAT/Pixel/health_GermlineRNA/PASS.vcf
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/health_GermlineRNA/ result/health_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/health_GermlineRNA/PASS.site.pos.txt scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/





