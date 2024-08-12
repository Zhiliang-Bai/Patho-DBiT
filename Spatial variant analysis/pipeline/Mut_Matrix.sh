###########
#cat result/LM0623_map/SampleBam.list
#result/LM0623_map/DNA/Normal/dedup_nor.bam
#result/LM0623_map/DNA/Tumor/dedup_tmr.bam
#result/LM0623_map/RNA/splitedBams/normal.bam
#result/LM0623_map/RNA/splitedBams/tumor.bam
###########


#筛选dna算出来的somatic mutation并且证明它在rna存在有意义的spatial information
zcat result/Strelka/DNA/somatic/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq > result/VarMAT/Sample/DNAvar/PASS.site.pos.txt
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Sample/DNAvar/ result/LM0623_map/SampleBam.list result/VarMAT/Sample/DNAvar/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterSomaticMut-DNA.R
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/FiltDNAvar/ result/LM0623_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/FiltDNAvar/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
conda activate perlr
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2,".",$3,$4}' result/VarMAT/Pixel/FiltDNAvar/PASS.site.pos.txt | sort -k1,1 -k2,2n > result/VarMAT/Pixel/FiltDNAvar/PASS.site.pos.txt.t.vcf
vep --dir data/vepcache/ --dir_cache data/vepcache/ --cache --offline --use_given_ref --fasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --symbol --force_overwrite --dont_skip --use_given_ref --cache_version 112 -i result/VarMAT/Pixel/FiltDNAvar/PASS.site.pos.txt.t.vcf -o result/VarMAT/Pixel/FiltDNAvar/SiteAnno.txt

Rscript scripts/StrictFilterSomaticMut-DNA.R
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/StrictFiltDNAvar/ result/LM0623_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/StrictFiltDNAvar/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
conda activate perlr
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2,".",$3,$4}' result/VarMAT/Pixel/StrictFiltDNAvar/PASS.site.pos.txt | sort -k1,1 -k2,2n > result/VarMAT/Pixel/StrictFiltDNAvar/PASS.site.pos.txt.t.vcf
vep --dir data/vepcache/ --dir_cache data/vepcache/ --cache --offline --use_given_ref --fasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --symbol --force_overwrite --dont_skip --use_given_ref --cache_version 112 -i result/VarMAT/Pixel/StrictFiltDNAvar/PASS.site.pos.txt.t.vcf -o result/VarMAT/Pixel/StrictFiltDNAvar/SiteAnno.txt

Rscript scripts/ChisqFilterSomaticMut-DNA.R
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/ChisqFiltDNAvar/ result/LM0623_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
conda activate perlr
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2,".",$3,$4}' result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt | sort -k1,1 -k2,2n > result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt.t.vcf
vep --dir data/vepcache/ --dir_cache data/vepcache/ --cache --offline --use_given_ref --fasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --symbol --force_overwrite --dont_skip --use_given_ref --cache_version 112 -i result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt.t.vcf -o result/VarMAT/Pixel/ChisqFiltDNAvar/SiteAnno.txt




#筛选rna算出来的somatic mutation并且证明它在dna也存在
zcat result/Strelka/RNA/somatic/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq > result/VarMAT/Sample/RNAvar/PASS.site.pos.txt
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Sample/RNAvar/ result/LM0623_map/SampleBam.list result/VarMAT/Sample/RNAvar/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterSomaticMut-RNA.R
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/FiltRNAvar/ result/LM0623_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/FiltRNAvar/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
conda activate perlr
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2,".",$3,$4}' result/VarMAT/Pixel/FiltRNAvar/PASS.site.pos.txt | sort -k1,1 -k2,2n > result/VarMAT/Pixel/FiltRNAvar/PASS.site.pos.txt.t.vcf
vep --dir data/vepcache/ --dir_cache data/vepcache/ --cache --offline --use_given_ref --fasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --symbol --force_overwrite --dont_skip --use_given_ref --cache_version 112 -i result/VarMAT/Pixel/FiltRNAvar/PASS.site.pos.txt.t.vcf -o result/VarMAT/Pixel/FiltRNAvar/SiteAnno.txt



####这个是看通过整个rna bam 文件进行germline call得到 variation 再通过这个variation提取四个bam样本里的mutation分布矩阵
####根据这个矩阵我们可以进一步筛选可信度高的rna derived germline mutation
sbatch pipeline/VCF2Mut/VCF2mut.sh bcd result/VarMAT/Sample/RNAgermline/ result/LM0623_map/SampleBam.list result/Strelka/RNA/germline/results/variants/variants.vcf.gz /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterGermlineMut-RNA.R
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/FiltGermlineRNAvar/ result/LM0623_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/FiltGermlineRNAvar/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2,".",$3,$4}' result/VarMAT/Pixel/FiltGermlineRNAvar/PASS.site.pos.txt | sort -k1,1 -k2,2n > result/VarMAT/Pixel/FiltGermlineRNAvar/PASS.site.pos.txt.t.vcf
vep --dir data/vepcache/ --dir_cache data/vepcache/ --cache --offline --use_given_ref --fasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --symbol --force_overwrite --dont_skip --use_given_ref --cache_version 112 -i result/VarMAT/Pixel/FiltGermlineRNAvar/PASS.site.pos.txt.t.vcf -o result/VarMAT/Pixel/FiltGermlineRNAvar/SiteAnno.txt



#利用RNA germline mmutation 得到的cluster来 算rna的somatic mutation并且证明它在dna也存在
zcat result/Strelka/RNA/somatic_SNVcluster/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq > result/VarMAT/Sample/RNAvar_SNVcluster/PASS.site.pos.txt
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Sample/RNAvar_SNVcluster/ result/LM0623_map/SampleBam.list result/VarMAT/Sample/RNAvar_SNVcluster/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterSomaticMut-RNA_SNVcluster.R

sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/FiltRNAvar_SNVcluster/ result/LM0623_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/FiltRNAvar_SNVcluster/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
conda activate perlr 
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2,".",$3,$4}' result/VarMAT/Pixel/FiltRNAvar_SNVcluster/PASS.site.pos.txt | sort -k1,1 -k2,2n > result/VarMAT/Pixel/FiltRNAvar_SNVcluster/PASS.site.pos.txt.t.vcf
vep --dir data/vepcache/ --dir_cache data/vepcache/ --cache --offline --use_given_ref --fasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --symbol --force_overwrite --dont_skip --use_given_ref --cache_version 112 -i result/VarMAT/Pixel/FiltRNAvar_SNVcluster/PASS.site.pos.txt.t.vcf -o result/VarMAT/Pixel/FiltRNAvar_SNVcluster/SiteAnno.txt





#这是算hd837的mutation matrix
STAR --genomeDir data/Genomes/GRCh38/StarIndex/ --readFilesIn /home/dz287/Mark/WHOLEres/HD837/oldmir/combine.fq --outFileNamePrefix /home/dz287/Mark/WHOLEres/HD837/STAR/temp --sjdbGTFfile data/Genomes/GRCh38/ncGTFs/hsa.all.gtf --runThreadN 16 --outSAMattributes NH HI AS nM NM --genomeLoad NoSharedMemory --limitOutSAMoneReadBytes 200000000 --outFilterMultimapNmax -1 --outFilterMultimapScoreRange 0 --readMatesLengthsIn NotEqual --limitBAMsortRAM 0 --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMorder Paired --outSAMprimaryFlag AllBestScore --outSAMmultNmax -1 --outFilterType Normal --outFilterScoreMinOverLread 0 --alignSJDBoverhangMin 30 --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0 --outFilterMismatchNoverLmax 0.1 --outFilterMismatchNoverReadLmax 0.15 --alignIntronMin 20 --alignIntronMax 1000000 --alignEndsType Local --twopassMode Basic
module load miniconda
conda activate mut
configureStrelkaGermlineWorkflow.py --bam /home/dz287/Mark/WHOLEres/HD837/STAR/tempAligned.sortedByCoord.out.bam --referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --rna --runDir result/Strelka/HD837_RNA/germline
result/Strelka/HD837_RNA/germline/runWorkflow.py -m local -j 16
bcftools filter -i 'max(AD[0:1-]) > 10' result/Strelka/HD837_RNA/germline/results/variants/variants.vcf.gz > result/VarMAT/Pixel/HD837_GermlineRNA/PASS.vcf
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/HD837_GermlineRNA/ result/HD837_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/HD837_GermlineRNA/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/











################################
################################
################################
################################
################################
################################
################################
################################
################################


####这个是用dna和rna进行somatic variant call得到的overlap的mutation list后再根据这个mutation list 提取 pixel mutation matrix
grep -f <(zcat result/Strelka/RNA/somatic/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq) <(zcat result/Strelka/DNA/somatic/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq) > result/VarMAT/Pixel/DNAoverlapRNA_SomMut/PASS.site.pos.txt
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/DNAoverlapRNA_SomMut/ result/LM0623_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/DNAoverlapRNA_SomMut/PASS.site.pos.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/


####这个提取了所有pass的dna somatic mutation
sbatch pipeline/VCF2Mut/VCF2mut.sh bcd result/VarMAT/Pixel/DNA_SomMut/ result/LM0623_map/RNA/splitedBams/pixels/ result/Strelka/DNA/somatic/results/variants/somatic.snvs.vcf.gz /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/

#这个好像没有实装
sbatch pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Sample/germline_SNVclusteris/ result/LM0623_map/SampleBam_SNVclusters.list result/VarMAT/Sample/DNAvar/PASS.site.pos.txt result/VarMAT/Sample/DNAvar/pixelPile/ pipeline/VCF2Mut/

