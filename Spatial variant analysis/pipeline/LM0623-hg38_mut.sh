### those software were used in the environment:  Strelka vep SAMtools STAR bcftools


# remove multiple mapping reads and unmapped reads
perl scripts/bam4mut.pl -o bam4mut.bam -i STAR_get.bam
wait

# An example of Cluster_and_spatial.csv
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

perl scripts/split_fa_cluster.pl -r Cluster_and_spatial.csv -i bam4mut.bam -o ./ -n 0:1:2:3:4:5:8:9:10:11:12:13:15:16:17:19 -t 6:7:14:18
mv tumor.bam* Sample_Tumor/splitedBams/
mv normal.bam* Sample_Tumor/splitedBams/
perl pipeline/VCF2Mut/split_fa_pixel.pl -i bam4mut.bam -o Sample_Tumor/splitedBams/pixels/


configureStrelkaGermlineWorkflow.py --bam bam4mut.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --rna --runDir RNA_germline/
RNA_germline/runWorkflow.py -m local -j 16


configureStrelkaSomaticWorkflow.py --tumorBam dedup_tmr.bam --normalBam dedup_nor.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --runDir DNA_somatic/
DNA_somatic/runWorkflow.py -m local -j 16

#SampleBam.list includes the names of tumor bam and normal bam from DNAseq and RNAseq
# like:
#DNA_Normal.bam
#DNA_Tumor.bam
#RNA_normal.bam
#RNA_tumor.bam

bash pipeline/VCF2Mut/VCF2mut.sh bcd VarMAT_RNAgermline/ SampleBam.list RNA_germline/results/variants/variants.vcf.gz Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterGermlineMut-RNA.R
bash pipeline/VCF2Mut/VCF2mut.sh cd VarMAT_FiltGermlineRNAvar/ Sample_Tumor/splitedBams/pixels/ VarMAT_FiltGermlineRNAvar/PASS.site.pos.txt Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
# unsupervised clustering using matrix in VarMAT_FiltGermlineRNAvar/ and get cluster_and_position_by_snv.csv
perl scripts/split_fa_cluster.pl -r cluster_and_position_by_snv.csv -i bam4mut.bam -o Sample_Tumor/splitedBams/SNVcluster
# get normal_SNVcluster.bam and tumor_SNVcluster.bam from normal clusters in Sample_Tumor/splitedBams/SNVcluster
configureStrelkaGermlineWorkflow.py --bam normal_SNVcluster.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --rna --runDir normal_SNVcluster_germline/
normal_SNVcluster_germline/runWorkflow.py -m local -j 16

configureStrelkaSomaticWorkflow.py --tumorBam tumor_SNVcluster.bam --normalBam normal_SNVcluster.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --runDir somatic_SNVcluster/
somatic_SNVcluster/runWorkflow.py -m local -j 16

for i in {0..8};
do
    configureStrelkaGermlineWorkflow.py --bam Sample_Tumor/splitedBams/SNVcluster/${i}.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --rna --runDir SNVclusters/${i}/
    SNVclusters/${i}/runWorkflow.py -m local -j 16
done

for i in {0..8};
do
    VCFin=SNVclusters/${i}/results/variants/variants.vcf.gz
    VCFout=result/Figures/VCFs4Figs/filteredSNVclusters${i}.vcf
    bcftools filter -i 'max(AD[0:1-]) > 10' $VCFin > $VCFout &
done
wait

bcftools filter -i 'max(AD[0:1-]) > 10' result/Strelka/RNA/normal_SNVcluster_germline/results/variants/variants.vcf.gz > result/Figures/VCFs4Figs/filteredSNVcluster_Normal.vcf




zcat result/Strelka/RNA/somatic_SNVcluster/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq > result/VarMAT/Sample/RNAvar_SNVcluster/PASS.site.pos.txt
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Sample/RNAvar_SNVcluster/ SampleBam.list result/VarMAT/Sample/RNAvar_SNVcluster/PASS.site.pos.txt Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterSomaticMut-RNA_SNVcluster.R
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/FiltRNAvar_SNVcluster/ result/Cancer_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/FiltRNAvar_SNVcluster/PASS.site.pos.txt Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/


# to get variant-sample/pixel count matrix
zcat result/Strelka/DNA/somatic/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq > result/VarMAT/Sample/DNAvar/PASS.site.pos.txt
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Sample/DNAvar/ SampleBam.list result/VarMAT/Sample/DNAvar/PASS.site.pos.txt scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/ChisqFilterSomaticMut-DNA.R
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/ChisqFiltDNAvar/ result/Cancer_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
conda activate perlr
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2,".",$3,$4}' result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt | sort -k1,1 -k2,2n > result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt.t.vcf
vep --dir data/vepcache/ --dir_cache data/vepcache/ --cache --offline --use_given_ref --fasta scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --symbol --force_overwrite --dont_skip --use_given_ref --cache_version 112 -i result/VarMAT/Pixel/ChisqFiltDNAvar/PASS.site.pos.txt.t.vcf -o result/VarMAT/Pixel/ChisqFiltDNAvar/SiteAnno.txt






# to get health variant-sample/pixel count matrix
STAR --genomeDir data/Genomes/GRCh38/StarIndex/ --readFilesIn health_combine.fq --outFileNamePrefix WHOLEres/health/STAR/temp --sjdbGTFfile data/Genomes/GRCh38/ncGTFs/hsa.all.gtf --runThreadN 16 --outSAMattributes NH HI AS nM NM --genomeLoad NoSharedMemory --limitOutSAMoneReadBytes 200000000 --outFilterMultimapNmax -1 --outFilterMultimapScoreRange 0 --readMatesLengthsIn NotEqual --limitBAMsortRAM 0 --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMorder Paired --outSAMprimaryFlag AllBestScore --outSAMmultNmax -1 --outFilterType Normal --outFilterScoreMinOverLread 0 --alignSJDBoverhangMin 30 --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0 --outFilterMismatchNoverLmax 0.1 --outFilterMismatchNoverReadLmax 0.15 --alignIntronMin 20 --alignIntronMax 1000000 --alignEndsType Local --twopassMode Basic
configureStrelkaGermlineWorkflow.py --bam WHOLEres/health/STAR/tempAligned.sortedByCoord.out.bam --referenceFasta scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --rna --runDir result/Strelka/health_RNA/germline
result/Strelka/health_RNA/germline/runWorkflow.py -m local -j 16
bcftools filter -i 'max(AD[0:1-]) > 10' result/Strelka/health_RNA/germline/results/variants/variants.vcf.gz > result/VarMAT/Pixel/health_GermlineRNA/PASS.vcf
bash pipeline/VCF2Mut/VCF2mut.sh cd result/VarMAT/Pixel/health_GermlineRNA/ result/health_map/RNA/splitedBams/pixels/ result/VarMAT/Pixel/health_GermlineRNA/PASS.site.pos.txt scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
















