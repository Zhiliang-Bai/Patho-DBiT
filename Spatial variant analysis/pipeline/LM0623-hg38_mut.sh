### those software were used in the environment:  Strelka vep SAMtools STAR bcftools
### pipeline/VCF2Mut/VCF2mut.sh is a pipeline to get variation-mutation matrix from bam files

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
perl pipeline/VCF2Mut/split_fa_pixel.pl -i bam4mut.bam -o Sample_name/splited_RNA_Bams/pixels/


configureStrelkaGermlineWorkflow.py --bam bam4mut.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --rna --runDir RNA_germline/
RNA_germline/runWorkflow.py -m local -j 16


configureStrelkaSomaticWorkflow.py --tumorBam dedup_tmr.bam --normalBam dedup_nor.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --runDir DNA_somatic/
DNA_somatic/runWorkflow.py -m local -j 16

#SampleBam.list includes the names of tumor bam and normal bam from DNAseq and RNAseq
# like:
#result/DNA_Normal.bam
#result/DNA_Tumor.bam
#result/RNA_normal.bam
#result/RNA_tumor.bam

bash pipeline/VCF2Mut/VCF2mut.sh bcd VarMAT_RNAgermline/ SampleBam.list RNA_germline/results/variants/variants.vcf.gz Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterGermlineMut-RNA.R
bash pipeline/VCF2Mut/VCF2mut.sh cd VarMAT_FiltGermlineRNAvar/ Sample_name/splited_RNA_Bams/pixels/ VarMAT_FiltGermlineRNAvar/PASS.site.pos.txt Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
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

bcftools filter -i 'max(AD[0:1-]) > 10' normal_SNVcluster_germline/results/variants/variants.vcf.gz > result/Figures/VCFs4Figs/filteredSNVcluster_Normal.vcf




zcat somatic_SNVcluster/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq > VarMAT_RNAvar_SNVcluster/PASS.site.pos.txt
bash pipeline/VCF2Mut/VCF2mut.sh cd VarMAT_RNAvar_SNVcluster/ SampleBam.list VarMAT_RNAvar_SNVcluster/PASS.site.pos.txt Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/FilterSomaticMut-RNA_SNVcluster.R
bash pipeline/VCF2Mut/VCF2mut.sh cd VarMAT_FiltRNAvar_SNVcluster/ Sample_name/splited_RNA_Bams/pixels/ VarMAT_FiltRNAvar_SNVcluster/PASS.site.pos.txt Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/




# to get variant-sample/pixel count matrix
zcat DNA_somatic/results/variants/somatic.snvs.vcf.gz | egrep -v '^#' |  awk '{print $1 ":" $2 ":" $4 ":" $5}' | sort | uniq > VarMAT_DNAvar/PASS.site.pos.txt
bash pipeline/VCF2Mut/VCF2mut.sh cd VarMAT_DNAvar/ SampleBam.list VarMAT_DNAvar/PASS.site.pos.txt Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
Rscript scripts/ChisqFilterSomaticMut-DNA.R
bash pipeline/VCF2Mut/VCF2mut.sh cd VarMAT_ChisqFiltDNAvar/ Sample_name/splited_RNA_Bams/pixels/ VarMAT_ChisqFiltDNAvar/PASS.site.pos.txt Homo_sapiens.GRCh38.dna.fa pipeline/VCF2Mut/
conda activate perlr
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2,".",$3,$4}' VarMAT_ChisqFiltDNAvar/PASS.site.pos.txt | sort -k1,1 -k2,2n > VarMAT_ChisqFiltDNAvar/PASS.site.pos.txt.t.vcf
vep --dir data/vepcache/ --dir_cache data/vepcache/ --cache --offline --use_given_ref --fasta Homo_sapiens.GRCh38.dna.fa --symbol --force_overwrite --dont_skip --use_given_ref --cache_version 112 -i VarMAT_ChisqFiltDNAvar/PASS.site.pos.txt.t.vcf -o VarMAT_ChisqFiltDNAvar/SiteAnno.txt



