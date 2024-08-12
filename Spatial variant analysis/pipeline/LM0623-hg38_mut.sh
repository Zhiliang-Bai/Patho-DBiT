#make sure SAMtools Strelka STAR bcftools are available


perl scripts/bam4mut.pl -o result/LM0623_map/RNA/bam4mut.bam -i /home/dz287/Mark/WHOLEres/LM0623P100/STAR/tempfiltered.bam
wait

perl scripts/split_fa_cluster.pl -r /home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/Cluster_and_spatial.csv -i result/LM0623_map/RNA/bam4mut.bam -o ./ -n 0:1:2:3:4:5:8:9:10:11:12:13:15:16:17:19 -t 6:7:14:18
mv tumor.bam* result/LM0623_map/RNA/splitedBams/
mv normal.bam* result/LM0623_map/RNA/splitedBams/
perl pipeline/VCF2Mut/split_fa_pixel.pl -i result/LM0623_map/RNA/bam4mut.bam -o result/LM0623_map/RNA/splitedBams/pixels/


configureStrelkaGermlineWorkflow.py --bam result/bam4mut.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --rna --runDir result/LM0623_map/RNA/germline
result/LM0623_map/RNA/germline/runWorkflow.py -m local -j 16
mv result/LM0623_map/RNA/germline/ result/Strelka/RNA/germline/runWorkflow.py

configureStrelkaGermlineWorkflow.py --bam result/LM0623_map/DNA/Tumor/dedup.bam --bam result/LM0623_map/DNA/Normal/dedup.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --runDir result/Strelka/DNA/germline/
result/Strelka/DNA/germline/runWorkflow.py -m local -j 16

configureStrelkaSomaticWorkflow.py --tumorBam result/LM0623_map/RNA/splitedBams/tumor.bam --normalBam result/LM0623_map/RNA/splitedBams/normal.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --runDir result/Strelka/RNA/somatic
result/Strelka/RNA/somatic/runWorkflow.py -m local -j 16
configureStrelkaSomaticWorkflow.py --tumorBam result/LM0623_map/DNA/Tumor/dedup.bam --normalBam result/LM0623_map/DNA/Normal/dedup.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --runDir result/Strelka/DNA/somatic
result/Strelka/DNA/somatic/runWorkflow.py -m local -j 16





perl scripts/split_fa_cluster.pl -r cluster_and_position_by_cnv.csv -i result/LM0623_map/RNA/bam4mut.bam -o result/LM0623_map/RNA/splitedBams/SNVclusters/SNVcluster
configureStrelkaGermlineWorkflow.py --bam result/LM0623_map/RNA/splitedBams/normal_SNVcluster.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --rna --runDir result/Strelka/RNA/normal_SNVcluster_germline/
result/Strelka/RNA/normal_SNVcluster_germline/runWorkflow.py -m local -j 16

configureStrelkaSomaticWorkflow.py --tumorBam result/LM0623_map/RNA/splitedBams/tumor_SNVcluster.bam --normalBam result/LM0623_map/RNA/splitedBams/normal_SNVcluster.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --runDir result/Strelka/RNA/somatic_SNVcluster
result/Strelka/RNA/somatic_SNVcluster/runWorkflow.py -m local -j 16

for i in {0..8};
do
    configureStrelkaGermlineWorkflow.py --bam result/LM0623_map/RNA/splitedBams/SNVclusters/SNVcluster${i}.bam --referenceFasta Homo_sapiens.GRCh38.dna.fa --rna --runDir result/Strelka/RNA/SNVclusters/${i}/
    result/Strelka/RNA/SNVclusters/${i}/runWorkflow.py -m local -j 16
done

for i in {0..8};
do
    VCFin=result/Strelka/RNA/SNVclusters/${i}/results/variants/variants.vcf.gz
    VCFout=result/Figures/VCFs4Figs/filteredSNVclusters${i}.vcf
    bcftools filter -i 'max(AD[0:1-]) > 10' $VCFin > $VCFout &
done
wait

bcftools filter -i 'max(AD[0:1-]) > 10' result/Strelka/RNA/normal_SNVcluster_germline/results/variants/variants.vcf.gz > result/Figures/VCFs4Figs/filteredSNVcluster_Normal.vcf












