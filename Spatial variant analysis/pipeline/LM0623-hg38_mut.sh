#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=hg38LM0623mut
#SBATCH --output=cache/hg38LM0623mut.txt
#SBATCH --ntasks=1 --cpus-per-task=16
#SBATCH --mem=55g
#SBATCH --time=23:00:00



module load SAMtools
module load miniconda
conda activate mut

#这是之前没有对dedup.bam重命名的时候

perl scripts/bam4mut.pl -o result/LM0623_map/RNA/bam4mut.bam -i /home/dz287/Mark/WHOLEres/LM0623P100/STAR/tempfiltered.bam
samtools view -@ 16 -F 1024 -b result/LM0623_map/DNA/Tumor/marked.bam >   result/LM0623_map/DNA/Tumor/dedup.bam  &
samtools view -@ 16 -F 1024 -b result/LM0623_map/DNA/Normal/marked.bam >  result/LM0623_map/DNA/Normal/dedup.bam
wait
samtools index result/LM0623_map/DNA/Tumor/dedup.bam  &
samtools index result/LM0623_map/DNA/Normal/dedup.bam
wait

perl scripts/split_fa_cluster.pl -r /home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/Cluster_and_spatial.id_LM0623_Large.csv -i result/LM0623_map/RNA/bam4mut.bam -o ./ -n 0:1:2:3:4:5:8:9:10:11:12:13:15:16:17:19 -t 6:7:14:18
mv tumor.bam* result/LM0623_map/RNA/splitedBams/
mv normal.bam* result/LM0623_map/RNA/splitedBams/
perl pipeline/VCF2Mut/split_fa_pixel.pl -i result/LM0623_map/RNA/bam4mut.bam -o result/LM0623_map/RNA/splitedBams/pixels/


configureStrelkaGermlineWorkflow.py --bam result/bam4mut.bam --referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --rna --runDir result/LM0623_map/RNA/germline
result/LM0623_map/RNA/germline/runWorkflow.py -m local -j 16
mv result/LM0623_map/RNA/germline/ result/Strelka/RNA/germline/runWorkflow.py

configureStrelkaGermlineWorkflow.py --bam result/LM0623_map/DNA/Tumor/dedup.bam --bam result/LM0623_map/DNA/Normal/dedup.bam --referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --runDir result/Strelka/DNA/germline/
result/Strelka/DNA/germline/runWorkflow.py -m local -j 16

configureStrelkaSomaticWorkflow.py --tumorBam result/LM0623_map/RNA/splitedBams/tumor.bam --normalBam result/LM0623_map/RNA/splitedBams/normal.bam --referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --runDir result/Strelka/RNA/somatic
result/Strelka/RNA/somatic/runWorkflow.py -m local -j 16
configureStrelkaSomaticWorkflow.py --tumorBam result/LM0623_map/DNA/Tumor/dedup.bam --normalBam result/LM0623_map/DNA/Normal/dedup.bam --referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --runDir result/Strelka/DNA/somatic
result/Strelka/DNA/somatic/runWorkflow.py -m local -j 16


#现在对dedup.bam重命名了已经
mv result/LM0623_map/DNA/Tumor/dedup.bam result/LM0623_map/DNA/Tumor/dedup_tmr.bam
mv result/LM0623_map/DNA/Normal/dedup.bam result/LM0623_map/DNA/Normal/dedup_nor.bam
mv result/LM0623_map/DNA/Tumor/dedup.bam.bai result/LM0623_map/DNA/Tumor/dedup_tmr.bam.bai
mv result/LM0623_map/DNA/Normal/dedup.bam.bai result/LM0623_map/DNA/Normal/dedup_nor.bam.bai


###这回是用snv分的tumor和normal cluster
perl scripts/split_fa_cluster.pl -r /home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/cluster_and_position_LM0623_Large_RNAsomatic.csv -i result/LM0623_map/RNA/bam4mut.bam -o ./ -n 0:2:4:5:6:7:8 -t 1:3
mv tumor.bam result/LM0623_map/RNA/splitedBams/tumor_SNVcluster.bam
mv tumor.bam.bai result/LM0623_map/RNA/splitedBams/tumor_SNVcluster.bam.bai
mv normal.bam result/LM0623_map/RNA/splitedBams/normal_SNVcluster.bam
mv normal.bam.bai result/LM0623_map/RNA/splitedBams/normal_SNVcluster.bam.bai
perl scripts/split_fa_cluster.pl -r /home/dz287/gibbs/scRNAmir/data/FFPE/Clusters/cluster_and_position_LM0623_Large_RNAsomatic.csv -i result/LM0623_map/RNA/bam4mut.bam -o result/LM0623_map/RNA/splitedBams/SNVclusters/SNVcluster

configureStrelkaGermlineWorkflow.py --bam result/LM0623_map/RNA/splitedBams/normal_SNVcluster.bam --referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --rna --runDir result/Strelka/RNA/normal_SNVcluster_germline/
result/Strelka/RNA/normal_SNVcluster_germline/runWorkflow.py -m local -j 16

configureStrelkaSomaticWorkflow.py --tumorBam result/LM0623_map/RNA/splitedBams/tumor_SNVcluster.bam --normalBam result/LM0623_map/RNA/splitedBams/normal_SNVcluster.bam --referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --runDir result/Strelka/RNA/somatic_SNVcluster
result/Strelka/RNA/somatic_SNVcluster/runWorkflow.py -m local -j 16

for i in {0..8};
do
    configureStrelkaGermlineWorkflow.py --bam result/LM0623_map/RNA/splitedBams/SNVclusters/SNVcluster${i}.bam --referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --rna --runDir result/Strelka/RNA/SNVclusters/${i}/
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

module load miniconda
conda activate mut
configureStrelkaGermlineWorkflow.py --bam /home/dz287/Mark/WHOLEres/HD837/STAR/tempAligned.sortedByCoord.out.bam --referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa --rna --runDir result/Strelka/HD837_RNA/germline
result/Strelka/HD837_RNA/germline/runWorkflow.py -m local -j 16
bcftools filter -i 'max(AD[0:1-]) > 10' result/Strelka/HD837_RNA/germline/results/variants/variants.vcf.gz > result/VarMAT/Pixel/HD837_GermlineRNA/PASS.vcf
perl pipeline/VCF2Mut/split_fa_pixel_different_format.pl -i /home/dz287/Mark/WHOLEres/HD837/STAR/tempAligned.sortedByCoord.out.bam -o result/HD837_map/RNA/splitedBams/pixels/


















exit







#configureStrelkaGermlineWorkflow.py \
#--bam result/LM0623_map/RNA/bam4mut.bam \
#--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
#--rna  --runDir result/Strelka/RNA/germline
#result/Strelka/RNA/germline/runWorkflow.py -m local -j 16

#configureStrelkaGermlineWorkflow.py \
#--bam result/LM0623_STAR/splitedBams/tumor.bam \
#--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
#--rna  --runDir result/Strelka/RNA/tumor_germline
#result/Strelka/RNA/tumor_germline/runWorkflow.py -m local -j 16

#configureStrelkaGermlineWorkflow.py \
#--bam result/LM0623_STAR/splitedBams/normal.bam \
#--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
#--rna  --runDir result/Strelka/RNA/normal_germline
#result/Strelka/RNA/normal_germline/runWorkflow.py -m local -j 16

###############
###############



#configureStrelkaGermlineWorkflow.py \
#--bam result/LM0623/marked.bam \
#--bam result/Normal/marked.bam \
#--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
#--runDir result/Strelka/DNA/germline/
#result/Strelka/DNA/germline/runWorkflow.py -m local -j 16

#cat data/Homo_sapiens_assembly38.dbsnp138.vcf | sed 's/^chr//g' > data/Homo_sapiens_assembly38.dbsnp138.nochr.vcf
#configureStrelkaGermlineWorkflow.py \
#--bam result/LM0623/marked.bam \
#--bam result/Normal/marked.bam \
#--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
#--callRegions sorted_snv_hg38.bed.gz \
#--runDir result/Strelka/DNA/Candidate_Strelka
#result/Strelka/DNA/Candidate_Strelka/runWorkflow.py -m local -j 16

#configureStrelkaSomaticWorkflow.py --tumorBam result/LM0623/marked.bam --normalBam result/Normal/marked.bam  \
#--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
#--runDir result/Strelka/DNA/somatic
#result/Strelka/DNA/somatic/runWorkflow.py -m local -j 16

#configureStrelkaGermlineWorkflow.py \
#--bam result/LM0623/marked.bam \
#--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
#--runDir result/Strelka/DNA/tumor_germline/
#result/Strelka/DNA/tumor_germline/runWorkflow.py -m local -j 16

#configureStrelkaGermlineWorkflow.py \
#--bam result/Normal/marked.bam \
#--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
#--runDir result/Strelka/DNA/normal_germline/
#result/Strelka/DNA/normal_germline/runWorkflow.py -m local -j 16


configureStrelkaGermlineWorkflow.py \
--bam /home/dz287/Mark/WHOLEres/LM0623P100/STAR/tempAligned.sortedByCoord.out.bam \
--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
--rna  --runDir result/Strelka/RNA/germline1
result/Strelka/RNA/germline1/runWorkflow.py -m local -j 16

configureStrelkaGermlineWorkflow.py \
--bam /home/dz287/Mark/WHOLEres/LM0623P100/STAR/tempfiltered.bam \
--referenceFasta /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
--rna  --runDir result/Strelka/RNA/germline0
result/Strelka/RNA/germline0/runWorkflow.py -m local -j 16

