#!/bin/bash
#SBATCH -p pi_gerstein -A gerstein
#SBATCH --job-name=DNAmutation
#SBATCH -o cache/DNAmutation.%j.out
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=1-20:00:00


module load BWA GATK/3.8-1-0-gf15c1c3ef-Java-1.8 SAMtools



#cat data/Homo_sapiens_assembly38.dbsnp138.vcf | sed 's/^chr//g' > data/Homo_sapiens_assembly38.dbsnp138.nochr.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T MuTect2 \
-R /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
-I:tumor result/LM0623/marked.bam \
-I:normal result/Normal/marked.bam \
--dbsnp data/Homo_sapiens_assembly38.dbsnp138.nochr.vcf \
--contamination_fraction_to_filter 0.02 \
-o cache/mutect_variants.vcf \
--output_mode EMIT_VARIANTS_ONLY \
--disable_auto_index_creation_and_locking_when_reading_rods



module load BWA GATK/3.8-1-0-gf15c1c3ef-Java-1.8 SAMtools


#select potenital snv sites from RNAseq
#cat data/Homo_sapiens_assembly38.dbsnp138.vcf | sed 's/^chr//g' > data/Homo_sapiens_assembly38.dbsnp138.nochr.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T MuTect2 \
-R /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa \
-I:tumor result/LM0623/marked.bam \
-I:normal result/Normal/marked.bam \
--dbsnp data/Homo_sapiens_assembly38.dbsnp138.nochr.vcf \
--contamination_fraction_to_filter 0.02 \
-o cache/mutect_variants_select.vcf \
--output_mode EMIT_VARIANTS_ONLY \
--disable_auto_index_creation_and_locking_when_reading_rods \
-L snv_hg38.bed
