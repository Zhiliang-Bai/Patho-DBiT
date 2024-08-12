#!/bin/bash
#SBATCH -p pi_gerstein -A gerstein
#SBATCH --job-name=DNAmutation
#SBATCH -o cache/DNAmutation.%j.out
#SBATCH --ntasks=1 --cpus-per-task=16
#SBATCH --mem-per-cpu=80G
#SBATCH --time=1-20:00:00


module load BWA GATK SAMtools
{
  inputfq1=data/DNAseq/LM0623/LM0623_1.fq.gz
  inputfq2=data/DNAseq/LM0623/LM0623_2.fq.gz
  outputfolder=result/LM0623_map/DNA/Tumor/
  groupname='@RG\tID:LM0623\tSM:LM0623\tPL:illumina\tLB:LM0623'
  
  indexfq=/home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa
  known_indels=data/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf
  
  mappedbam=$outputfolder/'mapped.bam'
  sortedbam=$outputfolder/'sorted.bam'
  markedbam=$outputfolder/'marked.bam'
  
  
  bwa mem -t 16 -T 0 -R $groupname $indexfq $inputfq1 $inputfq2 | samtools view -Shb -o $mappedbam -
  samtools sort -o $sortedbam $mappedbam
  samtools index $sortedbam
  
  
  gatk MarkDuplicates CREATE_INDEX=true INPUT=$sortedbam O=$markedbam M=${markedbam}.txt VALIDATION_STRINGENCY=STRICT
}


{
  inputfq1=data/DNAseq/Normal/Normal_1.fq.gz
  inputfq2=data/DNAseq/Normal/Normal_2.fq.gz
  outputfolder=result/LM0623_map/DNA/Normal/
  groupname='@RG\tID:Normal\tSM:Normal\tPL:illumina\tLB:Normal'
  
  indexfq=/home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa
  known_indels=data/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf
  
  mappedbam=$outputfolder/'mapped.bam'
  sortedbam=$outputfolder/'sorted.bam'
  markedbam=$outputfolder/'marked.bam'
  
  bwa mem -t 16 -T 0 -R $groupname $indexfq $inputfq1 $inputfq2 | samtools view -Shb -o $mappedbam -
  
  samtools sort -o $sortedbam $mappedbam
  samtools index $sortedbam
  gatk MarkDuplicates CREATE_INDEX=true INPUT=$sortedbam O=$markedbam M=${markedbam}.txt VALIDATION_STRINGENCY=STRICT
}

samtools view -@ 16 -F 1024 -b result/LM0623_map/DNA/Normal/marked.bam -o result/LM0623_map/DNA/Normal/dedup_nor.bam
samtools view -@ 16 -F 1024 -b result/LM0623_map/DNA/Tumor/marked.bam  -o result/LM0623_map/DNA/Tumor/dedup_tmr.bam




