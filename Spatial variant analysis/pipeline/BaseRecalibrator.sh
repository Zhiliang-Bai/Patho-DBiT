#!/bin/bash
#SBATCH -p pi_gerstein -A gerstein
#SBATCH --job-name=BaseRecalibrator
#SBATCH -o cache/BaseRecalibrator.%j.out
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --time=6-23:00:00


module load BWA GATK/3.8-1-0-gf15c1c3ef-Java-1.8 SAMtools
inputfq1=data/Normal/Normal_1.fq.gz
inputfq2=data/Normal/Normal_2.fq.gz
outputfolder=result/Normal/
groupname='@RG\tID:Normal\tSM:Normal\tPL:illumina\tLB:Normal'


indexfq=/home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa
known_indels=data/Homo_sapiens_assembly38.known_indels.vcf.gz
dbsnp_vcf=data/Homo_sapiens_assembly38.dbsnp138.nochr.vcf

markedbam1=result/Normal/marked.bam
markedbam2=result/LM0623/marked.bam
realign_target_intervals='realign_target.intervals'
bqsr_grp='bqsr_grp.txt'
aligned1=Normal_marked_Realigned.bam
aligned2=LM0623_marked_Realigned.bam
refinedbam1='refined_normal.bam'
refinedbam2='refined_tumor.bam'

#java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#-T RealignerTargetCreator \
#-R $indexfq \
#-known $known_indels \
#-I $markedbam1  \
#-I $markedbam2  \
#-o $realign_target_intervals 

#markedbam01=result/Normal/Normal_marked.bam
#markedbam02=result/LM0623/LM0623_marked.bam

#mv $markedbam1 $markedbam01
#mv $markedbam2 $markedbam02



#java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#-T IndelRealigner \
#-R $indexfq \
#-known $known_indels \
#-targetIntervals $realign_target_intervals \
#--noOriginalAlignmentTags \
#-I $markedbam01  \
#-I $markedbam02  \
#-nWayOut _Realigned.bam

echo '33'
java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $indexfq \
-I $markedbam1 \
-I $markedbam2 \
-knownSites $dbsnp_vcf \
-o $bqsr_grp


echo '441'
java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T PrintReads \
-R $indexfq \
-I $markedbam1 \
--BQSR $bqsr_grp \
-o $refinedbam1


echo '442'
java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T PrintReads \
-R $indexfq \
-I $markedbam2 \
--BQSR $bqsr_grp \
-o $refinedbam2