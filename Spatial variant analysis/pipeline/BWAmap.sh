inputfq1=Sample1_Normal/r1.fq.gz
inputfq2=Sample1_Normal/r2.fq.gz
outputfolder=result/Sample1_Normal/
groupname='@RG\tID:LM0623\tSM:LM0623\tPL:illumina\tLB:Sample1_Normal'
indexfq=Homo_sapiens.GRCh38.dna.fa

# make sure BWA GATK SAMtools are all available
{
  mappedbam=$outputfolder/'mapped.bam'
  sortedbam=$outputfolder/'sorted.bam'
  markedbam=$outputfolder/'marked.bam'
  
  
  bwa mem -t 16 -T 0 -R $groupname $indexfq $inputfq1 $inputfq2 | samtools view -Shb -o $mappedbam -
  samtools sort -o $sortedbam $mappedbam
  samtools index $sortedbam
  
  
  gatk MarkDuplicates CREATE_INDEX=true INPUT=$sortedbam O=$markedbam M=${markedbam}.txt VALIDATION_STRINGENCY=STRICT
}

inputfq1=Sample2_Tumor/Normal_1.fq.gz
inputfq2=Sample2_Tumor/Normal_2.fq.gz
outputfolder=result/Sample2_Tumor/
groupname='@RG\tID:Normal\tSM:Normal\tPL:illumina\tLB:Sample2_Tumor'
  

{


  
  mappedbam=$outputfolder/'mapped.bam'
  sortedbam=$outputfolder/'sorted.bam'
  markedbam=$outputfolder/'marked.bam'
  
  bwa mem -t 16 -T 0 -R $groupname $indexfq $inputfq1 $inputfq2 | samtools view -Shb -o $mappedbam -
  
  samtools sort -o $sortedbam $mappedbam
  samtools index $sortedbam
  gatk MarkDuplicates CREATE_INDEX=true INPUT=$sortedbam O=$markedbam M=${markedbam}.txt VALIDATION_STRINGENCY=STRICT
}

samtools view -@ 16 -F 1024 -b result/Sample2_Tumor/marked.bam -o dedup_nor.bam
samtools view -@ 16 -F 1024 -b result/Sample2_Tumor/marked.bam  -o dedup_tmr.bam




