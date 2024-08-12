# the version is not clear, but the files are downloaded in 4.20.2024
wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/mmu.gold.fa.gz
wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/bed/mmu.align.bed.gz
gunzip mmu.align.bed.gz 
gunzip mmu.gold.fa.gz 
module load BEDTools # make sure you could use bedtools
python3 filter_pi2gtf.py -o mmu.piRNA.bed -r mmu.gold.fa -i mmu.align.bed -b 
# need to change genome version for mmu.piRNA.bed
# use UCSC Genome browser (liftover  set "Minimum ratio of bases that must remap" as 1) change the genome version form GRCm38 to GRCm39/mm39  and the file mmu.piRNA_mm39.bed is gotten
bedtools getfasta -fi mm39.fa  -fo mmu.piRNA_mm39.fa -bed mmu.piRNA_mm39.bed -name -s
python3 examineFAs.py -i mmu.piRNA_mm39.fa -r mmu.gold.fa -o mmu.piRNA_mm39.gtf  #output final gtf and report errors
#examineFAs.py examined results and showed that no piRNAs should be deleted in this new gtf

