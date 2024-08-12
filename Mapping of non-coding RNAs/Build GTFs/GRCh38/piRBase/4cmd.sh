# the version is not clear, but the files are downloaded in 4.20.2024
wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/hsa.gold.fa.gz
wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/bed/hsa.align.bed.gz
gunzip hsa.align.bed.gz 
gunzip hsa.gold.fa.gz 


python3 filter_pi2gtf.py -o hsa.piRNA.gtf -r hsa.gold.fa -i hsa.align.bed
