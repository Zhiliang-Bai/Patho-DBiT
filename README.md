# Patho-DBiT: Spatially exploring RNA biology in archival formalin-fixed paraffin-embedded (FFPE) tissues

![Figure 1](https://github.com/user-attachments/assets/c136fb07-9dc1-43d0-a77e-4b1b5de72e38)

## 1. Sequence alignment and generation of mRNA expression matrix
Libraries were sequenced on an Illumina NovaSeq 6000 Sequencing System with a paired-end 150bp read length.

Scripts are included in the "Sequence alignment" folder.
### Patho-DBiT raw FASTQ file
Read 1: Contains the cDNA sequence

Read 2: Contains the spatial Barcode A, Barcode B and UMIs
### Reformat FASTQ Read 2 file
The Read 2 sequence needs to be reformated to run [ST Pipeline](https://github.com/SpatialTranscriptomicsResearch/st_pipeline), as explained in the following figure. Due to experimental design of Patho-DBiT, the Read 2 is equal to the "Read 1" in ST pipeline, while Read 1 will be the "Read 2".

To reformat the Read 2, run 'FASTQ_process.sh'.

All raw Read 1 and reformatted Read 2 files can be downloaded from the NCBI Gene Expression Omnibus (GEO) under the accession number XXXXXXXX 

### Run ST Pipeline
Run 'ST_Pipeline.sh' with Processed_R2.fastq.gz and Raw_R1.fastq.gz as inputs.

The pipeline requires a spatial barcode index file to decode spatial locations. Two files are provided in the folder: one for Patho-DBiT 50x50 datasets and another for the expanded 100x100 dataset.

The Mouse GRCm38-mm10 or human GRCh38 reference genome was used with STAR v2.7.7a.

### Convert ENSEMBL ID to gene name
Run 'Convert_ENSEMBL_to_gene_name.sh' to annotate the matrix output from the ST pipeline.

## 2. Read mapping of non-coding RNA species


## 3. Spatial alternative splicing analysis

## 4. Spatial adenosine-to-inosine (A-to-I) RNA editing analysis

