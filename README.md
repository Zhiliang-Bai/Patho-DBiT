# Patho-DBiT: Spatially exploring RNA biology in archival formalin-fixed paraffin-embedded (FFPE) tissues

![Figure 1](https://github.com/user-attachments/assets/c136fb07-9dc1-43d0-a77e-4b1b5de72e38)

## 1. Sequence alignment and generation of mRNA expression matrix
Libraries were sequenced on an Illumina NovaSeq 6000 Sequencing System with a paired-end 150bp read length.

Scripts are included in the "Sequence alignment" folder.
### Patho-DBiT raw FASTQ file
Read 1: Contains the cDNA sequence

Read 2: Contains the spatial Barcode A, Barcode B and UMIs
### Reformat FASTQ Read 2 file
The Read 2 sequence needs to be reformated to run [ST Pipeline](https://github.com/SpatialTranscriptomicsResearch/st_pipeline), as explained in the 'Reformat_Read2.md' file. 

To reformat the Read 2, run 'FASTQ_process.sh'.

All raw Read 1 and reformatted Read 2 files can be downloaded from the NCBI Gene Expression Omnibus (GEO) under the accession number [GSE274641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274641) 

### Run ST Pipeline
Run 'ST_Pipeline.sh' with Processed_R2.fastq.gz and Raw_R1.fastq.gz as inputs.

The pipeline requires a spatial barcode index file to decode spatial locations. Two files are provided in the folder: one for Patho-DBiT 50x50 datasets and another for the expanded 100x100 dataset.

The Mouse GRCm38-mm10 or human GRCh38 reference genome was used with STAR v2.7.7a.

### Convert ENSEMBL ID to gene name
Run 'Convert_ENSEMBL_to_gene_name.sh' to annotate the matrix output from the ST pipeline.

## 2. Read mapping of non-coding RNA species
Scripts are included in the "Mapping of non-coding RNAs" folder.
### Build reference genomic data of various non-coding RNA species
Run 'Build_GTF_Pipeline.sh' in the "GRCh38" or "mm39" folder to build the genomic reference for human or mouse.
### Get expression profile of non-coding RNAs
Run 'count_ncRNAs.sh' in the "ncRNA types" folder.

## 3. Spatial alternative splicing analysis
Scripts are included in the "Spatial alternative splicing" folder.
Run 'rMATS_pipeline.sh' with BAM file generated by 'ST_Pipeline.sh'. 

The pipeline also requires the following input files: a GTF annotation file, a spatial barcode index file, a TSV file linking cluster IDs to spatial pixel positions, and a TSV file linking cluster IDs to region names.

## 4. Spatial adenosine-to-inosine (A-to-I) RNA editing analysis
Scripts are included in the "Spatial A-to-I RNA editing" folder.
Run 'RNA_editing_pipeline.sh' with BAM file generated by 'ST_Pipeline.sh'. 

The pipeline also requires the following input files: a reference genome fasta file, a GTF annotation file, a reference RNA editing sites file, a TSV file linking cluster IDs to spatial pixel positions, and a TSV file linking cluster IDs to region names.

## 5. Spatial single nucleotide variant (SNV) analysis
Scripts are included in the "Spatial variant analysis" folder.
### Map WGS data
Run 'BWAmap.sh' in the "pipeline" folder to perform whole-genome sequencing (WGS) data alignment and analysis.
### Generate spatial SNV matrix
Run 'analyzeVariant.sh' in the "pipeline" folder to generate mutation-by-pixel expression matrix.

The pipeline also requires a CSV file containing cluster IDs linked to spatial pixel positions.

## 6. Spatial RNA splicing dynamics analysis
Scripts are included in the "Spatial RNA dynamics" folder.
Run 'RNAdynamics.sh' with spliced and unspliced count matrices as input.

The pipeline also requires two input CSV files: one containing UMAP embeddings and another with cluster IDs, both linked to spatial pixel positions.

## 7. Inferring high-resolution tissue architecture using iStar
Super-resolved tissue architecture was generated by integrating the Patho-DBiT gene expression matrix with high-resolution histology using [iStar](https://github.com/daviddaiweizhang/istar)

## 8. Spatial data visualization
Scripts are included in the "Spatial data visualization" folder.
### Identifying useful pixels in tissue scanning images
Follow the steps listed in 'Image processing.md' for pixel identification, then run 'Pixel_identification.m' to generate a "position.txt" file containing the identified useful pixels from the image.
### Unsupervised clustering analysis and data visualization
Perform spatial unsupervised clustering analysis by executing 'Patho-DBiT_Clustering.Rmd'.



