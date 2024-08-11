#!/bin/bash
#SBATCH --partition='your partition name'
#SBATCH --job-name='your name'
#SBATCH --mail-type=END
#SBATCH -o test.%j.out
#SBATCH -e test.%j.err
#SBATCH --mail-user='your email'
#SBATCH --ntasks=1 --cpus-per-task=5
#SBATCH --mem=100g 
#SBATCH --time=10:00:00

module load miniconda

conda activate 'your conda environment'

# Do not add / or \ to the sample name
sample=Sample_name
tmp="your path"

# FASTQ reads
FW=$tmp/${sample}_R2_processed.fastq
RV=$tmp/${sample}_R1_filtered.final.fastq.gz

# References for mapping, annotation and nonRNA-filtering
MAP="Path to reference genome file"
ANN="Path to reference genome annotation.gtf"


# Barcodes settings
ID=/Path to spatial barcode files/Spatial_barcode_50x50.txt

# Output folder and experiment name
OUTPUT=/Path to output folder/
mkdir -p $OUTPUT

TMP_ST=$OUTPUT/tmp
mkdir -p $TMP_ST

# Running the pipeline
st_pipeline_run.py \
  --output-folder $OUTPUT \
  --ids $ID \
  --ref-map $MAP \
  --ref-annotation $ANN \
  --expName $sample \
  --htseq-no-ambiguous \
  --verbose \
  --log-file $OUTPUT/${sample}_log.txt \
  --allowed-kmer 5 \
  --mapping-threads 20 \
  --temp-folder $TMP_ST \
  --no-clean-up \
  --umi-start-position 16 \
  --umi-end-position 26 \
  --overhang 0 \
  --min-length-qual-trimming 18 \
  $FW $RV