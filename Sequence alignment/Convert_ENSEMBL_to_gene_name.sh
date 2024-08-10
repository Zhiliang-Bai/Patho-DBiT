#!/bin/bash
#SBATCH --partition='your partition name'
#SBATCH --job-name='your name'
#SBATCH --mail-type=END
#SBATCH -o test.%j.out
#SBATCH -e test.%j.err
#SBATCH --mail-user='your email'
#SBATCH --ntasks=1 --cpus-per-task=5
#SBATCH --mem=10g 
#SBATCH --time=1:00:00

module load miniconda

conda activate 'your conda environment'

sample=Sample_name
OUTPUT=/Path to matrix output from the ST pipeline/

tsv_E=$OUTPUT/${sample}_stdata.tsv
path_to_annotation_file=/Path to reference genome annotation.gtf

convertEnsemblToNames.py --annotation $path_to_annotation_file --output $OUTPUT/${sample}_stdata_names.tsv $tsv_E