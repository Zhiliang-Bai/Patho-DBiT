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
tmp="your path"

input=$tmp/${sample}_R2.fastq.gz
output=$tmp/${sample}_R2_processed.fastq

python fastq_process.py --input $input --output $output
