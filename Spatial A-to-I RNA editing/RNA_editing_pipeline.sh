#!/bin/bash
# scripts to run RNA editing pipeline
gatk="Path to gatk" # https://gatk.broadinstitute.org/hc/en-us
samtools="Path to samtools" # http://www.htslib.org/
conda install bioconda::cgranges # https://github.com/lh3/cgranges https://anaconda.org/bioconda/cgranges
collect_editing_site_py=collect_editing_sites.py

# input bam
in_bam="Path to input spatial RNA-seq bam" # input bam

# spatial information
pixel_to_clu_ID_tsv="/Path to spatial barcode files/pixel_to_clu_ID.tsv" # eg: 1x50 2
                                                                         #     10x19 1
                                                                         #     11x18 7
clu_ID_to_region_names="/Path to spatial barcode files/clu_ID_to_reg_name.tsv"  # eg: 9   Midbrain
                                                                                #     12  Midbrain
                                                                                #     4   Isocortex

# reference genome
ref_fa="Path to reference genome fasta file"
gtf="Path to GTF annotation file"

# reference editing sites
ref_edit_site_tsv="Path to reference editing sites tsv file" # http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_mm10.txt.gz

# 1. generate dict file for reference genome
ref_dict=${ref_fa}.dict

$gatk CreateSequenceDictionary -R $ref_fa -O $ref_dict

# 2. split bam based on N cigars
split_bam="Path to split bam file"

$gatk SplitNCigarReads -R $ref_fa -I $in_bam -O $split_bam

# 3. mpileup
min_map_qual=25
min_base_qual=25
mpileup_flags="--no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends -B -d 0 -Q $min_base_qual -q $min_map_qual"
editing_out_tsv="Path to output tsv file"

$samtools mpileup $mpileup_flags -l $ref_edit_site_tsv -f $ref_fa $split_bam --output-QNAME -o $editing_out_tsv

# 4. collect editing sites
min_read_depth=10
min_edit_read_count=1
edit_out_dir="Path to output dir" # output dir for editing sites

python $collect_editing_site_py -m $min_read_depth -c $min_edit_read_count $editing_out_tsv $pixel_to_clu_ID_tsv $clu_ID_to_region_names $gtf $edit_out_dir