#!/bin/bash
# script to run rMATS-turbo for spatial RNA-seq data
rmats="Path to rmats-turbo/rmats.py" # https://github.com/Xinglab/rmats-turbo
split_py=split_bam_by_clu.py

# input bam
in_bam="Path to input Patho-DBiT bam" # input bam
read_len=120 # read length

# spatial information
barcode_to_pixel_coor="/Path to spatial barcode files/Spatial_barcode_50x50.txt"
pixel_to_clu_ID_tsv="/Path to spatial barcode files/pixel_to_clu_ID.tsv" # eg: 1x50 2
                                                                         #     10x19 1
                                                                         #     11x18 7
clu_ID_to_region_names="/Path to spatial barcode files/clu_ID_to_reg_name.tsv"  # eg: 9   Midbrain
                                                                                #     12  Midbrain
                                                                                #     4   Isocortex
N_regions=7


# 1. split bam by region
bam_by_clu_dir="bam_by_clu_name" # output dir of region-wise bam files
python $split_py $in_bam $barcode_to_pixel_coor $pixel_to_clu_ID_tsv $clu_ID_to_region_names $bam_by_clu_dir

# 2. pairwise rMATS-turbo
threads=1 # number of threads to use
gtf="Path to GTF annotation file"
all_clu_bams=($(ls $bam_by_clu_dir/*.bam))
rmats_out_dir="rmats_out_dir" # root dir for rMATS output
for((clu_i=0;clu_i<$N_regions-1;++clu_i))
do
    for((clu_j=clu_i+1;clu_j<$N_regions;++clu_j))
    do
        clu_i_bam=${all_clu_bams[$clu_i]}
        clu_j_bam=${all_clu_bams[$clu_j]}
        basename_i=$(basename ${clu_i_bam%.bam})
        basename_j=$(basename ${clu_j_bam%.bam})
        od=${rmats_out_dir}/${basename_i}_${basename_j}
        mkdir -p $od
        echo "$clu_i_bam" > $od/b1.txt
        echo "$clu_j_bam" > $od/b2.txt
        python $rmats --gtf $gtf --b1 $od/b1.txt --b2 $od/b2.txt --od $od --tmp $od/tmp --nthread $threads -t single --readLength $read_len --allow-clipping --variable-read-length
    done
done
