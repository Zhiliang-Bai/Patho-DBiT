#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=VCF2mut
#SBATCH -o cache/VCF2mut.%j.out
#SBATCH --ntasks=1 --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=23:00:00
#SBATCH --mail-type=END

module load SAMtools miniconda parallel BCFtools
conda activate VCF


multi_pileup () {
    BAM_FILE=$1
    candidatelocationlist=$2
    referencegenome=$3
    outputfolder=$4
    thename=`basename $BAM_FILE`
    thename=${thename%.bam}
    mkdir -p $outputfolder
    #echo "samtools mpileup -d 100000000 -l $candidatelocationlist -f $referencegenome -o ${outputfolder}'/'${thename}.vcf $BAM_FILE"
    samtools mpileup -d 100000000 -l $candidatelocationlist -f $referencegenome -o ${outputfolder}'/'${thename}.vcf $BAM_FILE
}
export -f multi_pileup





# 初始化变量
Mode=""
Output=""
Bams=""
VCFin=""
FastaFile="/home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa"

# 参数列表
params=("$@")
num_params=${#params[@]}

script="${params[-1]}"
mkdir -p $2
logf=$2/log.txt
echo 'parameter' > $logf


for ((i=0; i<num_params-1; i++)); do
  case $i in
    0) Mode="${params[$i]}" ;;
    1) Output="${params[$i]}" ;;
    2) Bams="${params[$i]}" ;;
    3) VCFin="${params[$i]}" ;;
    4) FastaFile="${params[$i]}" ;;
  esac
  echo ${params[$i]} >> $logf
done


if [[ $Mode == *"a"* ]]; then
    pixelfolder=$Output/pixels/
    mkdir $pixelfolder
    perl ${script}/split_fa_pixel.pl -i $Bams -o $pixelfolder
else
    pixelfolder=$Bams
fi

if [[ $Mode == *"b"* ]]; then
    bcftools view -f PASS $VCFin | egrep -v '^#' > $Output/PASS.vcf
    awk '{print $1 "\t" $2}' $Output/PASS.vcf | sort | uniq > $Output/PASS.site
    awk '{print $1 ":" $2 ":" $4 ":" $5}' $Output/PASS.vcf | sort | uniq > $Output/PASS.site.pos.txt
    interestSite=$Output/PASS.site
    posinterestSite=$Output/PASS.site.pos.txt
else
    posinterestSite=$VCFin
    interestSite=$Output/PASS.site
    cut -f 1-2 -d : $posinterestSite  | sed 's/:/\t/' > $interestSite
fi


if [[ $Mode == *"c"* ]]; then
    Pilefolder=$Output/pixelPile/
    bamlist=$Output/pixelBam.list
    if [ -d "$pixelfolder" ]; then
        find $pixelfolder -name '*.bam' > $bamlist
    else
        bamlist=$pixelfolder
    fi
    echo "cat $bamlist | parallel -j $SLURM_CPUS_PER_TASK multi_pileup {} $interestSite $FastaFile $Pilefolder"
    cat $bamlist | parallel -j $SLURM_CPUS_PER_TASK multi_pileup {} $interestSite $FastaFile $Pilefolder
else
    Pilefolder=$FastaFile
fi

if [[ $Mode == *"d"* ]]; then
    echo "Rscript ${script}/VCF2mut.R $Pilefolder $Output/mut.mat $posinterestSite $SLURM_CPUS_PER_TASK"
    Rscript ${script}/VCF2mut.R $Pilefolder $Output/mut.mat $posinterestSite $SLURM_CPUS_PER_TASK
    Rscript ${script}/VCF2mut_ratio.R $Pilefolder $Output/mutratio.mat $posinterestSite $SLURM_CPUS_PER_TASK
    Rscript ${script}/VCF2mut_ratiostr.R $Pilefolder $Output/mutratiostr.mat $posinterestSite $SLURM_CPUS_PER_TASK
fi


PSGmpileup () {
    candidatelocationlist=$1
    referencegenome=$2
    bamlistfile=$3
    outputfile=$4
    tmpfilefolder=${5:-tmptmptmp/}
    tmpfilefolder="${tmpfilefolder%/}"
    tmpfilefolder=$tmpfilefolder'/'
    threadnum=${6:-16}
    
    bamfileone=`head $bamlistfile -n 1`
    rm $tmpfilefolder/tmp.*.vcf
    mkdir $tmpfilefolder/
    
    readarray -t chr_names < <(samtools view -H "$bamfileone" | grep "@SQ" | sed 's/^.*SN://g' | cut -f 1)
    printf "%s\n" "${chr_names[@]}" | xargs -I {} -n 1 -P $threadnum sh -c "samtools mpileup -d 100000000 -l $candidatelocationlist -f $referencegenome -r {} -b $bamlistfile -o ${tmpfilefolder}tmp.{}.vcf"
    cat ${tmpfilefolder}tmp.*.vcf > $outputfile
    rm ${tmpfilefolder}tmp.*.vcf
    rmdir $tmpfilefolder
}
#PSGmpileup cache/DNARNA_germline.txt /home/dz287/gibbs/scRNAmir/data/Genomes/GRCh38/Homo_sapiens.GRCh38.dna.fa cache/bamlist.txt result/comparegermline.txt tmptmp/ 16
