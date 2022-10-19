#!/bin/bash
homeroad=$1
data=$2

#chrLambda=/public1/chenwei/seq_index/lambda_DNA_r.fa
#out_dir=/public1/chenwei/methylation/$resultfile/
fastqc=/mnt/data/chenwei/software/FastQC/fastqc
#cutadapt=/data1/yuanpeng/software/cutadapt-1.9.1/bin/cutadapt
#trim_g=/media/data1/yanzhiqiang/software/trim_galore_zip/trim_galore
#qualimap=/public1/chenwei/software/qualimap_v2.2/qualimap
mkdir -p $homeroad/$data/rawQC

for id in  $homeroad/$data/*.gz
do
echo `date` "start do QC for " $id
echo $id
$fastqc -q $id -o $homeroad/$data/rawQC
echo `date` "end do QC for " $id
done

#remove adapter & discard poor quality reads
output_road=/mnt/data/chenwei/huahua/3.mapped_data/project_two_RNA
mkdir -p $output_road/$data
trim_galore    --quality 20 --phred33 --stringency 3 --gzip --length 36 --paired --fastqc  -output_dir $output_road/$data   $homeroad/$data/${data}_1.fq.gz $homeroad/$data/${data}_2.fq.gz

#mkdir $homeroad/QC_results
#mv *zip *html $homeroad/QC_results/

