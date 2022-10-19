#!/bin/bash
homeroad=$1
#rawdata=$2
#rawdata2=$3
samtools=/mnt/data/chenwei/software/miniconda2_new/bin/samtools
#for rawdata  in  $homeroad/*
#for rawdata  in `ls  $homeroad/`
for rawdata  in `ls $homeroad`
do
$samtools sort  --threads 10   $homeroad/$rawdata/${rawdata}_1_val_1_bismark_bt2_pe.bam -o $homeroad/$rawdata/${rawdata}_1_val_1_bismark_bt2_pe_sort.bam
$samtools index $homeroad/$rawdata/${rawdata}_1_val_1_bismark_bt2_pe_sort.bam
#rm  $homeroad/$rawdata/${rawdata}_1_val_1_bismark_bt2_pe.bam 
done
