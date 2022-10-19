#!/bin/bash
homeroad=/mnt/data/chenwei/huahua/1.rawdata/AMA_term_vill_RNA/Base_Peri_Bulk_RNAseq
name1=$1
name2=$2
cd $homeroad
mkdir  $name2

mv ${name1}_1.fq.gz  ${name2}_1.fq.gz
mv ${name1}_2.fq.gz  ${name2}_2.fq.gz
mv  ${name2}_1.fq.gz ${name2}_2.fq.gz  $name2


