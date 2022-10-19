#!/bin/bash
#homeroad=/mnt/data/chenwei/huahua/1.rawdata/AMA_term_vill_meth
homeroad=/mnt/data/chenwei/huahua/1.rawdata/YZ_Rawdata
name1=$1
name2=$2
cd $homeroad
mv $name1 $name2
cd $name2
#mv ${name1}_2.clean.fq.gz  ${name2}_2.fq.gz
#mv ${name1}_1.clean.fq.gz  ${name2}_1.fq.gz

#mv ${name1}.R2.fastq.gz  ${name2}_2.fq.gz
#mv ${name1}.R1.fastq.gz  ${name2}_1.fq.gz

mv ${name1}_R2.fq.gz  ${name2}_2.fq.gz
mv ${name1}_R1.fq.gz  ${name2}_1.fq.gz

#mv ${name1}_2.fq.gz  ${name2}_2.fq.gz
#mv ${name1}_1.fq.gz  ${name2}_1.fq.gz


