#!/bin/bash
#list="11-2-3E  12-1-3G  12-3-4A  13-2-1H  15-1-2A  15-3-2C  16-2-4C  17-1-2D  17-3-2F  19-2-2H  1F  1Q  20M-4F  2-1-3A  2-3-3C  3M  4F-5B  4Q-5C  5-2  6F  6Q-5D  7M  8-1  8-3  9-2-1B  10-2-1D  11-1-3D  11-3-3F  12-2-3H  13-1-1F  13-3-1G  15-2-2B  16-1-4B  16-3-4D  17-2-2E  19-1-2G  19-3-4E  1M  20F-4G  20Q-4H  2-2-3B  3F  3Q  4M-5A  5-1 5-3  6M  7F 7Q-5E  8-2  9-1-1A  9-3-1C"

list="10-1-in-9"
for sample in $list;
 do
  mkdir /mnt/data/chenwei/huahua/1.rawdata/merge_fastq/$sample
  cat /mnt/data/chenwei/huahua/1.rawdata/first_seq/Rawdata/${sample}/${sample}_R1.fq.gz  /mnt/data/chenwei/huahua/1.rawdata/huahua_add/${sample}/${sample}_R1.fq.gz  > /mnt/data/chenwei/huahua/1.rawdata/merge_fastq/${sample}/${sample}_1.fq.gz
  cat /mnt/data/chenwei/huahua/1.rawdata/first_seq/Rawdata/${sample}/${sample}_R2.fq.gz  /mnt/data/chenwei/huahua/1.rawdata/huahua_add/${sample}/${sample}_R2.fq.gz  > /mnt/data/chenwei/huahua/1.rawdata/merge_fastq/${sample}/${sample}_2.fq.gz
  rm -rf  /mnt/data/chenwei/huahua/1.rawdata/first_seq/Rawdata/$sample  /mnt/data/chenwei/huahua/1.rawdata/huahua_add/$sample
 done

