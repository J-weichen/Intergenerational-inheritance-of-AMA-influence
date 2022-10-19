#!/bin/bash
sampledir=$1

#samtools=/media/data1/yanzhiqiang/software/samtools-1.5 
chrom_bed=/mnt/data/chenwei/qinmen_BR/0.meng_script/1.meth_script/CNV_1Mb_script/2.1Mb.24chrom.bed
#A25_1_val_1_bismark_bt2_pe_sort.bam
for bam in $sampledir/*/*_1_val_1_bismark_bt2_pe_sort.bam
#for bam in $sampledir/F/F.bam $sampledir/M/M.bam $sampledir/C_test/C_test.unsplit.bam $sampledir/D_test/D_test.unsplit.bam 
do
	temp=`dirname $bam`
	sample=`basename $temp`
	echo $sample
	samtools bedcov $chrom_bed  $bam > ${temp}/${sample}.1M.cov 
        #mv ${temp}/${sample}.1M.cov /media/data2/ermeng/200119_5G/1.rawdata/CNV_result/${sample}
#        nohup $samtools bedcov $chrom_bed  $bam > ${temp}/${sample}.10M.cov 2> ${sample}.log &
done



