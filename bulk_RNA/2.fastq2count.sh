#!/bin/bash
##https://blog.csdn.net/u012110870/article/details/102804525
#mapping
STAR=/mnt/data/chenwei/software/anaconda3/bin/STAR
qualimap=/mnt/data/chenwei/software/anaconda3/bin/qualimap

genemome=/mnt/data/yanzhiqiang/138/DB/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
gtf=/mnt/data/yanzhiqiang/138/DB/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf

featureCounts=/mnt/data/chenwei/software/anaconda3/bin/featureCounts

##make index
#$STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /mnt/data/chenwei/software/STARIndex --genomeFastaFiles $genemome --sjdbGTFfile $gtf --sjdbOverhang 99


#do the alignment
#list="A21"
##fastq_road=/mnt/data/chenwei/huahua/2.cleandata/project_two_RNA
##Road=/mnt/data/chenwei/huahua/3.mapped_data/project_two_RNA
##list="A_B1  A_B2  A_B3  A_B4  A_B5  A_P1  A_P2  A_P3  A_P4  A_P5  Y_B1  Y_B2  Y_B3  Y_B4  Y_B5  Y_P1  Y_P2  Y_P3  Y_P4  Y_P5"

#fastq_road=/mnt/data/chenwei/huahua/2.cleandata/AMA_term_vill_RNA
#Road=/mnt/data/chenwei/huahua/3.mapped_data/project_two_RNA/term_vill_RNAdata

fastq_road=/mnt/data/chenwei/huahua/2.cleandata/YZ_RNAdata
Road=/mnt/data/chenwei/huahua/3.mapped_data/YZ_data

#list="AMA_1  AMA_10  AMA_11  AMA_2  AMA_3  AMA_4  AMA_5  AMA_6  AMA_7  AMA_8  AMA_9  YMA_1  YMA_10  YMA_2  YMA_3  YMA_4  YMA_5  YMA_6  YMA_7  YMA_8  YMA_9"

#list="CTRL_1  CTRL_2  OE_1  OE_2"
list="CTRL_2  OE_1  OE_2"
#list="A3M  A3Q  A4M  A4Q  A5M  A5Q  A6M  A6Q  A7M  A7Q  Y1M  Y1Q  Y2M  Y2Q  Y4M  Y4Q  Y7M  Y7Q  Y8M  Y8Q"
#list="A2M"
for sample in $list;
do
 $STAR --runThreadN 20 --genomeDir  /mnt/data/chenwei/software/STARIndex  --readFilesIn $fastq_road/$sample/${sample}_1_val_1.fq.gz  $fastq_road/$sample/${sample}_2_val_2.fq.gz  --outFileNamePrefix $sample --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10
 mv ${sample}Aligned.sortedByCoord.out.bam ${sample}.bam
done

#Trim reads map to multiple regions. #提取比对到参考基因组上的数据 ##https://www.jianshu.com/p/667ad2230a8c
#mkdir $Road/3.mapdata
mv ./*.bam $Road

for i in $list;
do
$qualimap --java-mem-size=20G  bamqc -bam $Road/${i}.bam  -outfile   ${i}.mapreport.pdf
#$qualimap  --java-mem-size=20G  bamqc -bam $Road/3.mapdata/${i}Aligned.sortedByCoord.out.bam  -outfile   ${i}_aftrmdup.mapreport.pdf
done
mkdir $Road/qualimap_stat
mv ./*_stats $Road/qualimap_stat
#htseq-count -r name -f bam -s no -m union ${bam}.name.sort.bam $gtf > 2.htseq/${sample}.count 
#htseq-count -f bam -s no Hs_treat3_sort.bam ~/Ref/UCSC_hg19/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > sample.count

#随后用featurecounts计数 https://www.bioinfo-scrounger.com/archives/407/

mkdir  $Road/4.count
for  i in $list;
do 
$featureCounts -T 10 -p -t exon -g gene_id  -a $gtf -o $Road/4.count/${i}.counts.txt  $Road/${i}.bam
#$featureCounts -T 10 -p -t exon -g gene_id  -a $gtf -o $Road/4.count/${i}.counts.txt  $Road/3.mapdata/${i}Aligned.sortedByCoord.out.bam 
done

##

for  i in $list;
do
 cut -f 1,7 $Road/4.count/${i}.counts.txt |grep -v '^#' >  $Road/4.count/${i}_feactureCounts.txt
 sed -i '1d'  $Road/4.count/${i}_feactureCounts.txt
 perl -alne '$sum += $F[1]; END {print $sum}'  $Road/4.count/${i}_feactureCounts.txt
done

mkdir $Road/4.count/feactureCounts 
mv $Road/4.count/*_feactureCounts.txt $Road/4.count/feactureCounts
mkdir $Road/4.count/raw_Counts_files
mv $Road/4.count/*counts.txt*  $Road/4.count/raw_Counts_files

