#!/bin/bash
homeroad=$1
rawdata1=$2
rawdata2=$3


#Version=2
Bisulfite_Genome_lambda=/mnt/data/yanzhiqiang/138/DB/lambda/bismarkIndex
chrLambda=/mnt/data/yanzhiqiang/138/DB/lambda/bismarkIndex/lambda.fa
fastqc=/mnt/data/chenwei/software/FastQC/fastqc
#trim_galore=/mnt/data/chenwei/software/TrimGalore-0.6.6/trim_galore
bismark=/mnt/data/yanzhiqiang/software/anaconda3/bin/bismark
samtools=/mnt/data/chenwei/software/miniconda2_new/bin/samtools
## step1 
for id in  $homeroad/*.gz
do
echo $id
$fastqc -q $id
done


##step2
#remove adapter & discard poor quality reads
trim_galore  --quality 20 --phred33 --stringency 3 --gzip --length 36 --rrbs --paired  --trim1 -output_dir $homeroad/  $homeroad/${rawdata1}.fq.gz $homeroad/${rawdata2}.fq.gz 

##step 3

$fastqc -q $homeroad/${rawdata1}_val_1.fq.gz
$fastqc -q  $homeroad/${rawdata2}_val_2.fq.gz

#prepare Genome
#bismark_genome_preparation --path_to_bowtie /public1/test/software/bowtie2/ --bowtie2 --verbose /public1/chenwei/seq_index/Bowtie_2_Bisulfite_Genome/

#

$bismark   -path_to_bowtie /mnt/data/yanzhiqiang/software/anaconda3/bin  --bowtie2  --samtools_path /mnt/data/chenwei/software/miniconda2_new/bin   --parallel 4  -o $homeroad --temp_dir  $homeroad/temp  $Bisulfite_Genome_lambda  -1 $homeroad/${rawdata1}_val_1.fq.gz   -2   $homeroad/${rawdata2}_val_2.fq.gz \
&& java -jar /mnt/data/chenwei/software/picard-tools-1.141/picard.jar SortSam  I=$homeroad/${rawdata1}_val_1_bismark_bt2_pe.bam   O=$homeroad/${rawdata1}_Lambda.sort.bam   SORT_ORDER=coordinate   VALIDATION_STRINGENCY=LENIENT   VERBOSITY=ERROR \
&& $samtools index $homeroad/${rawdata1}_Lambda.sort.bam \
&& $samtools view -h $homeroad/${rawdata1}_Lambda.sort.bam > $homeroad/${rawdata1}_chrLambda.sam \
&& perl /mnt/data/chenwei/qinmen_BR/0.meng_script/1.meth_script/MethylExtractBSCR.pl seqFile=$chrLambda  inFile=$homeroad/${rawdata1}_chrLambda.sam  flagW=99,147  flagC=83,163  >  $homeroad/${rawdata1}_Bisulfy_Conversion_rate_report

rm $homeroad/*.bam   $homeroad/*.sam $homeroad/*.bai 
#
