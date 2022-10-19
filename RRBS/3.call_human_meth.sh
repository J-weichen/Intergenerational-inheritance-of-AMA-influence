#!/bin/bash
homeroad=$1
rawdata1=$2
rawdata2=$3
sample_name=$4


output_road=/mnt/data/chenwei/huahua/3.mapped_data
#Bisulfite_Genome=/media/data2/chenwei_RRBS_project/meth_DB/bismark2_Bowtie_2_Bisulfite_Genome
Bisulfite_Genome=/mnt/data/yanzhiqiang/138/DB/Homo_sapiens/UCSC/hg38/Sequence/bismarkIndex
bismark=/mnt/data/yanzhiqiang/software/anaconda3/bin/bismark
samtools=/mnt/data/chenwei/software/miniconda2_new/bin/samtools
#step1
$bismark   -path_to_bowtie /mnt/data/yanzhiqiang/software/anaconda3/bin  --bowtie2  --samtools_path /mnt/data/chenwei/software/miniconda2_new/bin   --parallel 8  -quiet -o $homeroad --temp_dir  $homeroad/temp $Bisulfite_Genome  -1 $homeroad/${rawdata1}_val_1.fq.gz   -2   $homeroad/${rawdata2}_val_2.fq.gz
#chrsome_list="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chrX chrY"
#$samtools view $homeroad/${rawdata1}_val_1_bismark_bt2_pe.bam chrsome_list > $homeroad/${sample_name}_24.bam

#step2
/mnt/data/yanzhiqiang/software/anaconda3/bin/bismark_methylation_extractor --parallel 8 --comprehensive --cytosine_report   --paired-end  --no_overlap  --buffer_size 30G  --genome_folder  $Bisulfite_Genome/  $homeroad/${rawdata1}_val_1_bismark_bt2_pe.bam  -o $homeroad
#/mnt/data/yanzhiqiang/software/anaconda3/bin/bismark_methylation_extractor --parallel 4 --CX --comprehensive --no_overlap --counts --buffer_size 30G --report --cytosine_report --genome_folder  $Bisulfite_Genome/  $homeroad/${rawdata1}_val_1_bismark_bt2_pe.bam  -o $homeroad
awk  '{OFS="\t";$1=$1;if($6~/CG/) {print $0} }'  $homeroad/${rawdata1}_val_1_bismark_bt2_pe.CpG_report.txt|grep -E -v "chrUn_*|chrEBV|chrLambda|chrM|*_random" - |awk '{OFS="\t"; {print  $1,$2,$2,$4,$5,$4+$5}}' - |awk  '{OFS="\t";$1=$1;if($6>0) {print $1,$2,$3,($4*100)/$6,$4,$5} }' - >  $homeroad/${rawdata1}_CpG.bismark.cov

#awk  '{OFS="\t";$1=$1;if($6~/CHH/) {print $0} }'  $homeroad/${rawdata1}_val_1_bismark_bt2_pe.CX_report.txt|grep -E -v "chrUn_*|chrEBV|chrLambda|chrM|*_random" - |awk '{OFS="\t"; {print  $1,$2,$2,$4,$5,$4+$5}}' - |awk  '{OFS="\t";$1=$1;if($6>0) {print $1,$2,$3,($4*100)/$6,$4,$5} }' - >  $homeroad/${rawdata1}_CHH.bismark.cov
#awk  '{OFS="\t";$1=$1;if($6~/CHG/) {print $0}}'  $homeroad/${rawdata1}_val_1_bismark_bt2_pe.CX_report.txt|grep -E -v "chrUn_*|chrEBV|chrLambda|chrM|*_random" - |awk '{OFS="\t"; {print  $1,$2,$2,$4,$5,$4+$5}}' - |awk  '{OFS="\t";$1=$1;if($6>0) {print $1,$2,$3,($4*100)/$6,$4,$5} }' - >  $homeroad/${rawdata1}_CHG.bismark.cov

awk '{OFS="\t";$1=$1;if($4 == 0) print $1,$2,$3+1,"-100"; else print $1,$2,$3+1,$4}' $homeroad/${rawdata1}_CpG.bismark.cov >  $homeroad/${rawdata1}_CpG.bismark.bedGraph
#awk '{OFS="\t";$1=$1;if($4 == 0) print $1,$2,$3+1,"-100"; else print $1,$2,$3+1,$4}' $homeroad/${rawdata1}_CHG.bismark.cov >  $homeroad/${rawdata1}_CHG.bismark.bedGraph
#awk '{OFS="\t";$1=$1;if($4 == 0) print $1,$2,$3+1,"-100"; else print $1,$2,$3+1,$4}' $homeroad/${rawdata1}_CHH.bismark.cov >  $homeroad/${rawdata1}_CHH.bismark.bedGraph

awk '{OFS="\t";$1=$1;if(($5+$6)>= 10){print $0}}'  $homeroad/${rawdata1}_CpG.bismark.cov |awk '{OFS="\t";$1=$1;if($4 == 0) print $1,$2,$3+1,"-100"; else print $1,$2,$3+1,$4}' - >  $homeroad/${rawdata1}_CpG.bismark_X10.bedGraph
awk '{OFS="\t";$1=$1;if(($5+$6)>= 6){print $0}}'  $homeroad/${rawdata1}_CpG.bismark.cov |awk '{OFS="\t";$1=$1;if($4 == 0) print $1,$2,$3+1,"-100"; else print $1,$2,$3+1,$4}' - >  $homeroad/${rawdata1}_CpG.bismark_X6.bedGraph

#awk '{OFS="\t";$1=$1;if(($5+$6)>= 3){print $0}}'  $homeroad/${rawdata1}_CHG.bismark.cov |awk '{OFS="\t";$1=$1;if($4 == 0) print $1,$2,$3+1,"-100"; else print $1,$2,$3+1,$4}' - >  $homeroad/${rawdata1}_CHG.bismark_X3.bedGraph
#awk '{OFS="\t";$1=$1;if(($5+$6)>= 3){print $0}}'  $homeroad/${rawdata1}_CHH.bismark.cov |awk '{OFS="\t";$1=$1;if($4 == 0) print $1,$2,$3+1,"-100"; else print $1,$2,$3+1,$4}' - > $homeroad/${rawdata1}_CHH.bismark_X3.bedGraph

awk  '{OFS="\t";$1=$1;if($6~/CG/) {print $0} }' $homeroad/${rawdata1}_val_1_bismark_bt2_pe.CpG_report.txt|grep -E -v "chrUn_*|chrEBV|chrLambda|chrM|*_random" - |awk  '{OFS="\t";$1=$1;if(($4+$5)>=1) {print $0} }' - > $homeroad/${rawdata1}_CpG.scale_region.txt
#awk  '{OFS="\t";$1=$1;if($6~/CHG/) {print $0}}' $homeroad/${rawdata1}_val_1_bismark_bt2_pe.CX_report.txt|grep -E -v "chrUn_*|chrEBV|chrLambda|chrM|*_random" - |awk  '{OFS="\t";$1=$1;if(($4+$5)>=1) {print $0} }' - > $homeroad/${rawdata1}_CHG.scale_region.txt
#awk  '{OFS="\t";$1=$1;if($6~/CHH/) {print $0}}' $homeroad/${rawdata1}_val_1_bismark_bt2_pe.CX_report.txt|grep -E -v "chrUn_*|chrEBV|chrLambda|chrM|*_random" - |awk  '{OFS="\t";$1=$1;if(($4+$5)>=1) {print $0} }' - > $homeroad/${rawdata1}_CHH.scale_region.txt

perl /mnt/data/chenwei/huahua/0.hua_script/genebody_around_meth/01.AbsoluteDistance_no_direction_CpG10.pl /mnt/data/chenwei/huahua/0.hua_script/genebody_around_meth/UCSC.hg38.gene.coordination.txt $homeroad/${rawdata1}_CpG.scale_region.txt > $homeroad/${rawdata1}_CpG_X10_meth_level_in_scaleR.txt 
perl /mnt/data/chenwei/huahua/0.hua_script/genebody_around_meth/01.AbsoluteDistance_no_direction_CpG6.pl /mnt/data/chenwei/huahua/0.hua_script/genebody_around_meth/UCSC.hg38.gene.coordination.txt $homeroad/${rawdata1}_CpG.scale_region.txt > $homeroad/${rawdata1}_CpG_X6_meth_level_in_scaleR.txt 

#perl /media/data2/lucunlin/0.lu_meth_script/genebody_around_meth/01.AbsoluteDistance_no_direction_CHG3.pl /media/data2/lucunlin/0.lu_meth_script/genebody_around_meth/UCSC.hg38.gene.coordination.txt $homeroad/${rawdata1}_CHG.scale_region.txt > $homeroad/${rawdata1}_CHG3_meth_level_in_scaleR.txt &
#perl /media/data2/lucunlin/0.lu_meth_script/genebody_around_meth/01.AbsoluteDistance_no_direction_CHH3.pl /media/data2/lucunlin/0.lu_meth_script/genebody_around_meth/UCSC.hg38.gene.coordination.txt $homeroad/${rawdata1}_CHH.scale_region.txt > $homeroad/${rawdata1}_CHH3_meth_level_in_scaleR.txt &



<<block
&&  java -jar /media/data1/yanzhiqiang/software/picard-tools-1.141/picard.jar SortSam  I=$homeroad/${rawdata1}_val_1_bismark_bt2_pe.bam  O=$homeroad/${rawdata1}_val_1_bismark_bt2_pe.sort.bam   SORT_ORDER=coordinate   VALIDATION_STRINGENCY=LENIENT   VERBOSITY=ERROR \
&& samtools index $homeroad/${rawdata1}_val_1_bismark_bt2_pe.sort.bam \
&& samtools view -h $homeroad/${rawdata1}_val_1_bismark_bt2_pe.sort.bam  chrLambda > $homeroad/${rawdata1}_val_1_bismark_bt2_pe_chrLambda.sam \
&& perl /media/data2/chenwei_RRBS_project/MethylExtract_1.8.1/MethylExtractBSCR.pl seqFile=$chrLambda  inFile=$homeroad/${rawdata1}_val_1_bismark_bt2_pe_chrLambda.sam  flagW=99,147  flagC=83,163  >  $homeroad/${rawdata1}_val_1_bismark_bt2_pe_Bisulfy_Conversion_rate_report 
#&& $qualimap bamqc -bam $homeroad/${rawdata1}_val_1.fq.gz_bismark_bt2_pe.sort.bam  -outfile  $homeroad/${rawdata1}_1.mapreport.pdf 

#
block
rm $homeroad/*_context_*_1_val_1_bismark_bt2_pe.txt    $homeroad/*_1_val_1_bismark_bt2_pe.bedGraph.gz $homeroad/*_1_val_1_bismark_bt2_pe.bismark.cov.gz   $homeroad/*_1_val_1_bismark_bt2_pe_splitting_report.txt $homeroad/*_1_val_1_bismark_bt2_pe.cytosine_context_summary.txt # $homeroad/*_1_val_1.fq.gz $homeroad/*_2_val_2.fq.gz
mkdir $output_road/${sample_name}
mv   $homeroad/*bismark*  $homeroad/*.txt $output_road/${sample_name}
