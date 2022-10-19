sampledir=$1
window=$2

Counter=/mnt/data/chenwei/software/HMMcopy/readCounter
#fot first
#outdir=/mnt/data/chenwei/qinmen_BR/3.methy_result/AMA_project/CNV_result/HMM_copy_result
#fir varification
outdir=/mnt/data/chenwei/huahua/4.methy_result/4.CNV
for bam in $sampledir/*/*_1_val_1_bismark_bt2_pe_sort.bam
do
	temp=`dirname $bam`
        sample=`basename $temp`
        mkdir -p ${outdir}/${sample}
        echo $sample
        $Counter -w $window -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $bam > ${outdir}/${sample}/${sample}.window${window}.readcounts.wig
done
