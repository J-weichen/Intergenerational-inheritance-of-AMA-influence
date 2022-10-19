#!/bin/bash
#A1  A10  A11  A21  A22  A23  A24  A25  A3  A31  A32  A33  A34  A35  A36  A37  A4  A6  A8  A9  P11  P12  P13  P14  P15  P16  P17  P18  P19  P20  P21  P22  P23  P24  P25  P31  P32  P33  P34  P35  P36  P37  P38

#list="14-1 15-2-2B  16-3-4D  19-1-2G  1M      20Q-4H  3F     4M-5A  5-3    7F     8-2     9-3-1C 10-2-1D    11-3-3F  13-1-1F  14-2     15-3-2C  17-1-2D  19-2-2H  1Q      2-1-3A  3M     4Q-5C  6F     7M 8-3 10-3-1E    12-1-3G  13-2-1H  14-3     16-1-4B  17-2-2E  19-3-4E  20F-4G  2-2-3B  3Q     5-1    6M     7Q-5E  9-1-1A 11-1-3D    12-2-3H  13-3-1G  15-1-2A  16-2-4C  17-3-2F  1F       20M-4F  2-3-3C  4F-5B  5-2    6Q-5D  8-1    9-2-1B"

list="14-1  14-2  14-3  18-1  18-2  18-3"
road=/mnt/data/chenwei/huahua/1.rawdata/merge_fastq/huhua_anno_0413

mv ${road}/*/*_Conversion_rate_report ${road}/*/*_bismark_bt2_PE_report.txt  /mnt/data/chenwei/huahua/1.rawdata/0.QC_file/Bisulfy_Conversion_file
mv ${road}/*/*_fastqc.html ${road}/*/*_fastqc.zip  /mnt/data/chenwei/huahua/1.rawdata/0.QC_file/fastqc_report_file
mv ${road}/*/*_trimming_report.txt  /mnt/data/chenwei/huahua/1.rawdata/0.QC_file/trimming_report_file

for sample in $list;
do
   mkdir /mnt/data/chenwei/huahua/2.cleandata/${sample}
   mv ${road}/${sample}/*_val_*.fq.gz /mnt/data/chenwei/huahua/2.cleandata/${sample}
done
