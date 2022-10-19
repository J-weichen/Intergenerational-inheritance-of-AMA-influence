#!/bin/bash
homeroad=$1
mnb=`ls $homeroad`
echo $mnb
for id in $mnb
do
 echo $id
# cat ${id}|grep "Bisulfite conversion rate" |awk '{if($4 < 0.9885){print $0} }'
 zcat $homeroad/${id}/${id}_1_val_1.fq.gz|head -n1 |awk '{print $0}'
done

