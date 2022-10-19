#!/bin/bash
homeroad=$1
mnb=`ls $homeroad`
echo $mnb
for i in $mnb
do
 cd $homeroad/$i
 md5sum -c  md5.txt
done

