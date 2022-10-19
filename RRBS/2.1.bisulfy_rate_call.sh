#!/bin/bash
homeroad=$1


#Version=2
for id in    $homeroad/*
do
 echo $id
# cat ${id}|grep "Bisulfite conversion rate" |awk '{if($4 < 0.9885){print $0} }'
 cat ${id}|grep "Bisulfite conversion rate" |awk '{print $0}'

done

