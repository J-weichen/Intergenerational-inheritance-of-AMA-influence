#!/bin/bash
homeroad=$1
name1=$2
name2=$3
cd $homeroad
mv ${name1}_1_CpG.bismark.cov.gz  ${name2}_CpG.bismark.cov.gz

