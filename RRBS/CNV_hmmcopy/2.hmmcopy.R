##注意：不同批次的需要单一批次单独进行矫正
rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

setwd("/mnt/data/chenwei/huahua/3.mapped_data/huahua_HMM_copy_result/")
library(HMMcopy)
samples<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/huahua_HMM_copy_result")
#samples =c("E16C","E16F","E16M")
gfile <- "/mnt/data/chenwei/software/HMMcopy/hg38/hg38.gc.1000000.wig"
mfile <- "/mnt/data/chenwei/software/HMMcopy/hg38/hg38.map.1000000.bw"

for (name in samples){
  sample.name<-name
  #sample.name<-"10-1-in-9"#test line
  print(sample.name)
  rfile <-paste("/mnt/data/chenwei/huahua/3.mapped_data/huahua_HMM_copy_result/",sample.name,"/",sample.name,".window1000000.readcounts.wig",sep="")

  normal_reads <- wigsToRangedData(rfile, gfile, mfile)
  normal_copy <- correctReadcount(normal_reads)
  normal_segment <- HMMsegment(normal_copy)
  
  ############################################################################
  pdf(file = paste0("/mnt/data/chenwei/huahua/3.mapped_data/huahua_HMM_copy_result/",sample.name,"_chrs.pdf"))
  
  for (i in c(1:22,"X")){
    par(mar = c(4, 4, 2, 0))
    plotCorrection(normal_copy, chr = paste0("chr",i), pch = ".")
  }
  
  dev.off()
  
  write.csv(normal_copy,file = paste0("/mnt/data/chenwei/huahua/3.mapped_data/huahua_HMM_copy_result/",sample.name,".hmmcopy.win1000000.corrected.csv"),row.names = F,quote = F)
  write.csv(normal_segment$segs,file = paste0("/mnt/data/chenwei/huahua/3.mapped_data/huahua_HMM_copy_result/",sample.name,".hmmcopy.win1000000.segment.csv"),row.names = F,quote = F)
} 
