#! /media/data4/chenwei/script_10X/software_10X/R-3.5.2/bin/Rscript --vanilla
rm(list = ls())
options(bitmapType="cairo")

setwd("/mnt/data/chenwei/qinmen_BR/3.methy_result/AMA_project/CNV_result")
chrlen <- read.table("/mnt/data/chenwei/qinmen_BR/0.meng_script/1.meth_script/CNV_1Mb_script/hg38.new_chrlength2.txt",header=T,row.names=1)
head(chrlen)
tlength <- sum(chrlen[1,])

BIN="1Mb"
LBIN=1000000

sbin=0 # sep the chromosome in fig
bin=BIN
lbin=LBIN

filenum=43
sample_road="/mnt/data/chenwei/qinmen_BR/2.pre_data/2.pre_RRBS_data/3.mapped_data/bam_file"
pdf(file="/mnt/data/chenwei/qinmen_BR/3.methy_result/AMA_project/CNV_result/AMA_and_Young.cnv.1M.pdf",width=20,height=filenum)

par(mfrow=c(1,1),mar=c(3,1,1,1))
xmax = sum(chrlen[rownames(chrlen)==bin,]) + 23*sbin
ymax = 4
par(mar=c(1,1,1,1))
plot(0,0,ylim=c(-2,(ymax + 1) * filenum),pch=21,cex=0.00001, xlim = c(-xmax/32,xmax),ann=FALSE,axes=FALSE)
start= 0; col=1; color = c("midnightblue","indianred")
for(jj in c(1:22,"X","Y")){
  j=paste("chr",jj,sep="")
  segments(start,-1,start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j],-1,col=color[col], lwd = 4)
  col = 3-col
  text(start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j]/2, -2, jj , adj = 0.5,cex=0.6)
  start = start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j] + sbin
}

row = 0

samples = c("A1","A10","A11","A21","A22","A23","A24","A25","A3","A31","A32","A33","A34","A35","A36","A37","A4","A6","A8","A9","P11","P12","P13","P14","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24","P25","P31","P32","P33","P34","P35","P36","P37","P38")
for (embryo in samples){
    
  s = paste(sample_road,embryo, "chrs_1M/all.chrs.normalized.win1M.hiddenseq", sep = "/")
  #s = paste(embryo, "chrs/all.chrs.no.normalised.depth", sep = "/")
  
  cnv <- read.table(s, header = F)

  col = 1; start = 0; column = 4; zoom = 0.3; ymax = 4
  for(j in c(paste("chr",c(1:22,"X","Y"), sep=""))){
    i = which(((cnv[,1]==j) * (cnv[,column] != 0)) == 1)
    points(cnv[i,3] + start, cnv[i,column] + 1 + (ymax + 1)* row,pch=16,cex=zoom,col = color[col])
    
    #points(cnv[i,5] + start, cnv[i,column] + (ymax + 1)* row,pch=16,cex=zoom,col = color[col])
    #lines(cnv[i,3] + start, cnv[i,5] + 1 + (ymax + 1)* row, lwd = 1.5, col = "black")
    col = 3 - col
    start = start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j] + sbin
  } 
  axis(side=2,at=c(1 + (ymax + 1)* row, 3 + (ymax + 1)* row, 5 + (ymax + 1)* row), tcl=-0.3,labels=c("0","2","4"),pos=0, cex.axis=0.6,mgp=c(0.6,0.2,0))
  segments(0,2 + (ymax + 1)* row,xmax,2 + (ymax + 1)* row,lty="dotted")
  segments(0,3 + (ymax + 1)* row,xmax,3 + (ymax + 1)* row,lty="dotted")
  segments(0,4 + (ymax + 1)* row,xmax,4 + (ymax + 1)* row,lty="dotted")
  rect(0,1 - 0.1 + (ymax + 1)* row,xmax,5 + 0.1 + (ymax + 1)* row)
  
  text(-xmax/24,3 + (ymax + 1)* row, embryo, cex=0.6)
  row = row + 1
}

dev.off()


