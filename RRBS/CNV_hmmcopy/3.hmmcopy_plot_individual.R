rm(list = ls())
options(bitmapType="cairo")
setwd("/mnt/data/chenwei/qinmen_BR/3.methy_result/AMA_project/CNV_result/HMM_copy_result")
chrlen <- read.table("/mnt/data/chenwei/qinmen_BR/0.meng_script/1.meth_script/CNV_1Mb_script/hg38.new_chrlength2.txt",header=T,row.names=1)
head(chrlen)

tlength <- sum(chrlen[1,])

BIN="1Mb"
LBIN=1000000

sbin=0 # sep the chromosome in fig
bin=BIN
lbin=LBIN


pdf(file="A_P_cnv.hmmcopy.1M.individual.pdf",width=10,height=2)

#samples = read.table("../samples.txt",header=F)
#samples = samples$V1
#samples = rev(c("A1","A10","A11","A21","A22","A23","A24","A25","A3","A31","A32","A33","A34","A35","A36","A37","A4","A6","A8","A9","P11","P12","P13","P14","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24","P25","P31","P32","P33","P34","P35","P36","P37","P38"))
samples<-c("A1")
for (embryo in samples){
  embryo<-"A1"
  par(mfrow=c(1,1),mar=c(3,1,1,1))
  xmax = sum(chrlen[rownames(chrlen)==bin,]) + 23*sbin
  ymax = 2
  par(mar=c(1,1,1,1))
  plot(0,0,ylim=c(-4,(ymax + 1) * 2),pch=21,cex=0.00001, xlim = c(-xmax/32,xmax),ann=FALSE,axes=FALSE)
  start= 0; col=1; color = c("midnightblue","indianred")
  for(jj in c(1:22,"X","Y")){
    j=paste("chr",jj,sep="")
    segments(start,-2,start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j],-2,col=color[col], lwd = 4)
    col = 3-col
    text(start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j]/2, -3, jj , adj = 0.5,cex=0.8)
    start = start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j] + sbin
  }
  
  s = paste0(embryo,".hmmcopy.win1000000.corrected.csv")
  
  cnv <- read.csv(s, header = T)
  cnv$color <- cnv$copy
  color = c("forestgreen",'dodgerblue4',"brown2","white")
  col = 1; start = 0; column = 12; zoom = 0.3; ymax = 4
  for(j in c(paste("chr",c(1:22,"X","Y"), sep=""))){
    i <- which(cnv$valid & cnv$ideal & (cnv$space ==j) & !is.na(cnv$copy))
    d <- cnv[i,]
    
    d[d$copy > 0.44, 13] <- 3
    d[d$copy < -0.57, 13] <- 1
    d[d$copy >= -0.57 & d$copy <= 0.44, 13] <- 2
    
    points(d[,2] + start, d[,column] + 1,pch=16,cex=zoom,col = color[d[,13]])
    
    #lines(cnv[i,3] + start, cnv[i,5] + 1 + (ymax + 1)* row, lwd = 1.5, col = "black")
    
    start = start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j] + sbin
  } 
  axis(side=2,at=c(0, 1 , 1.5 ), tcl=-0.3,labels=c( "-1","0","0.5"),pos=0, cex.axis=0.5,mgp=c(0.6,0.2,0))
  segments(0, -1, xmax,-1, lty="dotted")
  segments(0,  0, xmax,0, lty="dotted")
  segments(0,  1, xmax, 1,lty="dotted")
  segments(0, 1.5, xmax, 1.5, lty="dotted")
  segments(0, 2, xmax, 2, lty="dotted")
  
  rect(0, -1 - 0.1, xmax,2 + 0.1)
  text(xmax/2, 4 , embryo, cex=1.2)
  text(-xmax/40, 0.5 , "Log2(copy ratio)", cex=0.8,xpd=TRUE, srt=90)
  
}

dev.off()

