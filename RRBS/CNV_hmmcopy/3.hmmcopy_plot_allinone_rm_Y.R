#function reference：http://developer.51cto.com/art/201305/393439.htm
rm(list = ls())
setwd("/mnt/data/chenwei/huahua/3.mapped_data/huahua_HMM_copy_result")
##ls()返回global environment 里面的所有object的名字。 是一个character vector。
#rm()清除全部对象
options(bitmapType="cairo")
#环境设置函数为options(),用options()命令可以设置一些环境变量
#在ui.R 或者server.R中添加 对bitmapType的设置
#library(shiny)
#options(bitmapType=”cairo”)
chrlen <- read.table("/mnt/data/chenwei/qinmen_BR/0.meng_script/1.meth_script/CNV_1Mb_script/hg38.new_chrlength2.txt",header=T,row.names=1)
head(chrlen)

tlength <- sum(chrlen[1,])
#sum(x):x中各元素的加和

BIN="1Mb"
LBIN=1000000

sbin=50000000 # number refer to seperating chromosomes in figure
bin=BIN
lbin=LBIN

#filenum=73
#pdf(file="/public2/chenwei/Chip-wei/wholeblood/Chip-seq1/CNV_practice/hmmcopy/cw.hmmcopy.cnv.1Mb.pdf",width=30,height=filenum*2)
#filenum=20
filenum=60
pdf(file="huahua_60_sample_cnv.hmmcopy.1M.all_in_one_remove_Y.pdf",width=30,height=filenum*0.5)
#pdf(file="Mengyang_20_sample_cnv.hmmcopy.1M.all_in_one.pdf",width=30,height=filenum*2)

#opar <- par(no.readonly = TRUE)
#par(opar)
#mai 以数值向量表示的边界大小，顺序为“下、左、上、右”，单位为英寸
#mar 以数值向量表示的边界大小，顺序为“下、左、上、右”，单位为英分*。
#默认值为c(5, 4, 4, 2) + 0.1 #*一英分等于十二分之一英寸。——译者注
par(mfrow=c(1,1),mar=c(3,1,1,1))
xmax = sum(chrlen[rownames(chrlen)==bin,]) + 22*sbin
ymax = 4
par(mar=c(1,1,1,1))
plot(0,0,ylim=c(-4,(ymax + 1) * filenum),pch=21,cex=0.00001, xlim = c(-xmax/32,xmax),ann=FALSE,axes=FALSE)
start= 0; col=1; color = c("midnightblue","indianred")
for(jj in c(1:22,"X")){
  j=paste("chr",jj,sep="")
  segments(start,-2,start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j],-2,col=color[col], lwd = 4)
  #segments(x0, y0,x1, y1)从(x0,y0)各点到(x1,y1)各点画线段
  col = 3-col
  text(start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j]/2, -3, jj , adj = 0.5,cex=0.6)
  #text函数用来在一张图表上添加文字，只需要指定对应的x和y坐标，以及需要添加的文字内容就可以了
  #text(x = 3, y = 3, labels = "text")
  #adj : 调整文字的位置，一个值时调整的是x轴的位置，如果为两个值时，第一个调整的是x轴的位置，第二个调整的是y轴的位置，可选范围为[0, 1]
  start = start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j] + sbin
}

row = 0
#kids_sample_list<-c("10-3-1E","16-3-4D","17-3-2F","18-3","19-3-4E","20Q-4H","6-3","7-3","8-3","E16C",
#                    "11-3-3F","12-3-4A","13-3-1G","14-3","15-3-2C","1Q","2-3-3C","3-3","4-3","5-3")
#father_sample_list<-c("10-2-1D","16-2-4C","17-2-2E","18-2","19-2-2H","20F-4G","6-2","7F","8-2","E16F",
#                      "11-2","12-2-3H","13-2-1H","14-2","15-2","1F","2-2","3F","4-2","5-2")
#mother_sample_list<-c("10-1-in-9","16-1","17-1-2D","18-1","19-1-2G","20M-4F","6-1","7-1","8-1","E16M",
#                      "11-1","12-1-3G","13-1-1F","14-1","15-1","1M","2-1","3M","4-1","5-1")
huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
huahua_meta$analysis_name<-factor(huahua_meta$analysis_name,levels=c(paste0("AMA_K_",1:10),paste0("YOUNG_K_",1:10),
                                                                     paste0("AMA_M_",1:10),paste0("YOUNG_M_",1:10),
                                                                     paste0("AMA_F_",1:10),paste0("YOUNG_F_",1:10)))
huahua_meta<-huahua_meta[order(huahua_meta$analysis_name),]
##samples<-rev(c(kids_sample_list,mother_sample_list,father_sample_list))
samples<-rev(huahua_meta$library_code)
length(samples)

for (CWname in samples){
  
  s = paste0(CWname, ".hmmcopy.win1000000.corrected.csv")
  #s="hmmcopy.win1000000.corrected.csv"
  #s="/public2/chenwei/Chip-wei/wholeblood/Chip-seq1/CNV_practice/hmmcopy/D.hmmcopy.win10M.corrected.csv"
  cnv <- read.csv(s, header = T)
  cnv$color <- cnv$copy
  color = c("forestgreen",'dodgerblue4',"brown2","white")
  col = 1; start = 0; column = 12; zoom = 0.5; ymax = 4
  for(j in c(paste("chr",c(1:22,"X"), sep=""))){
    
    i <- which(cnv$valid & cnv$ideal & (cnv$space ==j) & !is.na(cnv$copy))
    #valid Bins with valid GC and average mappability and non-zero read
    #ideal Valid bins of high mappability and reads that are not outliers
    d <- cnv[i,]
    
    d[d$copy > 0.44, 13] <- 3
    d[d$copy < -0.57, 13] <- 1
    d[d$copy >= -0.57 & d$copy <= 0.44, 13] <- 2
    
    points(d[,2] + start, d[,column] + 1 + (ymax + 1)* row,pch=16,cex=zoom,col = color[d[,13]])
    
    #lines(cnv[i,3] + start, cnv[i,5] + 1 + (ymax + 1)* row, lwd = 1.5, col = "black")
    
    start = start + chrlen[rownames(chrlen)==bin,colnames(chrlen)==j] + sbin
  } 
  axis(side=2,at=c((ymax + 1)* row, 1 + (ymax + 1)* row, 2+ (ymax + 1)* row), tcl=-0.3,labels=c( "-1","0","1"),pos=0, cex.axis=0.5,mgp=c(0.6,0.2,0))
  #axis(side, vect)画坐标轴,side=1时画在下边,side=2时画在左边,side=3时画在上边,side=4时画在右边。
  #可选参数at指定画刻度线的位置坐标
  
  #mgp设置坐标轴标题，坐标标记和坐标轴边界宽度。默认值为c(3,1,0)，mgp[1]影响坐标轴标题，mgp[2,3]影响坐标标记和坐标轴。
  segments(0, -1 + (ymax + 1)* row,xmax,-1 + (ymax + 1)* row,lty="dotted")
  segments(0,  (ymax + 1)* row,xmax,(ymax + 1)* row,lty="dotted")
  segments(0, 0.5 + (ymax + 1)* row,xmax,0.5+(ymax + 1)* row,lty="dotted")
  segments(0,  1 + (ymax + 1)* row,xmax,1 + (ymax + 1)* row,lty="dotted")
  segments(0, 1.5 + (ymax + 1)* row,xmax,1.5 + (ymax + 1)* row,lty="dotted")
  segments(0, 2 + (ymax + 1)* row,xmax,2 + (ymax + 1)* row,lty="dotted")
  segments(0, 2.5 + (ymax + 1)* row,xmax,2.5 + (ymax + 1)* row,lty="dotted")
  
  rect(0, -1 - 0.1 + (ymax + 1)* row,xmax,2.5 + 0.1 + (ymax + 1)* row)
  # rect(x1, y1, x2, y2)绘制长方形,(x1, y1)为左下角,(x2,y2)为右上角
  text(-xmax/24-1,1 + (ymax + 1)* row, CWname, cex=1)
  row = row + 1
}
#par(opar)
dev.off()