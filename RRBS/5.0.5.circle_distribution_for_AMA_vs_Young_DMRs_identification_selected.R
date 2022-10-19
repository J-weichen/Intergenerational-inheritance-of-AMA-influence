#Circular Visualization in R
#https://www.jianshu.com/p/deb59290ad96?utm_campaign=maleskine&utm_content=note&utm_medium=seo_notes&utm_source=recommendation
#https://jokergoo.github.io/circlize_book/book/index.html
#美化：https://www.jianshu.com/p/a87bcc1cb67b
#其他参考：https://zhuanlan.zhihu.com/p/113758594
#https://zhuanlan.zhihu.com/p/103132265
rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(grid)
library(gridExtra)
library(circlize)
library(ggsci)
library(scales)
library(data.table)
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)
#FOR kids
#read meth data 
compare_name<-"kids_AMA_vs_Young";split_region<- "200bp";depth<-"6X"
#data preparation
kids_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
head(kids_DMR_all_big)
kids_DMR_ALL_bed_raw<-kids_DMR_all_big[,c("seqnames","start","end","meth.diff")]
kids_DMR_ALL_bed_raw$meth.diff<-as.numeric(kids_DMR_ALL_bed_raw$meth.diff)
kids_DMR_ALL_bed_raw$start<-as.numeric(kids_DMR_ALL_bed_raw$start)
kids_DMR_ALL_bed_raw$end<-as.numeric(kids_DMR_ALL_bed_raw$end)
colnames(kids_DMR_ALL_bed_raw)<-c("chr","start","end","value")
str(kids_DMR_ALL_bed_raw)
kids_DMR_ALL_bed <- data.frame(as.data.table(kids_DMR_ALL_bed_raw))
kids_DMR_ALL_bed$col<-kids_DMR_ALL_bed$chr
kids_DMR_ALL_bed$col<-ifelse((kids_DMR_ALL_bed$value>0), ppCor[[3]], ppCor[[4]])
head(kids_DMR_ALL_bed)

#for mother


#甲基化对应基因列表的读取
#FOR mom
#read meth data 
compare_name<-"mother_AMA_vs_Young";split_region<- "200bp";depth<-"6X"
#data preparation
mother_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
head(mother_DMR_all_big)
mother_DMR_ALL_bed_raw<-mother_DMR_all_big[,c("seqnames","start","end","meth.diff")]
mother_DMR_ALL_bed_raw$meth.diff<-as.numeric(mother_DMR_ALL_bed_raw$meth.diff)
mother_DMR_ALL_bed_raw$start<-as.numeric(mother_DMR_ALL_bed_raw$start)
mother_DMR_ALL_bed_raw$end<-as.numeric(mother_DMR_ALL_bed_raw$end)
colnames(mother_DMR_ALL_bed_raw)<-c("chr","start","end","value")
str(mother_DMR_ALL_bed_raw)
mother_DMR_ALL_bed <- data.frame(as.data.table(mother_DMR_ALL_bed_raw))
mother_DMR_ALL_bed$col<-mother_DMR_ALL_bed$chr
mother_DMR_ALL_bed$col<-ifelse((mother_DMR_ALL_bed$value>0), ppCor[[3]], ppCor[[4]])
head(mother_DMR_ALL_bed)

#for father
#read meth data 
compare_name<-"father_AMA_vs_Young";split_region<- "200bp";depth<-"6X"
#data preparation
father_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
head(father_DMR_all_big)
father_DMR_ALL_bed_raw<-father_DMR_all_big[,c("seqnames","start","end","meth.diff")]
father_DMR_ALL_bed_raw$meth.diff<-as.numeric(father_DMR_ALL_bed_raw$meth.diff)
father_DMR_ALL_bed_raw$start<-as.numeric(father_DMR_ALL_bed_raw$start)
father_DMR_ALL_bed_raw$end<-as.numeric(father_DMR_ALL_bed_raw$end)
colnames(father_DMR_ALL_bed_raw)<-c("chr","start","end","value")
str(father_DMR_ALL_bed_raw)
father_DMR_ALL_bed <- data.frame(as.data.table(father_DMR_ALL_bed_raw))
father_DMR_ALL_bed$col<-father_DMR_ALL_bed$chr
father_DMR_ALL_bed$col<-ifelse((father_DMR_ALL_bed$value>0), ppCor[[3]], ppCor[[4]])
head(father_DMR_ALL_bed)


pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/circle_distribution_for_kids_mom_dad_DMRs.pdf",width =5,height =5)
grid.newpage(); 
circos.clear()# 结束上次绘图绘图，否则会提示警告信息
circos.par("start.degree" = 80)#设置1号染色体起始角度，然后顺时针
circos.par("gap.degree" = c(rep(c(2, 4), 11),2)) 
circos.par("gap.degree" = c(rep(2,21),20)) 
#circos.par(gap.after=c(rep(c(2, 4), 11),2),"start.degree" = 80)

circos.initializeWithIdeogram(species = "hg38", chromosome.index = paste0("chr", c(1:22)))#除去contig
#circos.genomicDensity(kids_DMR_ALL_bed, col = c(ppCor[5]), track.height = 0.05)#密度图
circos.genomicTrack(kids_DMR_ALL_bed, panel.fun = function(region,value, ...) {circos.genomicPoints(region, value,pch = 16, cex = abs(value$value/100),col = value$col, ...)})#col = 'red',在这个轨道里画散点，散点图只需要一列value
circos.genomicTrack(mother_DMR_ALL_bed, panel.fun = function(region,value, ...) {circos.genomicPoints(region, value,pch = 16, cex = abs(value$value/100),col = value$col, ...)})
circos.genomicTrack(father_DMR_ALL_bed, panel.fun = function(region,value, ...) {circos.genomicPoints(region, value,pch = 16, cex = abs(value$value/100),col = value$col, ...)})
dev.off() 
