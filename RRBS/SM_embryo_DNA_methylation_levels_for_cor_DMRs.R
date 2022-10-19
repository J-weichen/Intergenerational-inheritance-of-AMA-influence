rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(Seurat)
library(ggplot2)
library(ggpubr) 
library(methylKit)
library(data.table)
library(scales)
library(ggsci)
library(stringr)
library(pheatmap)
##提供自定义函数

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # 计算长度
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # 以 groupvars 为组,计算每组的长度,均值,以及标准差
  # ddply 就是 dplyr 中的 group_by + summarise
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # 重命名  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  # 计算标准偏差
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  # 计算置信区间
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#step2: read parental kids correlated DMRs
#for father 
data_analysis <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_father.txt",header=T,sep="\t")
data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region"),]
plot_cordata<-data_analysis_sig[,c("DMR_region","kids","parent","Family")]
father_cor_DMRs<-unique(as.character(plot_cordata$DMR_region))

#for mother 
data_analysis <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_mother.txt",header=T,sep="\t")
data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region"),]
plot_cordata<-data_analysis_sig[,c("DMR_region","kids","parent","Family")]
mother_cor_DMRs<-unique(as.character(plot_cordata$DMR_region))
length(father_cor_DMRs);length(mother_cor_DMRs)#52 54

#all_cor_DMRs<-unique(mother_cor_DMRs)
all_cor_DMRs<-unique(c(father_cor_DMRs,mother_cor_DMRs))
length(all_cor_DMRs)#102
all_cor_DMRs_bed<-gsub("\t","_",all_cor_DMRs)
FF_bed<-as.data.frame(str_split_fixed(all_cor_DMRs, "_", n = 3))

colnames(FF_bed)<-c("chr","start","end")
FF_bed<-FF_bed[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",FF_bed$chr,invert=TRUE),]
FF_bed <- data.frame(FF_bed,V4=paste(FF_bed$chr,FF_bed$start,FF_bed$end,sep="_"))
FF_bed$chr<-factor(FF_bed$chr,levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
FF_bed<-FF_bed[order(FF_bed$chr,decreasing = F),]
FF_bed$chr<-as.character(FF_bed$chr)
FF_bed$start<-as.numeric(as.character(FF_bed$start))
FF_bed$end <-as.numeric(as.character(FF_bed$end ))
#FF_bed$V4<-as.character(FF_bed$V4)
head(FF_bed);dim(FF_bed)

#extend 200 bins to 500 bins
FF_bed$start<-FF_bed$start-150;FF_bed$end<- FF_bed$end + 150
FF_bed$V4 <- as.character(paste(FF_bed$chr,FF_bed$start,FF_bed$end,sep="_"))

head(FF_bed);dim(FF_bed)
str(FF_bed)

#step3 reading single cell methylation files
SM_files <-as.data.frame(read.csv(file="/mnt/data/chenwei/huahua/6.SM_embryo_meth/SM_sample_list.csv",header=T,sep=","))
SM_files2 <-SM_files$data_name
##step4: calculate DNA methylation level for selected regions 
sample_names<-c();interSE<-list()
SM_files2[which(SM_files$sample_name == "blastocyst-ZHWX-E1-043")]

for (file in SM_files2) {
  # file<-"/mnt/data/kongsiming/zhaifan/covFiles/blastocyst/WYQL-E2/WYQL-E2-101.CG.txt.gz"                                                         
  #print(file)
  SM_names<-SM_files[which(SM_files$data_name == file),]$sample_name
  print(as.character(SM_names))
  
  #data_meth1<-read.table(gzfile(file),header = F,fill=TRUE, na.strings = "",stringsAsFactors=F)
  data_meth1 = fread(file,header = F,fill=TRUE, na.strings = "",stringsAsFactors=F)
  head(data_meth1);dim(data_meth1)
  data_meth1$meth<-data_meth1$V4/(data_meth1$V4+data_meth1$V5)
  data_meth1<-data_meth1[,c("V1","V2","meth")]
  colnames(data_meth1)<-c("chr","start","meth")
  data_meth1$end<-data_meth1$start
  head(data_meth1)
  data_meth2<-data_meth1[,.(chr,start,end,meth)]
  #data_meth2<-as.data.table(data_meth2)
  data_meth2<-na.omit(data_meth2)
  str(data_meth2)
  data_meth2$start<-as.numeric(data_meth2$start)
  data_meth2$end<-as.numeric(data_meth2$end)
  data_meth2$meth<-as.numeric(data_meth2$meth)
  head(FF_bed)
  FF2<-as.data.table(FF_bed[,c(1:3)])
  colnames(FF2)<-c("chr","start","end")
  dim(data_meth2);dim(FF2)
  # head(data_meth2)
  setkey(FF2, chr, start, end)
  data_meth3<-as.data.frame(foverlaps(data_meth2,FF2,nomatch=0))
  dim(data_meth3);head(data_meth3)
  data_meth3<-data_meth3[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",data_meth3$chr,invert=TRUE),]
  #dim(data_meth3);head(data_meth3)
  data_meth3$ID <- paste(data_meth3$chr,data_meth3$start,data_meth3$end,sep = "_")
  ##select  region with at least 3 CpGs
  CpG_num<-table(data_meth3$ID)
  data_meth4<-data_meth3[which(data_meth3$ID %in% names(CpG_num[CpG_num>=3])),]
  head(data_meth4)
  if(nrow(data_meth4)==0){next;}
  data_meth5 <-aggregate(data_meth4[,6],list(data_meth4[,7]),FUN=mean, na.rm=TRUE, na.action=NULL)
  colnames(data_meth5)<-c("region","meth")
  interSE<-c(interSE,list(as.data.frame(data_meth5)))
  sample_names<-c(sample_names,SM_names)
  data_meth1<-data_meth2<-data_meth3<-data_meth4<-data_meth5<-0
} 

length(interSE);length(sample_names)
names(interSE)<- sample_names
saveRDS(interSE,file = "/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/SM_single_Cell_meth_data_in_two_corDMRs_extend_500bp.rds")

mymergeddata<-Reduce(function(x,y) {merge(x,y,all=T,by="region")}, interSE)
colnames(mymergeddata)<-c("region",sample_names)
head(mymergeddata);dim(mymergeddata)
saveRDS(mymergeddata,file = "/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/data_merge_SM_single_Cell_meth_data_in_two_corDMRs_extend_500bp.rds")
#rm(list = ls())
mymergeddata<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/data_merge_SM_single_Cell_meth_data_in_two_corDMRs_extend_500bp.rds")
#which(colnames(mymergeddata)=="Reference_CN")
colnames(mymergeddata)
rownames(mymergeddata)<-mymergeddata$region
mymergeddata<-mymergeddata[,-1]

rownames(mymergeddata) ## "chr22_20988051_20988550"

#show methlation level in heatmap in global
SM_meta <-SM_files
head(SM_meta);dim(SM_meta)
rownames(SM_meta)<-SM_meta$sample

#colData构建
colData <-SM_meta[colnames(mymergeddata),]
dim(colData);length(colnames(mymergeddata))
#colData<-colData[-c(which(rownames(colData)=="Reference_CN")),]
head(colData)
table(colData$stage)
colData$stage<-factor(colData$stage,levels = c("2cell","4cell","8cell","morula","blastocyst"))
colData<-colData[order(colData$stage,decreasing = T),]
mymergeddata<-mymergeddata[,rownames(colData)]

head(mymergeddata)
range(mymergeddata)

mymergeddata<-mymergeddata*100
dim(mymergeddata);dim(colData)
pheatmap(mymergeddata, cluster_rows=F, cluster_cols =F,scale="row")
pheatmap(mymergeddata, cluster_rows=F, cluster_cols =F,scale="none")

plot_sample<-data.frame(colMeans(mymergeddata,na.rm = T))
plot_sample$sample<-rownames(plot_sample);colnames(plot_sample)<-c("meth","sample")
plot_sample2<-merge(plot_sample,colData,by=0)
head(plot_sample2)
plot_sample2<-plot_sample2[order(plot_sample2$stage,decreasing = F),]
table(plot_sample2$stage)
plot_sample2$sample<-as.character(plot_sample2$sample)
plot_sample2$sample<-factor(plot_sample2$sample,levels = as.character(plot_sample2$sample))

ggplot(data = plot_sample2, mapping = aes(x = sample, y = meth, fill = stage))+ 
  geom_bar(stat ='identity', position = 'dodge')+scale_color_manual(values=c(ppCor))+NoLegend()+
  theme_classic()+labs(x="stage",y="DNA meth level",title="Early embryo DNA methlation level in cor_DMRs")+
  theme(axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1,angle = 90))


colData2<-colData[,c("stage","stage")];colnames(colData2)<-c("stage1","stage2")
pheatmap(mymergeddata, cluster_rows=F, cluster_cols =F, annotation_col=colData2)

Stages<-c("2cell","4cell","8cell","morula","blastocyst")
colData2$stage1<-factor(colData2$stage1,levels=Stages)
colData2<-colData2[order(colData2$stage1,decreasing = F),]
mymergeddata<-mymergeddata[,rownames(colData2)]

#设定 annotations
# 生成行 列的注释
#annotation_col<-data.frame(Treatment = factor(c(rep("Young", 12),rep("AMA", 10))))
#rownames(annotation_col) = colnames(heat_plot)
#annotation_row = data.frame(Class = factor(rep(c("UP_DMRs","Down_DMRs"), c(length(which(myDiff.all$meth.diff>0)), length(which(myDiff.all$meth.diff<0))))))
#rownames(annotation_row) = rownames(myDiff.all)

## 自定义分组的颜色
#ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(UP_DMRs=ppCor[3],Down_DMRs=ppCor[4]))

p1<-pheatmap(mymergeddata, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = colData2, #annotation_row=annotation_row,
             #  annotation_colors = ann_colors, 
             #  gaps_row = length(which(data_region_merge_DMR_all$meth.diff>0)),
             #   gaps_col =(length(colnames(heat_plot))+2)/2,
             # cutree_col = 2,treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","orange","firebrick3"))(50),
             main ="DNA methylation level of DMRs in early embryo development(original value)",angle_col ="90")
#scale by row
p2<-pheatmap(mymergeddata, cluster_row =T,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = colData2,#annotation_row=annotation_row,
             # annotation_colors = ann_colors, 
             # gaps_row = length(which(data_region_merge_DMR_all$meth.diff>0)),gaps_col =(length(colnames(heat_plot))+2)/2,cutree_col = 2,
             # treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             scale ="row", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(50),
             main ="DNA methylation level of DMRs in early embryo development:scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/SM_DNA_methylation_level_of_DMRs_in_early_embryo_development_heatmap_500bp_extend_1.pdf",width=16,height=8)
print(p1)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/SM_DNA_methylation_level_of_DMRs_in_early_embryo_development_heatmap_500bp_extend_2.pdf",width=16,height=8)
print(p2)
dev.off()
dev.new()

#step four : re-orgenization meth data 
my_data<-as.data.frame(t(mymergeddata),stringsAsFactors = F)
#my_data<-my_data[-c(which(rownames(my_data)=="Reference_CN")),]
my_data$sample<-rownames(my_data)
my_data2 <- melt(my_data)
head(my_data2)
colData2$sample<-rownames(colData2)
my_data3<-merge(my_data2,colData2,by="sample")
head(my_data3)
head(my_data3,n=7)
dim(my_data3)#85578     5
write.csv(my_data3,"/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/data_merge_SM_single_Cell_meth_data_in_two_corDMRs_add_meta_500bp_extend.csv")

##step five:  plot target regions DNA methylation
##for all region
plot_data <-read.csv("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/data_merge_SM_single_Cell_meth_data_in_two_corDMRs_add_meta_500bp_extend.csv",row.names = 1)
head(plot_data,n=7)
names(table(plot_data$stage2))

Stage<-c("2cell","4cell","8cell","morula","blastocyst")
plot_data$stage1<-factor(plot_data$stage1,levels=Stages)
plot_data$stage2<-factor(plot_data$stage2,levels=Stages)
head(plot_data)
plot_data[which(plot_data$variable == "chr1_11971658_11971857"),]
#plot_data$meth_mean<-plot_data$value*100
ggboxplot(plot_data, x = "stage1", y = "value",add = "jitter",size = 0.1)

## selectd DMRs corverage CpGs in every stage of early stages
plot_data_select<-plot_data[,c("sample","variable","value","stage1","stage2")]
plot_data_select2<-na.omit(plot_data_select)
dim(plot_data_select2)
head(plot_data_select2)
plot_data_select2$value<-1
plot_data_select3<-dcast(plot_data_select2[,c("variable","value","stage1")],variable~Stage)
rownames(plot_data_select3)<-plot_data_select3$variable
plot_data_select3<-plot_data_select3[,-1]
plot_data_select3[plot_data_select3<3]<-0
plot_data_select3[plot_data_select3>=3]<-1
head(plot_data_select3)

plot_data_select3$sum<-rowSums(plot_data_select3)
DMR_region<-rownames(plot_data_select3[which(plot_data_select3$sum>=4),])

plot_data_select2<-na.omit(plot_data_select)
plot_data_select3<-plot_data_select2[which(plot_data_select2$variable %in% DMR_region),]
head(plot_data_select3)
#plot_data_select3$value<-plot_data_select3$value*100
# summarySE 计算标准差和标准误差以及95%的置信区间.
tgc <- summarySE(plot_data_select3, measurevar="value", groupvars= c("stage1","variable"))
head(tgc);tail(tgc)
#https://blog.csdn.net/zhouhucheng00/article/details/88764082
embryo_plot<-ggplot() + 
  geom_errorbar(data=tgc,aes(x = stage1,ymin=(value-se), ymax=(value+se)),colour="black",size=.5,width=.2,position=position_dodge(.9)) +
  geom_point(data=tgc, aes(x = stage1,y = value),size = 1, shape = 21, fill = "red") + 
  geom_line(data=tgc, aes(x = stage1,y = value,colour= variable,group=variable)) + 
  stat_boxplot(geom = "errorbar",width=0.2)+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(data=plot_data_select3, aes(x = stage1,y = value,fill=variable),alpha = .5,width=0.5,size=0.2) +
  geom_jitter(width =0.2,shape = 21,size=0.1)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
  # scale_fill_manual(values = c("#E69F00", "#0072B2","#F0E442"))+  #设置填充的颜色
  #  scale_color_manual(values=c("black","black","black"))+ #设置散点图的圆圈的颜色为黑色
  theme_bw()+ ggtitle("methylation level for each DMRs")+ # ylab("Miles Per Gallon")+xlab("Number of Cylinders")+
  theme(legend.position="none", #不需要图例
        plot.title = element_text(hjust=0.5,size=5,vjust=0.5),
        axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),
        axis.line = element_line(colour="black"))+
  theme( axis.title.y = element_text(size=10,colour = "black",face = "bold"),
         axis.text.x = element_text(size=10),
         axis.text.y = element_text(size=10,colour = "black"))

embryo_plot2<-embryo_plot+facet_wrap(.~variable,ncol = 4)
ggsave("/mnt/data/chenwei/intergeneration_epigenetics/early_embryo/SM_box_meth_in_each_DMRs_for_early_embryo_500bp_extend.pdf",embryo_plot2,width=12, height=40)


##plot global DNA methylation of each sample in every stages
plot_data2 <-aggregate(x=plot_data$value,by=list(plot_data$sample,plot_data$stage1),FUN=mean, na.rm=TRUE, na.action=NULL)
head(plot_data2)
colnames(plot_data2)<-c("sample","stage","meth_mean")
#plot_data2$meth_mean<-plot_data2$meth_mean*100
ggboxplot(plot_data2, x = "stage", y = "meth_mean",color = "stage", add = c("mean"),add.params = list(color = c("red")),size = 0.5)
ggboxplot(plot_data2, x = "stage", y = "meth_mean",color = "stage", add = c("jitter"),add.params = list(color = c("grey")),size = 0.5)

plot_data2 <-aggregate(x=plot_data$value,by=list(plot_data$sample,plot_data$stage),FUN=mean, na.rm=TRUE, na.action=NULL)
head(plot_data2)
colnames(plot_data2)<-c("sample","stage","meth_mean")
#plot_data2$meth_mean<-plot_data2$meth_mean*100
ggboxplot(plot_data2, x = "stage", y = "meth_mean",add = "jitter",size = 0.5)+ylim(0,100)+
  scale_y_continuous(breaks = pretty_breaks(n = 10))


#分组求均值
head(plot_data_select3)
#趋势分析(STC, Series Test of Cluster)#ref:https://www.jianshu.com/p/ba78fbd05f3f
#https://www.jianshu.com/p/2067864259c5
library(TCseq)#1）调整聚类方法 2）调整类别数目
plot_data <-aggregate(x=plot_data_select3$value,by=list(plot_data_select3$variable,plot_data_select3$stage1),FUN=mean, na.rm=TRUE, na.action=NULL)
head(plot_data)
colnames(plot_data)<-c("variable","stage1","value")
plot_data<-dcast(plot_data,variable~stage1)
rownames(plot_data)<-plot_data$variable
plot_data<-plot_data[,-(1)]
head(plot_data)

#构建对象
data_expr3a<-as.matrix(plot_data)
set.seed(19921010)
cluster_num <- 6
tcseq_cluster_cm <- timeclust(data_expr3a, algo = 'cm', k = cluster_num, standardize = TRUE)
#timeclust()是一个整合函数，可执行数据标准化、聚类等多步操作，将上述输入数据中具有相似的时间表达特征的蛋白聚在一类。#standardize 用于 z-score 标准化变量
#基于模糊c均值聚类（timeclust()参数algo='cm'）的原理对蛋白质表达值的时间序列进行了聚类。timeclust()还提供了其它的聚类算法，如层次聚类（参数algo='hc'）、k均值划分（参数algo='km'）、围绕中心点划分（参数algo='pam'）等。
#颜色、线宽、坐标轴、字体等细节可以在函数中调整，具体参数详见函数帮助
timecplot_cm <- timeclustplot(tcseq_cluster_cm, value = 'z-score', cols = 3, 
                              axis.line.size = 0.6, axis.title.size = 8, axis.text.size = 8, 
                              title.size = 8, legend.title.size = 8, legend.text.size = 8)


tcseq_cluster_hc <- timeclust(data_expr3a, algo = 'hc', k = cluster_num, standardize = TRUE)
#timeclust()是一个整合函数，可执行数据标准化、聚类等多步操作，将上述输入数据中具有相似的时间表达特征的蛋白聚在一类。#standardize 用于 z-score 标准化变量
#基于模糊c均值聚类（timeclust()参数algo='cm'）的原理对蛋白质表达值的时间序列进行了聚类。timeclust()还提供了其它的聚类算法，如层次聚类（参数algo='hc'）、k均值划分（参数algo='km'）、围绕中心点划分（参数algo='pam'）等。
#颜色、线宽、坐标轴、字体等细节可以在函数中调整，具体参数详见函数帮助
timecplot_hc <- timeclustplot(tcseq_cluster_hc, value = 'z-score', cols = 3,
                              axis.line.size = 0.6, axis.title.size = 8, axis.text.size = 8, 
                              title.size = 8, legend.title.size = 8, legend.text.size = 8)

#上述获得了 10 组聚类群 #如果绘制单个的聚类群
timecplot_cm[2]

tcseq_cluster<-tcseq_cluster_cm
#查看每个蛋白所属的聚类群，展示前几个为例
head(tcseq_cluster@cluster)
#统计每个聚类群中各自包含的蛋白数量
table(tcseq_cluster@cluster)

#聚类过程通过计算 membership 值判断蛋白质所属的聚类群，以最大的 membership 值为准
head(tcseq_cluster@membership)

#标准化后的表达值(绘制曲线图值）
head(tcseq_cluster@data)

#最后，提取所有蛋白所属的聚类群，并和它们的原始表达值整合在一起
express_cluster <- tcseq_cluster@cluster
data_expr3b <- cbind(data_expr3a[names(express_cluster), ], express_cluster)
head(data_expr3b);dim(data_expr3b)
write.table(data_expr3b,'/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents.pdf', sep = '\t', col.names = NA, quote = FALSE)