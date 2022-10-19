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

##for down_stream
#extend 200 bins to 500 bins
FF_bed$start_raw<-FF_bed$start; FF_bed$end_raw <- FF_bed$end

FF_bed$start<-FF_bed$start-150;FF_bed$end<-FF_bed$end + 150
FF_bed$V4 <- as.character(paste(FF_bed$chr,FF_bed$start,FF_bed$end,sep="_"))

head(FF_bed);dim(FF_bed)
str(FF_bed)

#step3 reading single cell methylation files
JTA_files<-list.files("/mnt/data/anjianting/project/devomics/human/oocyte/methy/WCG.bed")
length(JTA_files)
JTA_files2<-as.list(paste0("/mnt/data/anjianting/project/devomics/human/oocyte/methy/WCG.bed/",JTA_files))
JTA_files2[[1]]
sample_names<-c();interSE<-list()
for (file in JTA_files2) {
  
  #file<-"/mnt/data/anjianting/project/devomics/human/oocyte/methy/WCG.bed/FGO_sc1.ACG.TCG.bed.gz"                                                
  #print(file)
  Sample_Name<-unlist(lapply(strsplit(file,"/"), function(x) x[11]))
  # print(as.character(Sample_Name))
  JTA_names<- as.character(lapply(strsplit(Sample_Name,"[.]"), function(x) x[1]))
  print(as.character(JTA_names))
  
  #data_meth1<-read.table(gzfile(file),header = F,fill=TRUE, na.strings = "",stringsAsFactors=F)
  data_meth1 = fread(file,header = F,fill=TRUE, na.strings = "",stringsAsFactors=F)
  head(data_meth1);dim(data_meth1)
  data_meth1$meth<-data_meth1$V4/(data_meth1$V4+data_meth1$V5)
  data_meth1<-data_meth1[,c("V1","V2","meth")]
  colnames(data_meth1)<-c("chr","start","meth")
  data_meth1$end<-data_meth1$start
  # head(data_meth1)
  data_meth2<-data_meth1[,.(chr,start,end,meth)]
  #data_meth2<-as.data.table(data_meth2)
  data_meth2<-na.omit(data_meth2)
  # str(data_meth2)
  data_meth2$start<-as.numeric(data_meth2$start)
  data_meth2$end<-as.numeric(data_meth2$end)
  data_meth2$meth<-as.numeric(data_meth2$meth)
  # head(FF_bed)
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
  # head(data_meth4)
  if(nrow(data_meth4)==0){next;}
  data_meth5 <-aggregate(data_meth4[,6],list(data_meth4[,7]),FUN=mean, na.rm=TRUE, na.action=NULL)
  colnames(data_meth5)<-c("region","meth")
  interSE<-c(interSE,list(as.data.frame(data_meth5)))
  sample_names<-c(sample_names,JTA_names)
  data_meth1<-data_meth2<-data_meth3<-data_meth4<-data_meth5<-0
} 

length(interSE);length(sample_names)
names(interSE)<- sample_names
saveRDS(interSE,file = "/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/JTA_single_Cell_human_oocyte_development_meth_data_in_two_corDMRs_extend_500bp.rds")

mymergeddata<-Reduce(function(x,y) {merge(x,y,all=T,by="region")}, interSE)
colnames(mymergeddata)<-c("region",sample_names)
head(mymergeddata);dim(mymergeddata)#96 364
saveRDS(mymergeddata,file = "/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/data_merge_JTA_single_Cell_human_oocyte_development_meth_data_in_two_corDMRs_extend_500bp.rds")
#rm(list = ls())
mymergeddata<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/data_merge_JTA_single_Cell_human_oocyte_development_meth_data_in_two_corDMRs_extend_500bp.rds")
#which(colnames(mymergeddata)=="Reference_CN")
colnames(mymergeddata)
rownames(mymergeddata)<-mymergeddata$region
mymergeddata<-mymergeddata[,-1]

rownames(mymergeddata) ## "chr22_20988051_20988550"

#show methlation level in heatmap in global
JTA_meta <-data.frame(sample=sample_names)
JTA_meta$sample<-as.character(JTA_meta$sample)
head(JTA_meta);dim(JTA_meta)
rownames(JTA_meta)<-JTA_meta$sample
JTA_meta$stage<- as.character(lapply(strsplit(JTA_meta$sample,"_"), function(x) x[1]))
table(JTA_meta$stage)
#FGO GO1 GO2  MI MII 
# 70  40  43 135  75 

#colData构建
colData <-JTA_meta[colnames(mymergeddata),]
dim(colData);length(colnames(mymergeddata))# 363
#colData<-colData[-c(which(rownames(colData)=="Reference_CN")),]
head(colData)
table(colData$stage)
#FGO GO1 GO2  MI MII 
# 70  40  43 135  75 
colData$stage<-factor(colData$stage,levels = c("GO1","GO2","FGO","MI","MII"))
colData<-colData[order(colData$stage,decreasing = T),]
mymergeddata<-mymergeddata[,rownames(colData)]

head(mymergeddata)
range(mymergeddata)

mymergeddata<-mymergeddata*100
dim(mymergeddata);dim(colData)
pheatmap(mymergeddata, cluster_rows=F, cluster_cols =F,scale="row")
pheatmap(mymergeddata, cluster_rows=F, cluster_cols =F,scale="none")

##去除导致无法聚类的行
#No rows with all NAs or zero variance like you said, but if you do the dist calculation, 
#there are NAs in some of the entries, indicating between some rows, it's not possible to calculate euclidean distances. 
#You need to the euclidean distance matrix to have no NAs to do clustering:
sum(is.na(as.matrix(dist(mymergeddata))))#4238
giveNAs = which(is.na(as.matrix(dist(mymergeddata))),arr.ind=TRUE)
#          row col

data_expr2[c(83,272),1:4]
##We get the rows out and start checking what to remove:
tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
  sum(is.na(as.matrix(dist(mymergeddata[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]

mymergeddata2 = mymergeddata[-as.numeric(rmv),]

head(mymergeddata2)
dim(mymergeddata2);dim(mymergeddata)#24 363 #  96 363
pheatmap(mymergeddata2, cluster_rows=T, cluster_cols =F,scale="row")
pheatmap(mymergeddata2, cluster_rows=T, cluster_cols =F,scale="none")

##正式绘制
##for no more remove uncluster region
plot_sample<-data.frame(colMeans(mymergeddata,na.rm = T))
plot_sample$sample<-rownames(plot_sample);colnames(plot_sample)<-c("meth","sample")
plot_sample2<-merge(plot_sample,colData,by=0)
head(plot_sample2)
plot_sample2<-plot_sample2[order(plot_sample2$stage,decreasing = F),]
table(plot_sample2$stage)
plot_sample2$sample<-as.character(plot_sample2$sample.x)
plot_sample2$sample<-factor(plot_sample2$sample,levels = as.character(plot_sample2$sample))

ggplot(data = plot_sample2, mapping = aes(x = sample, y = meth, fill = stage))+ 
  geom_bar(stat ='identity', position = 'dodge')+scale_color_manual(values=c(ppCor))+NoLegend()+
  theme_classic()+labs(x="stage",y="DNA meth level",title="Early embryo DNA methlation level in cor_DMRs")+
  theme(axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1,angle = 90))


colData2<-colData[,c("stage","stage")];colnames(colData2)<-c("stage1","stage2")
pheatmap(mymergeddata, cluster_rows=F, cluster_cols =F, annotation_col=colData2)

Stages<-c("GO1","GO2","FGO","MI","MII")
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
pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/JTA_single_Cell_human_oocyte_DNA_methylation_level_of_DMRs_in_heatmap_500bp_extend_1.pdf",width=16,height=8)
print(p1)
dev.off()


##for  more remove uncluster region
plot_sample<-data.frame(colMeans(mymergeddata2,na.rm = T))
plot_sample$sample<-rownames(plot_sample);colnames(plot_sample)<-c("meth","sample")
head(plot_sample)
plot_sample<-na.omit(plot_sample)
unique(plot_sample$sample)
colData1<-colData[unique(plot_sample$sample),]
plot_sample2<-merge(plot_sample,colData1,by=0)
head(plot_sample2)
plot_sample2<-plot_sample2[order(plot_sample2$stage,decreasing = F),]
table(plot_sample2$stage)
#GO1 GO2 FGO  MI MII 
# 36  41  61 107  67 
plot_sample2$sample<-as.character(plot_sample2$sample.x)
plot_sample2$sample<-factor(plot_sample2$sample,levels = as.character(plot_sample2$sample))

ggplot(data = plot_sample2, mapping = aes(x = sample, y = meth, fill = stage))+ 
  geom_bar(stat ='identity', position = 'dodge')+scale_color_manual(values=c(ppCor))+NoLegend()+
  theme_classic()+labs(x="stage",y="DNA meth level",title="Early embryo DNA methlation level in cor_DMRs")+
  theme(axis.text.x = element_text(size = 10,colour = 'black',vjust=1,hjust=1,angle = 90))


colData2<-colData1[,c("stage","stage")];colnames(colData2)<-c("stage1","stage2")

Stages<-c("GO1","GO2","FGO","MI","MII")
colData2$stage1<-factor(colData2$stage1,levels=Stages)
colData2<-colData2[order(colData2$stage1,decreasing = F),]
mymergeddata_3<-mymergeddata2[,rownames(colData2)]

p1<-pheatmap(mymergeddata_3, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
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
pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/JTA_single_Cell_human_oocyte_after_rm_DNA_methylation_level_of_DMRs_in_heatmap_500bp_extend_1.pdf",width=16,height=8)
print(p1)
dev.off()
#scale by row
p2<-pheatmap(mymergeddata_3, cluster_row =T,cluster_col =FALSE,na_col = "gray",
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

pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/JTA_single_Cell_human_oocyte_after_rm_DNA_methylation_level_of_DMRs_in_heatmap_500bp_extend_2.pdf",width=16,height=8)
print(p2)
dev.off()
dev.new()

head(mymergeddata);dim(mymergeddata)#101 221
number_row<-nrow(mymergeddata)
na_number<-data.frame(apply(mymergeddata, 2, function(x) sum(is.na(x))))
colnames(na_number)<-"freq"
na_number$freq2<-na_number$freq

table(as.numeric(apply(mymergeddata, 2, function(x) sum(is.na(x)))))
ggplot(na_number, aes(x = freq)) + geom_histogram(binwidth =2, fill = "lightblue", colour = "black")

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
dim(my_data3)# 29952     5
write.csv(my_data3,"/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/data_merge_JTA_single_Cell_human_oocyte_meth_data_in_two_corDMRs_add_meta_500bp_extend.csv")

##step five:  plot target regions DNA methylation
##for all region
plot_data <-read.csv("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/data_merge_JTA_single_Cell_human_oocyte_meth_data_in_two_corDMRs_add_meta_500bp_extend.csv",row.names = 1)
head(plot_data,n=7)
names(table(plot_data$stage2))

Stages<-c("GO1","GO2","FGO","MI","MII")
plot_data$stage1<-factor(plot_data$stage1,levels=Stages)
plot_data$stage2<-factor(plot_data$stage2,levels=Stages)
head(plot_data)
plot_data[which(plot_data$variable == "chr1_11971658_11971857"),]
#plot_data$meth_mean<-plot_data$value*100
ggboxplot(plot_data, x = "stage1", y = "value",add = "jitter",size = 0.1)


## selectd DMRs corverage CpGs in every stage of early stages
plot_data_select<-plot_data[,c("sample","variable","value","stage")]
plot_data_select2<-na.omit(plot_data_select)
dim(plot_data_select2)
head(plot_data_select2)
plot_data_select2$value<-1
plot_data_select3<-dcast(plot_data_select2[,c("variable","value","stage")],variable~stage)
rownames(plot_data_select3)<-plot_data_select3$variable
plot_data_select3<-plot_data_select3[,-1]

##修饰：对于特定区域，要求在该发育阶段至少有三个样本才对目标阶段进行平均值计算
plot_data_select3[plot_data_select3<3]<-0
plot_data_select3[plot_data_select3>=3]<-1
head(plot_data_select3)

##提取至少在一个阶段中存在甲基化值的区域
plot_data_select3$sum<-rowSums(plot_data_select3)
DMR_region<-rownames(plot_data_select3[which(plot_data_select3$sum>=1),])
length(DMR_region)#92

##重新抽取可以用于下游计算的数据
plot_data_select2<-na.omit(plot_data_select)
plot_data_select3<-plot_data_select2[which(plot_data_select2$variable %in% DMR_region),]
head(plot_data_select3)

###单一基因不同发育阶段表达动态绘制
#plot_data_select3$value<-plot_data_select3$value*100
# summarySE 计算标准差和标准误差以及95%的置信区间.
tgc <- summarySE(plot_data_select3, measurevar="value", groupvars= c("stage","variable"))
head(tgc);tail(tgc)
#https://blog.csdn.net/zhouhucheng00/article/details/88764082
embryo_plot<-ggplot() + 
  geom_errorbar(data=tgc,aes(x = stage,ymin=(value-se), ymax=(value+se)),colour="black",size=.5,width=.2,position=position_dodge(.9)) +
  geom_point(data=tgc, aes(x = stage,y = value),size = 1, shape = 21, fill = "red") + 
  geom_line(data=tgc, aes(x = stage,y = value,colour= variable,group=variable)) + 
  stat_boxplot(geom = "errorbar",width=0.2)+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(data=plot_data_select3, aes(x = stage,y = value,fill=variable),alpha = .5,width=0.5,size=0.2) +
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

embryo_plot2<-embryo_plot+facet_wrap(.~variable,ncol = 10)
embryo_plot2
ggsave("/mnt/data/chenwei/intergeneration_epigenetics/early_embryo/no_remove-box_meth_in_each_DMRs_for_JTA_single_Celloocyte_development_500bp_extend.pdf",embryo_plot,width=50, height=50,limitsize = FALSE)

##plot global DNA methylation of each sample in every stages
plot_data2 <-aggregate(x=plot_data_select3$value,by=list(plot_data_select3$sample,plot_data_select3$stage1),FUN=mean, na.rm=TRUE, na.action=NULL)
head(plot_data2)
colnames(plot_data2)<-c("sample","stage","meth_mean")
#plot_data2$meth_mean<-plot_data2$meth_mean*100
plot1<-ggboxplot(plot_data2, x = "stage", y = "meth_mean",color = "stage", add = c("mean"),add.params = list(color = c("red")),size = 0.5)
plot2<-ggboxplot(plot_data2, x = "stage", y = "meth_mean",color = "stage", add = c("jitter"),add.params = list(color = c("grey")),size = 0.5)

pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/JTA_single_Cell_oocyte_development_boxplot_DNA_methylation_level_of_DMRs_heatmap.pdf",width=16,height=8)
plot1
plot2
dev.off()

#选择区域的分组求均值后进行归一化 
head(plot_data_select3);dim(plot_data_select3)
plot_data2<-plot_data_select3[,c("variable","stage1","value")]
head(plot_data2)

plot_data3 <-aggregate(x=plot_data2$value,by=list(plot_data2$variable,plot_data2$stage1),FUN=mean, na.rm=TRUE, na.action=NULL)
head(plot_data3)
colnames(plot_data3)<-c("variable","stage1","value")
plot_data3<-dcast(plot_data3,variable~stage1)
rownames(plot_data3)<-plot_data3$variable
plot_data3<-plot_data3[,-(1)]
head(plot_data3)
stages<-c("GO1","GO2","FGO","MI","MII")
plot_data3<-plot_data3[,stages]
plot_data3
dim(plot_data3)

##归一化处理
plot_matrix <-data.frame(t(scale(t(plot_data3), center = TRUE, scale = TRUE)))
head(plot_matrix)
write.table(plot_matrix, '/mnt/data/chenwei/huahua/manuscript/figure/figure4/table/JTA_single_Celloocyte_development_all_cor_DMRs_zscore_data_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(plot_data3, '/mnt/data/chenwei/huahua/manuscript/figure/figure4/table/no_scale_JTA_single_Celloocyte_development_all_cor_DMRs_zscore_data_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents.txt', sep = '\t', col.names = NA, quote = FALSE)

##DMR根据发育阶段变化的聚类plot 
data_expr<-read.table('/mnt/data/chenwei/huahua/manuscript/figure/figure4/table/JTA_single_Celloocyte_development_all_cor_DMRs_zscore_data_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents.txt', header = T,row.names = 1)
head(data_expr)

##去除导致无法聚类的行
#No rows with all NAs or zero variance like you said, but if you do the dist calculation, 
#there are NAs in some of the entries, indicating between some rows, it's not possible to calculate euclidean distances. 
#You need to the euclidean distance matrix to have no NAs to do clustering:
sum(is.na(as.matrix(dist(data_expr))))# 134
giveNAs = which(is.na(as.matrix(dist(data_expr))),arr.ind=TRUE)
#          row col
##We get the rows out and start checking what to remove:
tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
  sum(is.na(as.matrix(dist(data_expr[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]

data_expr2 = data_expr[-as.numeric(rmv),]

head(data_expr2)
dim(data_expr2);dim(data_expr)# 67  5 # 68  5


#plot_matrix<-data_expr[,-ncol(data_expr)]
data_expr2[1:4,]
anno_cell<-c("GO1","GO2","FGO","MI","MII")
#anno_cell2<-c(rep("Zygote",1),rep("ZGA_pre",2),rep("ZGA_after",2),rep("Differention",2))
ann_colors = list(cell_stage=c(GO1=ppCor[1],GO2=ppCor[2],FGO=ppCor[3],MI=ppCor[4],MII=ppCor[5]))#,
#cell_stage2=c(Sperm=ppCor[1],Oocyte=ppCor[2],ZGA_pre=ppCor[3],ZGA_after=ppCor[4],Differention=ppCor[5]))
column_ha = HeatmapAnnotation(cell_stage = anno_cell,col = ann_colors)#cell_stage2 = anno_cell2,

htkm0 <- Heatmap(data_expr2,name= "z-score", border = TRUE,
                 top_annotation = column_ha,
                 cluster_column_slices = FALSE,column_title = "DNA methylation trend  in Human oocyte stages for correlated_DMRs_between_parent_kids", 
                 col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                 show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                 row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,
                 row_dend_width = unit(4, "cm")
)
print(htkm0)

pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/JTA_single_Celloocyte_development_all_cor_DMRs_zscore_data_gamete_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents_heatmaps.pdf",width = 6,height =6)
print(htkm0)
dev.off()

dend = as.dendrogram(hclust(dist(data_expr2),method="ward.D2"))
d_num<-5;dend = color_branches(dend, k =d_num)
htkm <- Heatmap(data_expr2,name= "z-score", border = TRUE,
                top_annotation = column_ha,
                cluster_column_slices = FALSE,column_title = "DNA methylation in Human early and early embryo stages for correlated_DMRs_between_parental_kids", 
                col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,#cluster_rows= FALSE,
                cluster_rows = dend, 
                row_split = d_num,row_gap = unit(rep(2, d_num-1), "mm"),
                row_dend_width = unit(4, "cm")
)
print(htkm)

pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/JTA_single_Celloocyte_development_split_mean_cluster_pattern_for_correlated_DMR_Kid_parents_heatmaps.pdf",width = 8,height =8)
print(htkm)
dev.off()

###step three :: plot trend line for gene expression of each cluster
#exact genes names in each cluster
clusterlist = row_order(htkm)
names(clusterlist)<-1:d_num
htkm_module <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(data_expr2[clusterlist[[i]],]),Cluster = paste0("cluster", i),stringsAsFactors = FALSE)
  return(out)}) %>%   do.call(rbind, .)

head(htkm_module);colnames(htkm_module)<-c("row_DMRs","Cluster")
table(htkm_module$Cluster)
#cluster1 cluster2 cluster3 cluster4 cluster5 
#     6       18       20       15        8 
head(data_expr2)
data_plot<-data_expr2[htkm_module$row_DMRs,]
cell_names<-colnames(data_plot)

data_plot$type<-as.character(htkm_module$Cluster)
data_plot$DMRs<-rownames(data_plot)

data_plot[1:4,]
data_preplot <- melt(data_plot,variable.name="cell",value.name = "expression",id.vars = c("type","DMRs"))
head(data_preplot)
data_preplot3<-data_preplot
range(data_preplot3$expression)# -2.297675  2.815053
data_preplot3<-data_preplot3[order(data_preplot3$cell),]
head(data_preplot3)
data_plot_A<-data_preplot3
#data_plot_B<-aggregate(expression ~ type+cell, data = data_preplot, mean)
data_plot_B<-aggregate(expression ~ type+cell, data = data_preplot, median)

data_plot_B<-data_plot_B[order(data_plot_B$cell),]
head(data_plot_A)
head(data_plot_B)
#colnames(data_plot_B)<-c("group","cell","expression")
trend_plot<-ggplot(data_plot_A,aes(x= cell,y=expression))+ 
  geom_point(aes(colour = cell),size=1,shape=19,alpha = 0.5)+
  geom_line(aes(group=DMRs), color="gray" ,position = position_dodge(0.02)) +
  geom_line(data=data_plot_B,aes(group=type),alpha=0.7,colour="steelblue",size=1.2,linetype=1)+
  scale_colour_manual(values=ppCor)+#ylim(-2,7)+
  theme_bw()+ theme(plot.title = element_text(hjust = 0.5),axis.text.y=element_blank(),
                    axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),
                    panel.grid=element_blank(),axis.ticks = element_blank())+
  xlab("pseudotime") + ylab("gene expression") 
trend_plot2<-trend_plot+ facet_wrap(~type,scales="free_y",ncol = 2) + stat_summary(fun = "median", color = "red", size = 1, geom = "point",alpha=0.7)
trend_plot2
trend_plot3<-trend_plot+ facet_wrap(~type,scales="free_y",ncol = 1) + stat_summary(fun = "median", color = "red", size = 1, geom = "point",alpha=0.7)

ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure4/JTA_single_Celloocyte_development_median_methy_trends_for_different_cluster_gamete_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents.pdf",trend_plot2,width =6, height =6,limitsize = FALSE)
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure4/JTA_single_Celloocyte_development_median_methy_trends_for_different_cluster_gamete_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents2.pdf",trend_plot3,width =4, height =15,limitsize = FALSE)


###为热图添加注释
data_plot[1:4,];dim(data_plot)
head(FF_bed)
FF_bed$V5<-paste(FF_bed$chr,FF_bed$start_raw,FF_bed$end_raw,sep = "_")
#DMR_anno$hg38_region<-paste(DMR_anno$chr,DMR_anno$start,DMR_anno$end,sep="_")

common_DMRs<-Reduce(intersect,list(father_cor_DMRs,mother_cor_DMRs))
MF_DMRs<- mother_cor_DMRs[!(mother_cor_DMRs %in% common_DMRs)]
FF_DMRs<- father_cor_DMRs[!(father_cor_DMRs %in% common_DMRs)]

#extend 200 bins to 500 bins
FF_bed$parental_region<-"none"
FF_bed$parental_region<-ifelse(FF_bed$V5 %in% FF_DMRs,"father",ifelse(FF_bed$V5 %in% MF_DMRs,"mother",ifelse(FF_bed$V5 %in% common_DMRs,"Common","no_DMRs")))
table(FF_bed$parental_region)
#ommon father mother 
#  6     52     60 
head(FF_bed);dim(FF_bed)
DMR_anno<-FF_bed
dim(DMR_anno)#118   8
DMR_anno2<-DMR_anno[which(DMR_anno$V4 %in% rownames(data_plot)),]
head(DMR_anno2);dim(DMR_anno2)# 67  8

#DMR_anno3<-merge(DMR_anno2,FF_bed[,c("V4","hg38_region","parental_region")],by="V4")
rownames(DMR_anno2)<-DMR_anno2$V4
head(DMR_anno2);dim(DMR_anno2)
DMR_anno2<-DMR_anno2[rownames(data_plot),]

anno_cell<-c("GO1","GO2","FGO","MI","MII")
#anno_cell2<-c(rep("Zygote",1),rep("ZGA_pre",2),rep("ZGA_after",2),rep("Differention",2))
ann_colors = list(cell_stage=c(GO1=ppCor[1],GO2=ppCor[2],FGO=ppCor[3],MI=ppCor[4],MII=ppCor[5]))#,
#cell_stage2=c(Sperm=ppCor[1],Oocyte=ppCor[2],ZGA_pre=ppCor[3],ZGA_after=ppCor[4],Differention=ppCor[5]))
column_ha = HeatmapAnnotation(cell_stage = anno_cell,col = ann_colors)#cell_stage2 = anno_cell2,

parental_origin<-as.character(DMR_anno2$parental_region)
ann_colors2 = list(parental_origin=c(mother=ppCor[2],father=ppCor[3],Common=ppCor[4]))
row_ha = rowAnnotation(origin = parental_origin,col = ann_colors2)

#row_ha = HeatmapAnnotation(parental_origin = parental_origin,col = ann_colors2)
groups<-names(table(data_plot$type))
gap_number0<-as.numeric(table(data_plot$type))
gap_number<-gap_number0[-length(gap_number0)]

data_plot2<-data_plot[,-c((ncol(data_plot)-1):ncol(data_plot))]#rownames(DMR_anno2)
htkm <- Heatmap(data_plot2,name= "z-score", border = TRUE,
                top_annotation = column_ha,left_annotation = row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                cluster_column_slices = FALSE,column_title = "DNA methylation in Human early and early embryo stages for correlated_DMRs_between_parental_kids", 
                col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,#cluster_rows= FALSE,
                # cluster_rows = dend, 
                row_dend_width = unit(4, "cm")
)
print(htkm)

#row_split = d_num,row_gap = unit(rep(2, d_num-1), "mm"),
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/JTA_single_Celloocyte_development_rearrange_split_mean_expression_cluster_pattern_for_correlated_DMR_Kid_parents_heatmaps.pdf",width = 8,height =8)
print(htkm)
dev.off()

##add gene to DMRs
#show genes for cor_DMRs
data_region_merge_DMR_all_big<-read.table("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/kids_AMA_vs_Young_200bp_DMR_myDiff15q005_merge_data.txt",header=T,sep="\t")
DMRs_unique<-unique(as.character(DMR_anno2$V5))
cor_DMR_gene<-data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$Row.names %in% DMRs_unique),c("bin_region","distanceToTSS","annotation","SYMBOL")]
write.table(cor_DMR_gene, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/JTA_single_Celloocyte_development_overlapped_cor_DMR_gene_kids_AMA_vs_Young_200bp_DMR_myDiff15q005.txt"),quote=F,row.names=F,col.names=T,sep = "\t")

head(cor_DMR_gene);dim(cor_DMR_gene)
cor_DMR_gene$V5<- cor_DMR_gene$bin_region
DMR_anno4<-merge(DMR_anno2,cor_DMR_gene[,c("distanceToTSS","SYMBOL","V5")],by="V5")
head(DMR_anno4)
rownames(DMR_anno4)<-DMR_anno4$V4

data_plot2<-data_plot[,-c((ncol(data_plot)-1):ncol(data_plot))]

plot_matrix2<-merge(data_plot2,DMR_anno4[,c("V4","SYMBOL")],by=0)
rownames(plot_matrix2)<-plot_matrix2$Row.names
plot_matrix2<-plot_matrix2[rownames(DMR_anno2),]
rownames(plot_matrix2)<-paste0(plot_matrix2$V4,"_",plot_matrix2$SYMBOL)
plot_matrix2<-plot_matrix2[,2:(ncol(plot_matrix2)-2)]
head(plot_matrix2);dim(plot_matrix2)#67  5

anno_cell<-c("GO1","GO2","FGO","MI","MII")
#anno_cell2<-c(rep("Zygote",1),rep("ZGA_pre",2),rep("ZGA_after",2),rep("Differention",2))
ann_colors = list(cell_stage=c(GO1=ppCor[1],GO2=ppCor[2],FGO=ppCor[3],MI=ppCor[4],MII=ppCor[5]))#,
#cell_stage2=c(Sperm=ppCor[1],Oocyte=ppCor[2],ZGA_pre=ppCor[3],ZGA_after=ppCor[4],Differention=ppCor[5]))
column_ha = HeatmapAnnotation(cell_stage = anno_cell,col = ann_colors)#cell_stage2 = anno_cell2,

parental_origin<-as.character(DMR_anno2$parental_region)
ann_colors2 = list(parental_origin=c(mother=ppCor[2],father=ppCor[3],Common=ppCor[4]))
row_ha = rowAnnotation(origin = parental_origin,col = ann_colors2)

#row_ha = HeatmapAnnotation(parental_origin = parental_origin,col = ann_colors2)
groups<-names(table(data_plot$type))
gap_number0<-as.numeric(table(data_plot$type))
gap_number<-gap_number0[-length(gap_number0)]

#plot_matrix<-plot_matrix[rownames(DMR_anno3),]
htkm <- Heatmap(plot_matrix2,name= "z-score", border = TRUE,
                top_annotation = column_ha,left_annotation = row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                cluster_column_slices = FALSE,column_title = "DNA methylation in Human oocyte stages for correlated_DMRs_between_parental_kids", 
                col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,cluster_rows= FALSE,
                # cluster_rows = dend, 
                row_dend_width = unit(4, "cm")
)#row_split = d_num,row_gap = unit(rep(2, d_num-1), "mm"),
print(htkm)

pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/Gene_add_rearrange_split_mean_JTA_single_Celloocyte_development_cluster_pattern_for_correlated_DMR_Kid_parents_heatmaps.pdf",width = 8,height =8)
print(htkm)
dev.off()

###使用原始表达矩阵
data_raw<-read.table('/mnt/data/chenwei/huahua/manuscript/figure/figure4/table/no_scale_JTA_single_Celloocyte_development_all_cor_DMRs_zscore_data_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents.txt', header = T,row.names = 1)
#colnames(data_raw)<-c("Zygote","Two_cell","Four_cell","Eight_cell","Morula","ICM","TE")
head(data_raw);head(DMR_anno4)
dim(DMR_anno4);dim(data_raw)

plot_matrix3<-merge(data_raw,DMR_anno4[,c("V4","SYMBOL")],by=0)
rownames(plot_matrix3)<-plot_matrix3$Row.names
plot_matrix3<-plot_matrix3[rownames(DMR_anno2),]
rownames(plot_matrix3)<-paste0(plot_matrix3$V4,"_",plot_matrix3$SYMBOL)
plot_matrix3<-plot_matrix3[,2:(ncol(plot_matrix3)-2)]
head(plot_matrix3);dim(plot_matrix3)# 67  5

anno_cell<-c("GO1","GO2","FGO","MI","MII")
#anno_cell2<-c(rep("Zygote",1),rep("ZGA_pre",2),rep("ZGA_after",2),rep("Differention",2))
ann_colors = list(cell_stage=c(GO1=ppCor[1],GO2=ppCor[2],FGO=ppCor[3],MI=ppCor[4],MII=ppCor[5]))#,
#cell_stage2=c(Sperm=ppCor[1],Oocyte=ppCor[2],ZGA_pre=ppCor[3],ZGA_after=ppCor[4],Differention=ppCor[5]))
column_ha = HeatmapAnnotation(cell_stage = anno_cell,col = ann_colors)#cell_stage2 = anno_cell2,

parental_origin<-as.character(DMR_anno2$parental_region)
ann_colors2 = list(parental_origin=c(mother=ppCor[2],father=ppCor[3],Common=ppCor[4]))
row_ha = rowAnnotation(origin = parental_origin,col = ann_colors2)

#row_ha = HeatmapAnnotation(parental_origin = parental_origin,col = ann_colors2)
groups<-names(table(data_plot$type))
gap_number0<-as.numeric(table(data_plot$type))
gap_number<-gap_number0[-length(gap_number0)]

#plot_matrix<-plot_matrix[rownames(DMR_anno3),]
htkm <- Heatmap(plot_matrix3,name= "DNA methylation level", border = TRUE,
                top_annotation = column_ha,left_annotation = row_ha,
                row_gap = unit(rep(2, length(gap_number)), "mm"), row_split =rep(groups,gap_number0),
                cluster_column_slices = FALSE,column_title = "DNA methylation in Human early and early embryo stages for correlated_DMRs_between_parental_kids", 
                #col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                col = colorRampPalette(c("navy","lightgoldenrod1","gold2","orange","darkorange","firebrick3"))(50),
                show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,cluster_rows= FALSE,
                # cluster_rows = dend, 
                row_dend_width = unit(4, "cm")
)#row_split = d_num,row_gap = unit(rep(2, d_num-1), "mm"),
print(htkm)

pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/raw_data_Gene_add_rearrange_split_mean_JTA_single_Celloocyte_development_cluster_pattern_for_correlated_DMR_Kid_parents_heatmaps.pdf",width = 8,height =8)
print(htkm)
dev.off()


#趋势分析(STC, Series Test of Cluster)#ref:https://www.jianshu.com/p/ba78fbd05f3f
#https://www.jianshu.com/p/2067864259c5
library(TCseq)#1）调整聚类方法 2）调整类别数目
#构建对象
data_expr3a<-as.matrix(plot_data3)
set.seed(19921010)
cluster_num <-5
tcseq_cluster_cm <- timeclust(data_expr3a, algo = 'cm', k = cluster_num, standardize = TRUE)
#timeclust()是一个整合函数，可执行数据标准化、聚类等多步操作，将上述输入数据中具有相似的时间表达特征的蛋白聚在一类。#standardize 用于 z-score 标准化变量
#基于模糊c均值聚类（timeclust()参数algo='cm'）的原理对蛋白质表达值的时间序列进行了聚类。timeclust()还提供了其它的聚类算法，如层次聚类（参数algo='hc'）、k均值划分（参数algo='km'）、围绕中心点划分（参数algo='pam'）等。
#颜色、线宽、坐标轴、字体等细节可以在函数中调整，具体参数详见函数帮助
timecplot_cm <- timeclustplot(tcseq_cluster_cm, value = 'z-score', cols = 3, 
                              axis.line.size = 0.6, axis.title.size = 8, axis.text.size = 8, 
                              title.size = 8, legend.title.size = 8, legend.text.size = 8)


tcseq_cluster_km <- timeclust(data_expr3a, algo = 'km', k = cluster_num, standardize = TRUE)
#timeclust()是一个整合函数，可执行数据标准化、聚类等多步操作，将上述输入数据中具有相似的时间表达特征的蛋白聚在一类。#standardize 用于 z-score 标准化变量
#基于模糊c均值聚类（timeclust()参数algo='cm'）的原理对蛋白质表达值的时间序列进行了聚类。timeclust()还提供了其它的聚类算法，如层次聚类（参数algo='hc'）、k均值划分（参数algo='km'）、围绕中心点划分（参数algo='pam'）等。
#颜色、线宽、坐标轴、字体等细节可以在函数中调整，具体参数详见函数帮助
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/timecplot_km_JTA_single_Celloocyte_development_expression_pattern_for_correlated_DMRs_Kid_parents.pdf", width =24,height = 12)
timecplot_km <- timeclustplot(tcseq_cluster_km, value = 'z-score', cols = 4,
                              axis.line.size = 0.6, axis.title.size = 8, axis.text.size = 8, 
                              title.size = 8, legend.title.size = 8, legend.text.size = 8)
dev.off()

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
write.table(data_expr3b,'/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/JTA_single_Celloocyte_development_cluster_pattern_for_correlated_DMR_Kid_parents.txt', sep = '\t', col.names = NA, quote = FALSE)

tcseq_cluster<-tcseq_cluster_km
express_cluster <- tcseq_cluster@cluster
zscore_data<-data.frame(tcseq_cluster@data)
data_expr3a <- cbind(zscore_data[names(express_cluster), ], express_cluster)
write.table(data_expr3a, '/mnt/data/chenwei/huahua/manuscript/figure/figure4/table/hc_km_zscore_data_JTA_single_Celloocyte_development_cluster_pattern_for_correlated_DMR_Kid_parents.txt', sep = '\t', col.names = NA, quote = FALSE)
