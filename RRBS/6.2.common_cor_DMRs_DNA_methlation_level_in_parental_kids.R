rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(ggplot2)
library(ggpubr) 
library(methylKit)
library(data.table)
library(scales)
library(ggsci)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(dendextend)
library(ComplexHeatmap)
library(circlize)
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
#read correlated DMRs
#for mother 
mother_cordata <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlated_pvalue_for_kid_AMA_DMR_between_kids_and_mother.txt",header=T,sep="\t")
MF_cor_DMR_sig<-as.character(mother_cordata[which(mother_cordata$group == "candidate_region"),]$DMR_region)
father_cordata <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlated_pvalue_for_kid_AMA_DMR_between_kids_and_father.txt",header=T,sep="\t")
FF_cor_DMR_sig<-as.character(father_cordata[which(father_cordata$group == "candidate_region"),]$DMR_region)
all_cor_DMRs<-unique(c(FF_cor_DMR_sig,MF_cor_DMR_sig))
length(all_cor_DMRs)#118

##read common DMRs
compare_name1<-"kids_AMA_vs_Young";split_region1<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region1,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper1<-unique(as.character(data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]$Row.names))
DMR_hypo1<-unique(as.character(data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]$Row.names))

compare_name2<-"mother_AMA_vs_Young";split_region2<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name2,"_",split_region2,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper3<-unique(as.character(data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]$Row.names))
DMR_hypo3<-unique(as.character(data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]$Row.names))

compare_name2<-"father_AMA_vs_Young";split_region2<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name2,"_",split_region2,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper4<-unique(as.character(data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]$Row.names))
DMR_hypo4<-unique(as.character(data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]$Row.names))

listinput<-c(list(DMR_hyper1),list(DMR_hyper3),list(DMR_hyper4),list(DMR_hypo1),list(DMR_hypo3),list(DMR_hypo4))
names(listinput)<-c("Kids_hyper","Mother_hyper","Father_hyper","Kids_hypo","Mother_hypo","Father_hypo" )

common_hyper<-Reduce(intersect,listinput[c("Kids_hyper","Mother_hyper","Father_hyper")])
common_hypo<-Reduce(intersect,listinput[c("Kids_hypo","Mother_hypo","Father_hypo")])
common_DMRs<-c(common_hyper,common_hypo)
commonlist_kid_mom_hyper<-Reduce(intersect,listinput[c("Kids_hyper","Mother_hyper")])
commonlist_kid_mom_hypo<-Reduce(intersect,listinput[c("Kids_hypo","Mother_hypo")])
mother_SF_DMRs_hyper<- commonlist_kid_mom_hyper[!(commonlist_kid_mom_hyper %in% common_hyper)]
mother_SF_DMRs_hypo<- commonlist_kid_mom_hypo[!(commonlist_kid_mom_hypo %in% common_hypo)]
MF_DMRs<-c(mother_SF_DMRs_hyper,mother_SF_DMRs_hypo)
commonlist_kid_dad_hyper<-Reduce(intersect,listinput[c("Kids_hyper","Father_hyper")])
commonlist_kid_dad_hypo<-Reduce(intersect,listinput[c("Kids_hypo","Father_hypo")])
father_SF_DMRs_hyper<- commonlist_kid_dad_hyper[!(commonlist_kid_dad_hyper %in% common_hyper)]

common_mix<-Reduce(intersect,listinput[c("Kids_hypo","Mother_hyper","Father_hypo")])
father_SF_DMRs_hypo<- commonlist_kid_dad_hypo[!(commonlist_kid_dad_hypo %in% c(common_hypo,common_mix))]
FF_DMRs<-c(father_SF_DMRs_hyper,father_SF_DMRs_hypo)

kid_mom_common_DMRs<-c(commonlist_kid_mom_hyper,commonlist_kid_mom_hypo)
kid_dad_common_DMRs<-c(commonlist_kid_dad_hyper,commonlist_kid_dad_hypo)
kid_parent_common_DMRs<-unique(c(kid_mom_common_DMRs,kid_dad_common_DMRs))

##filted commmon and correlated DMRs
commmon_cor_DMRs<-intersect(kid_parent_common_DMRs,all_cor_DMRs)
FF_bed<-as.data.frame(str_split_fixed(commmon_cor_DMRs, "_", n = 3))

colnames(FF_bed)<-c("chr","start","end")
FF_bed<-FF_bed[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",FF_bed$chr,invert=TRUE),]
FF_bed <- data.frame(FF_bed,V4=paste(FF_bed$chr,FF_bed$start,FF_bed$end,sep="_"))
FF_bed$chr<-factor(FF_bed$chr,levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
FF_bed<-FF_bed[order(FF_bed$chr,decreasing = F),]
FF_bed$chr<-as.character(FF_bed$chr);FF_bed$start<-as.numeric(as.character(FF_bed$start));FF_bed$end <-as.numeric(as.character(FF_bed$end ))
#extend 200 bins to 500 bins
FF_bed$start<-FF_bed$start-150;FF_bed$end<- FF_bed$end + 150
FF_bed$V4 <- as.character(paste(FF_bed$chr,FF_bed$start,FF_bed$end,sep="_"))
head(FF_bed);dim(FF_bed)
#get the hg19 bed file 
DMR_bed.hg19<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/DMR_bed_hg19_in_two_corDMRs_add_meta_500bp_extend.csv")
DMR_bed.hg19$hg19_region<-paste(DMR_bed.hg19$chr,DMR_bed.hg19$start,DMR_bed.hg19$end,sep="_")
head(DMR_bed.hg19)
DMR_bed.hg19$V4<-as.character(DMR_bed.hg19$V4)
DMR_bed.hg19_2<-DMR_bed.hg19[which(DMR_bed.hg19$V4 %in% FF_bed$V4),]
dim(DMR_bed.hg19_2)


mymergeddata<-read.table('/mnt/data/chenwei/huahua/manuscript/figure/figure4/table/na_rm_DNA_methylation_level_of_DMRs_in_early_embryo_development_for_correlated_DMR_Kid_parents.txt', header = T,check.names = F)
head(mymergeddata)
mymergeddata2<-mymergeddata[unique(DMR_bed.hg19_2$hg19_region),]
mymergeddata2[1:4,1:4]
dim(mymergeddata2)#42 178

##colData构建
ZP_meta <-read.csv(file="/mnt/data/chenwei/intergeneration_epigenetics/zhuping_sample_sc_meta.csv", header = T,row.names= 1)
ZP_meta$sample<-rownames(ZP_meta)
colData <-ZP_meta[colnames(mymergeddata2),]
dim(colData);length(colnames(mymergeddata2))
head(colData)
colData<-subset(colData,Stage != "Blastocyst (Mixed)")

colData2<-colData[,c("Stage","Developmental.Time")]
stages0<-c("Sperm","GV Oocyte","Mature Oocyte","early zygotic stage (10-11h after ICSI)", "mid- zygotic stage (22-23h after ICSI)","late zygotic stage (25h or even later after ICSI)","2-cell","4-cell","8-cell","Morula","Blastocyst")
colData2$Developmental.Time<-factor(colData2$Developmental.Time,levels=stages0)
stages<-c("Sperm","GV Oocyte","MII Oocyte","PN","2-cell","4-cell","8-cell","Morula","Blastocyst (ICM)","Blastocyst (TE)")
colData2$Stage<-factor(colData2$Stage,levels=stages)
colData2<-colData2[order(colData2$Stage,decreasing = F),]
colData2<-colData2[order(colData2$Developmental.Time,decreasing = F),]

mymergeddata2<-mymergeddata2[,rownames(colData2)]
mymergeddata2<-mymergeddata2*100
dim(mymergeddata2);dim(colData2)#42 178
pheatmap(mymergeddata2, cluster_rows=F, cluster_cols =F, annotation_col=colData2)


#step four : re-orgenization meth data 
my_data<-as.data.frame(t(mymergeddata2),stringsAsFactors = F)
#my_data<-my_data[-c(which(rownames(my_data)=="Reference_CN")),]
my_data$sample<-rownames(my_data)
my_data2 <- melt(my_data)
head(my_data2)

my_data3<-merge(my_data2,ZP_meta,by="sample")
head(my_data3)
head(my_data3,n=7)
dim(my_data3)#7476   11
write.csv(my_data3,"/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/common_cor_DMR_data_merge_ZP_NG_single_Cell_meth_data_in_two_corDMRs_add_meta_500bp_extend.csv")

##step five:  plot target regions DNA methylation
##for all region
plot_data <-read.csv("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/common_cor_DMR_data_merge_ZP_NG_single_Cell_meth_data_in_two_corDMRs_add_meta_500bp_extend.csv",row.names = 1)
head(plot_data,n=7)
names(table(plot_data$Developmental.Time))
names(table(plot_data$Stage))
stages0<-c("Sperm","GV Oocyte","Mature Oocyte","early zygotic stage (10-11h after ICSI)", "mid- zygotic stage (22-23h after ICSI)","late zygotic stage (25h or even later after ICSI)","2-cell","4-cell","8-cell","Morula","Blastocyst")
plot_data$Developmental.Time<-factor(plot_data$Developmental.Time,levels=stages0)
stages<-c("Sperm","GV Oocyte","MII Oocyte","PN","2-cell","4-cell","8-cell","Morula","Blastocyst (ICM)","Blastocyst (TE)")
plot_data$Stage<-factor(plot_data$Stage,levels=stages)
head(plot_data)
plot_data[which(plot_data$variable == "chr1_11971658_11971857"),]
#plot_data$meth_mean<-plot_data$value*100
ggboxplot(plot_data, x = "Stage", y = "value",add = "jitter",size = 0.1)

## selectd DMRs corverage CpGs in every stage of early stages
plot_data_select<-plot_data[,c("sample","variable","value","Stage","Developmental.Time")]
plot_data_select2<-na.omit(plot_data_select)
dim(plot_data_select2)
head(plot_data_select2)
plot_data_select2$value<-1
plot_data_select3<-dcast(plot_data_select2[,c("variable","value","Stage")],variable~Stage)
rownames(plot_data_select3)<-plot_data_select3$variable
plot_data_select3<-plot_data_select3[,-1]
plot_data_select3[plot_data_select3<3]<-0
plot_data_select3[plot_data_select3>=3]<-1
head(plot_data_select3)

plot_data_select3$sum<-rowSums(plot_data_select3)
DMR_region<-rownames(plot_data_select3[which(plot_data_select3$sum>=9),])

plot_data_select2<-na.omit(plot_data_select)
plot_data_select3<-plot_data_select2[which(plot_data_select2$variable %in% DMR_region),]
head(plot_data_select3)
#plot_data_select3$value<-plot_data_select3$value*100
# summarySE 计算标准差和标准误差以及95%的置信区间.
tgc <- summarySE(plot_data_select3, measurevar="value", groupvars= c("Stage","variable"))
head(tgc);tail(tgc)
#https://blog.csdn.net/zhouhucheng00/article/details/88764082
embryo_plot<-ggplot() + 
  geom_errorbar(data=tgc,aes(x = Stage,ymin=(value-se), ymax=(value+se)),colour="black",size=.5,width=.2,position=position_dodge(.9)) +
  geom_point(data=tgc, aes(x = Stage,y = value),size = 1, shape = 21, fill = "red") + 
  geom_line(data=tgc, aes(x = Stage,y = value,colour= variable,group=variable)) + 
  stat_boxplot(geom = "errorbar",width=0.2)+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
  geom_boxplot(data=plot_data_select3, aes(x = Stage,y = value,fill=variable),alpha = .5,width=0.5,size=0.2) +
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
ggsave("/mnt/data/chenwei/intergeneration_epigenetics/early_embryo/common_cor_DMRs_box_meth_in_each_DMRs_for_early_embryo_500bp_extend.pdf",embryo_plot,width=12, height=15)

##plot global DNA methylation of each sample in every stages

plot_data2 <-aggregate(x=plot_data$value,by=list(plot_data$sample,plot_data$Stage),FUN=mean, na.rm=TRUE, na.action=NULL)
head(plot_data2)
colnames(plot_data2)<-c("sample","stage","meth_mean")
#plot_data2$meth_mean<-plot_data2$meth_mean*100
plot1<-ggboxplot(plot_data2, x = "stage", y = "meth_mean",color = "stage", add = c("mean"),add.params = list(color = c("red")),size = 0.5)
plot2<-ggboxplot(plot_data2, x = "stage", y = "meth_mean",color = "stage", add = c("jitter"),add.params = list(color = c("grey")),size = 0.5)

plot_data2 <-aggregate(x=plot_data$value,by=list(plot_data$sample,plot_data$Developmental.Time),FUN=mean, na.rm=TRUE, na.action=NULL)
head(plot_data2)
colnames(plot_data2)<-c("sample","stage","meth_mean")
#plot_data2$meth_mean<-plot_data2$meth_mean*100
plot3<-ggboxplot(plot_data2, x = "stage", y = "meth_mean",add = "jitter",size = 0.5)

pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/common_cor_boxplot_DNA_methylation_level_of_DMRs_in_early_embryo_development_heatmap.pdf",width=16,height=8)
plot1
plot2
plot3
dev.off()

#分组求均值
head(plot_data)
plot_data2<-plot_data[,c("variable","Stage","Developmental.Time","value")]
head(plot_data2)
plot_data3 <-aggregate(x=plot_data2$value,by=list(plot_data2$variable,plot_data2$Stage),FUN=mean, na.rm=TRUE, na.action=NULL)
colnames(plot_data3)<-c("variable","stage1","value")
plot_data3<-dcast(plot_data3,variable~stage1)
rownames(plot_data3)<-plot_data3$variable
plot_data3<-plot_data3[,-(1)]
head(plot_data3)
stages<-c("Sperm","GV Oocyte","MII Oocyte","PN","2-cell","4-cell","8-cell","Morula","Blastocyst (ICM)","Blastocyst (TE)")
plot_data3<-plot_data3[,stages]
plot_matrix <-data.frame(t(scale(t(plot_data3), center = TRUE, scale = TRUE)))
plot_matrix
write.table(plot_matrix, '/mnt/data/chenwei/huahua/manuscript/figure/figure4/table/common_cor_DMRs_zscore_data_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents.txt', sep = '\t', col.names = NA, quote = FALSE)

##plot 
data_expr<-read.table('/mnt/data/chenwei/huahua/manuscript/figure/figure4/table/common_cor_DMRs_zscore_data_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents.txt', header = T,row.names = 1)
head(data_expr)
colnames(data_expr)<-c("Sperm","GV_Oocyte","MII_Oocyte","PN","Two_cell","Four_cell","Eight_cell","Morula","Blastocyst_ICM","Blastocyst_TE")
data_expr[1:4,]

anno_cell<-c("Sperm","GV_Oocyte","MII_Oocyte","PN","Two_cell","Four_cell","Eight_cell","Morula","Blastocyst_ICM","Blastocyst_TE")
anno_cell2<-c(rep("Sperm",1),rep("Oocyte",2),rep("ZGA_pre",3),rep("ZGA_after",2),rep("Differention",2))
ann_colors = list(cell_stage=c(Sperm=ppCor[1],GV_Oocyte=ppCor[2],MII_Oocyte=ppCor[3],PN=ppCor[4],Two_cell=ppCor[5],Four_cell=ppCor[6],Eight_cell=ppCor[7],
                               Morula=ppCor[8],Blastocyst_ICM=ppCor[9],Blastocyst_TE=ppCor[10]),
                  cell_stage2=c(Sperm=ppCor[1],Oocyte=ppCor[2],ZGA_pre=ppCor[3],ZGA_after=ppCor[4],Differention=ppCor[5]))
column_ha = HeatmapAnnotation(cell_stage = anno_cell,cell_stage2 = anno_cell2,col = ann_colors)

htkm0 <- Heatmap(data_expr,name= "z-score", border = TRUE,
                 top_annotation = column_ha,
                 cluster_column_slices = FALSE,column_title = "DNA methylation trend  in Human gamete and early embryo stages for correlated_DMRs_between_parent_kids", 
                 col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                 show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                 row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,
                 row_dend_width = unit(4, "cm")
)
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/common_cor_DMRs_zscore_data_gamete_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents_heatmaps.pdf",width = 6,height =6)
print(htkm0)
dev.off()

dend = as.dendrogram(hclust(dist(data_expr),method="ward.D2"))
d_num<-9;dend = color_branches(dend, k =d_num)
htkm <- Heatmap(data_expr,name= "z-score", border = TRUE,
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
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/split_mean_common_cor_DMRs_zscore_gamete_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents_heatmaps.pdf",width = 8,height =8)
print(htkm)
dev.off()

###step three :: plot trend line for gene expression of each cluster
#exact genes names in each cluster
clusterlist = row_order(htkm)
names(clusterlist)<-1:d_num
htkm_module <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(plot_matrix[clusterlist[[i]],]),Cluster = paste0("cluster", i),stringsAsFactors = FALSE)
  return(out)}) %>%   do.call(rbind, .)

head(htkm_module);colnames(htkm_module)<-c("row_DMRs","Cluster")
table(htkm_module$Cluster)
#cluster1 cluster2 cluster3 cluster4 cluster5 cluster6 cluster7 cluster8 cluster9 
#  5       10        5        4        7        4        3        1        3 
head(data_expr)
data_plot<-data.frame(data_expr[htkm_module$row_DMRs,])
data_plot$type<-as.character(htkm_module$Cluster)
cell_names<-colnames(data_plot)

data_plot[1:4,]
data_plot$DMRs<-rownames(data_plot)
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
  xlab("pseudotime") + ylab("gene expression") + facet_wrap(~type,scales="free_y",ncol = 1)
trend_plot2<-trend_plot + stat_summary(fun = "median", color = "red", size = 1, geom = "point",alpha=0.7)
trend_plot2
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure4/mean_methy_trends_for_different_cluster_common_cor_DMRs_zscore_gamete_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents.pdf",trend_plot2,width =4, height =14,limitsize = FALSE)
                                                                           
###为热图添加注释
data_plot[1:4,]
rownames(data_plot) %in% paste(DMR_bed.hg19$chr,DMR_bed.hg19$start,DMR_bed.hg19$end,sep="_")
DMR_anno<-DMR_bed.hg19
DMR_anno$hg19_region<-paste(DMR_anno$chr,DMR_anno$start,DMR_anno$end,sep="_")

#extend 200 bins to 500 bins
FF_bed$start_orig<-FF_bed$start+150;FF_bed$end_orig<- FF_bed$end-150
FF_bed$hg38_region <- as.character(paste(FF_bed$chr,FF_bed$start_orig,FF_bed$end_orig,sep="_"))
FF_bed$parental_region<-"none"
FF_bed$parental_region<-ifelse(FF_bed$hg38_region %in% FF_DMRs,"father",ifelse(FF_bed$hg38_region %in% MF_DMRs,"mother","Common"))
table(FF_bed$parental_region)
head(FF_bed);dim(FF_bed)

DMR_anno2<-DMR_anno[which(DMR_anno$hg19_region %in% rownames(data_plot)),]
head(DMR_anno2);dim(DMR_anno2)
DMR_anno3<-merge(DMR_anno2,FF_bed[,c("V4","hg38_region","parental_region")],by="V4")
rownames(DMR_anno3)<-DMR_anno3$hg19_region
DMR_anno3<-DMR_anno3[rownames(data_plot),]
head(DMR_anno3);dim(DMR_anno3)

anno_cell<-c("Sperm","GV_Oocyte","MII_Oocyte","PN","Two_cell","Four_cell","Eight_cell","Morula","Blastocyst_ICM","Blastocyst_TE")
anno_cell2<-c(rep("Sperm",1),rep("Oocyte",2),rep("ZGA_pre",3),rep("ZGA_after",2),rep("Differention",2))
ann_colors = list(cell_stage=c(Sperm=ppCor[1],GV_Oocyte=ppCor[2],MII_Oocyte=ppCor[3],PN=ppCor[4],Two_cell=ppCor[5],Four_cell=ppCor[6],Eight_cell=ppCor[7],
                               Morula=ppCor[8],Blastocyst_ICM=ppCor[9],Blastocyst_TE=ppCor[10]),
                  cell_stage2=c(Sperm=ppCor[1],Oocyte=ppCor[2],ZGA_pre=ppCor[3],ZGA_after=ppCor[4],Differention=ppCor[5]))
column_ha = HeatmapAnnotation(cell_stage = anno_cell,cell_stage2 = anno_cell2,col = ann_colors)
parental_origin<-as.character(DMR_anno3$parental_region)
ann_colors2 = list(parental_origin=c(mother=ppCor[1],father=ppCor[3],Common=ppCor[6]))
row_ha = rowAnnotation(origin = parental_origin,col = ann_colors2)

#row_ha = HeatmapAnnotation(parental_origin = parental_origin,col = ann_colors2)
groups<-names(table(data_plot$type))
gap_number0<-as.numeric(table(data_plot$type))
gap_number<-gap_number0[-length(gap_number0)]

data_plot2<-data_plot[rownames(DMR_anno3),-c((ncol(data_plot)-1):ncol(data_plot))]
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
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/commmon_cor_rearrange_split_mean_gamete_early_embryo_expression_time_cluster_pattern_for_correlated_DMR_Kid_parents_heatmaps1.pdf",width = 8,height =8)
print(htkm)
dev.off()

##add gene to DMRs
rownames(data_plot) %in% paste(DMR_bed.hg19$chr,DMR_bed.hg19$start,DMR_bed.hg19$end,sep="_")
DMR_anno<-DMR_bed.hg19
DMR_anno$hg19_region<-paste(DMR_anno$chr,DMR_anno$start,DMR_anno$end,sep="_")

#extend 200 bins to 500 bins
FF_bed$start_orig<-FF_bed$start+150;FF_bed$end_orig<- FF_bed$end-150
FF_bed$hg38_region <- as.character(paste(FF_bed$chr,FF_bed$start_orig,FF_bed$end_orig,sep="_"))
FF_bed$parental_region<-"none"
FF_bed$parental_region<-ifelse(FF_bed$hg38_region %in% FF_DMRs,"father",ifelse(FF_bed$hg38_region %in% MF_DMRs,"mother","Common"))
table(FF_bed$parental_region)
head(FF_bed);dim(FF_bed)

DMR_anno2<-DMR_anno[which(DMR_anno$hg19_region %in% rownames(data_plot)),]
head(DMR_anno2);dim(DMR_anno2)
DMR_anno3<-merge(DMR_anno2,FF_bed[,c("V4","hg38_region","parental_region")],by="V4")
rownames(DMR_anno3)<-DMR_anno3$hg19_region
DMR_anno3<-DMR_anno3[rownames(data_plot),]
head(DMR_anno3);dim(DMR_anno3)

#show genes for cor_DMRs
data_region_merge_DMR_all_big<-read.table("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/kids_AMA_vs_Young_200bp_DMR_myDiff15q005_merge_data.txt",header=T,sep="\t")
DMRs_unique<-unique(as.character(DMR_anno3$hg38_region))
cor_DMR_gene<-data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$Row.names %in% DMRs_unique),c("bin_region","distanceToTSS","annotation","SYMBOL")]
write.table(cor_DMR_gene, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/single_cell_overlapped_common_cor_DMR_gene_kids_AMA_vs_Young_200bp_DMR_myDiff15q005.txt"),quote=F,row.names=F,col.names=T,sep = "\t")

head(cor_DMR_gene)
cor_DMR_gene$hg38_region<- cor_DMR_gene$bin_region
DMR_anno4<-merge(DMR_anno3,cor_DMR_gene[,c("distanceToTSS","SYMBOL","hg38_region")],by="hg38_region")
head(DMR_anno4)
rownames(DMR_anno4)<-DMR_anno4$hg19_region

data_plot2<-data_plot[rownames(DMR_anno3),-c((ncol(data_plot)-1):ncol(data_plot))]
plot_matrix2<-merge(data_plot2,DMR_anno4[,c("hg19_region","SYMBOL")],by=0)
rownames(plot_matrix2)<-plot_matrix2$Row.names
plot_matrix2<-plot_matrix2[rownames(DMR_anno3),]
rownames(plot_matrix2)<-paste0(plot_matrix2$hg19_region,"_",plot_matrix2$SYMBOL)
plot_matrix2<-plot_matrix2[,2:(ncol(plot_matrix2)-2)]
head(plot_matrix2);dim(plot_matrix2)# 42  10

anno_cell<-c("Sperm","GV_Oocyte","MII_Oocyte","PN","Two_cell","Four_cell","Eight_cell","Morula","Blastocyst_ICM","Blastocyst_TE")
anno_cell2<-c(rep("Sperm",1),rep("Oocyte",2),rep("ZGA_pre",3),rep("ZGA_after",2),rep("Differention",2))
ann_colors = list(cell_stage=c(Sperm=ppCor[1],GV_Oocyte=ppCor[2],MII_Oocyte=ppCor[3],PN=ppCor[4],Two_cell=ppCor[5],Four_cell=ppCor[6],Eight_cell=ppCor[7],
                               Morula=ppCor[8],Blastocyst_ICM=ppCor[9],Blastocyst_TE=ppCor[10]),
                  cell_stage2=c(Sperm=ppCor[1],Oocyte=ppCor[2],ZGA_pre=ppCor[3],ZGA_after=ppCor[4],Differention=ppCor[5]))
column_ha = HeatmapAnnotation(cell_stage = anno_cell,cell_stage2 = anno_cell2,col = ann_colors)
parental_origin<-as.character(DMR_anno3$parental_region)
ann_colors2 = list(parental_origin=c(mother=ppCor[1],father=ppCor[3],Common=ppCor[6]))
row_ha = rowAnnotation(origin = parental_origin,col = ann_colors2)

#row_ha = HeatmapAnnotation(parental_origin = parental_origin,col = ann_colors2)
groups<-names(table(data_plot$type))
gap_number0<-as.numeric(table(data_plot$type))
gap_number<-gap_number0[-length(gap_number0)]
dim(plot_matrix2)#42 10
#plot_matrix<-plot_matrix[rownames(DMR_anno3),]
htkm <- Heatmap(plot_matrix2,name= "z-score", border = TRUE,
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
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/Gene_add_rearrange_split_mean_gamete_early_embryo_expression_time_cluster_pattern_for_common_correlated_DMR_Kid_parents_heatmaps2.pdf",width = 8,height =8)
print(htkm)
dev.off()
