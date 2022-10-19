rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ggpubr)
library(methylKit)
library(scatterplot3d)
library(pcaMethods)
library(scales)
library(ggsci)
library(grid)
library(reshape2)
library(stringr)
library(pheatmap)
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#############for father
#selection for DMRs FOR AMA vs Young
father_200bin_data<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv",row.names = 1,check.names=F)
head(father_200bin_data);dim(father_200bin_data)
#701937     26

##test evaluation
str(father_200bin_data)
myDiffq005_all<-father_200bin_data[father_200bin_data$qvalue <0.05,]
dim(myDiffq005_all)#211740     26
myDiffq005_all$group1<-as.factor(ifelse(myDiffq005_all$meth.diff>0, "hyper", "hypo"))
myDiffq005_all$group_mean<-as.factor(ifelse(myDiffq005_all$meth_mean_def>0, "hyper", "hypo"))
head(myDiffq005_all);dim(myDiffq005_all)

pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/father_AMA_vs_Young_200bin_q005_dif15_distribution.pdf")
# Histogram overlaid with kernel density curve ##https://www.jianshu.com/p/95f54a8b1269
ggplot(myDiffq005_all, aes(x=meth.diff)) +
  geom_histogram(aes(y=..density..,fill=group1), binwidth=1,colour="gray") + # 这一步很重要,使用density代替y轴
  #  geom_density(alpha=.2,aes(fill=group))+
  geom_density(alpha=.2, fill="#FF6666")+ # 重叠部分采用透明设置
  geom_vline(xintercept=c(-50,-25,-15,15,25,50),linetype="dashed",colour=c("blue","blue","blue","red","red","red"))+
  theme_bw()+scale_x_continuous(breaks=seq(-55, 55, 5))+
  theme(axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5,  hjust = 0.5))+
  theme(axis.title.y = element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5))+labs(title = "methykit difference: AMA vs Young")

ggplot(myDiffq005_all, aes(x=meth_mean_def)) +
  geom_histogram(aes(y=..density..,fill=group_mean), binwidth=1,colour="gray") + # 这一步很重要,使用density代替y轴
  #  geom_density(alpha=.2,aes(fill=group))+
  geom_density(alpha=.2, fill="#FF6666")+ # 重叠部分采用透明设置
  geom_vline(xintercept=c(-50,-25,-15,15,25,50),linetype="dashed",colour=c("blue","blue","blue","red","red","red"))+
  theme_bw()+scale_x_continuous(breaks=seq(-55, 55, 5))+
  theme(axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5,  hjust = 0.5))+
  theme(axis.title.y = element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5))+labs(title = "mean difference: : AMA vs Young")
dev.off()

summary(myDiffq005_all$meth.diff)
#      Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#    -32.019  -3.135  -1.753  -1.251   1.073  32.958 
 
summary(myDiffq005_all[myDiffq005_all$group1=="hyper",]$meth.diff)
#      Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.1487  1.2391  2.0257  2.6340  3.3259 32.9581

summary(myDiffq005_all[myDiffq005_all$group1=="hypo",]$meth.diff)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-32.019  -3.788  -2.539  -2.998  -1.669  -0.137

library(psych)
describe(myDiffq005_all$meth.diff)
##非缺失值的数量、平均数、标准差、中位数、截尾均值、绝对中位差、最小值、最大值、值域、偏度、峰度和平均值的标准误
##vars     n    mean    sd median trimmed   mad    min   max  range skew kurtosis   se
##X1    1 211740 -1.25 3.32  -1.75   -1.35 2.77 -32.02 32.96 64.98 0.37     3.38 0.01

describe(myDiffq005_all[myDiffq005_all$group1=="hyper",]$meth.diff)
#vars     n     mean    sd median trimmed  mad  min   max range skew kurtosis   se
##X1    1 65675 2.63 2.18   2.03    2.27 1.4 0.15 32.96 32.81  2.8    14.29 0.01

describe(myDiffq005_all[myDiffq005_all$group1=="hypo",]$meth.diff)
#vars        n   mean   sd median trimmed  mad    min   max range  skew kurtosis   se
##X1    1 146065   -3  2  -2.54   -2.72 1.48 -32.02 -0.14 31.88 -2.35     11.6 0.01

pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/father_AMA_vs_Young_Diff_200bp_DMR_pie50_25_20_10_5_def_plot.pdf")
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=5,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=10,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=15,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=20,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=25,]$group1))
dev.off()

#pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/AMA_arrest_vs_AMA_abortion_q005_Diff_200bp_DMR_vocano_plot.pdf")
#https://www.jianshu.com/p/41375176906e

father_200bin_data$threshold = as.factor(ifelse(father_200bin_data$qvalue < 0.05 & abs(father_200bin_data$meth.diff)>= 50, ifelse(father_200bin_data$meth.diff>= 50 ,'Hyper','Hypo'),'NoSignifi'))
table(father_200bin_data$threshold)
#NoSignifi 
#  701937    
father_200bin_data$threshold = as.factor(ifelse(father_200bin_data$qvalue < 0.05 & abs(father_200bin_data$meth.diff)>= 25, ifelse(father_200bin_data$meth.diff>= 25 ,'Hyper','Hypo'),'NoSignifi'))
table(father_200bin_data$threshold)
# Hyper      Hypo NoSignifi 
#  12        17    701908 
father_200bin_data$threshold = as.factor(ifelse(father_200bin_data$qvalue < 0.05 & abs(father_200bin_data$meth.diff)>= 20, ifelse(father_200bin_data$meth.diff>= 20 ,'Hyper','Hypo'),'NoSignifi'))
table(father_200bin_data$threshold)
#Hyper      Hypo NoSignifi 
#45        67    701825
father_200bin_data$threshold = as.factor(ifelse(father_200bin_data$qvalue < 0.05 & abs(father_200bin_data$meth.diff)>= 15, ifelse(father_200bin_data$meth.diff>= 15 ,'Hyper','Hypo'),'NoSignifi'))
table(father_200bin_data$threshold)
#  Hyper    Hypo  NoSignifi 
#  206       258    701473 
father_200bin_data$threshold = as.factor(ifelse(father_200bin_data$qvalue < 0.05 & abs(father_200bin_data$meth.diff)>= 10, ifelse(father_200bin_data$meth.diff>= 10 ,'Hyper','Hypo'),'NoSignifi'))
table(father_200bin_data$threshold)
#Hyper      Hypo    NoSignifi 
#906      1575    699456 

#Final usage:
compare_name<-"father_AMA_vs_Young"
father_200bin_data$threshold = as.factor(ifelse(father_200bin_data$qvalue < 0.05 & abs(father_200bin_data$meth.diff)>= 25, ifelse(father_200bin_data$meth.diff>= 25 ,'Hyper','Hypo'),'NoSignifi'))
table(father_200bin_data$threshold)
father_200bin_data$threshold<-factor(father_200bin_data$threshold,level=c("Hyper","Hypo","NoSignifi"))
p_df25<-ggplot(data = father_200bin_data, aes(x = meth.diff, y = -log10(qvalue), colour=threshold)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("red","blue","grey")) +
  scale_x_continuous(breaks=seq(-80, 80, 20))+
  ylim(0, 600)+
  geom_vline(xintercept=c(-25,25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="difference in DNA methlation level ",y="-log10 (q-value)",title=paste0("methykit difference: ",compare_name,"_Diff25q005_DMR")) +
  theme_bw()+theme(panel.border = element_rect(colour="black",fill=NA),panel.grid.major =element_blank())+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5), legend.position="right", legend.title = element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))

father_200bin_data$threshold = as.factor(ifelse(father_200bin_data$qvalue < 0.05 & abs(father_200bin_data$meth.diff)>= 20, ifelse(father_200bin_data$meth.diff>= 20 ,'Hyper','Hypo'),'NoSignifi'))
table(father_200bin_data$threshold)
father_200bin_data$threshold<-factor(father_200bin_data$threshold,level=c("Hyper","Hypo","NoSignifi"))
p_df20<-ggplot(data = father_200bin_data, aes(x = meth.diff, y = -log10(qvalue), colour=threshold)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("red","blue","grey")) +
  scale_x_continuous(breaks=seq(-80, 80, 20))+
  ylim(0, 600)+
  geom_vline(xintercept=c(-20,20),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="difference in DNA methlation level ",y="-log10 (q-value)",title=paste0("methykit difference: ",compare_name,"_Diff20q005_DMR")) +
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),panel.grid.major =element_blank())+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5), legend.position="right", legend.title = element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))

father_200bin_data$threshold = as.factor(ifelse(father_200bin_data$qvalue < 0.05 & abs(father_200bin_data$meth.diff)>= 15, ifelse(father_200bin_data$meth.diff>= 15 ,'Hyper','Hypo'),'NoSignifi'))
table(father_200bin_data$threshold)
father_200bin_data$threshold<-factor(father_200bin_data$threshold,level=c("Hyper","Hypo","NoSignifi"))
p_df15<-ggplot(data = father_200bin_data, aes(x = meth.diff, y = -log10(qvalue), colour=threshold)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("red","blue","grey")) +
  scale_x_continuous(breaks=seq(-80, 80, 20))+
  ylim(0, 600)+
  geom_vline(xintercept=c(-15,15),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="difference in DNA methlation level ",y="-log10 (q-value)",title=paste0("methykit difference: ",compare_name,"_Diff15q005_DMR")) +
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),panel.grid.major =element_blank())+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5), legend.position="right", legend.title = element_blank(),axis.line = element_line(colour="black"))+
  theme(axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
#print(p_df15);print(p_df20);print(p_df25)
ggsave(p_df15,file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_vocano_200bp_DMR_q005_df15.pdf"), width = 200, height = 200, units = "mm")
ggsave(p_df15,file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_vocano_200bp_DMR_q005_df15.png"), width = 8, height = 8)
ggsave(p_df20,file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_vocano_200bp_DMR_q005_df20.pdf"), width = 200, height = 200, units = "mm")
ggsave(p_df20,file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_vocano_200bp_DMR_q005_df20.png"), width = 8, height = 8)
ggsave(p_df25,file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_vocano_200bp_DMR_q005_df25.pdf"), width = 200, height = 200, units = "mm")
ggsave(p_df25,file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_vocano_200bp_DMR_q005_df25.png"), width = 8, height = 8)

#对DMRs的分布可视化
##difference=20,qvalue=0.05
#visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
myDiff15q005_all <- father_200bin_data[which(father_200bin_data$qvalue < 0.05 & abs(father_200bin_data$meth.diff)>= 15), ]

pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_Diff_200bp_DMR_percent_pie_plot.pdf"))
myDiff15q005_all$group<-as.factor(ifelse(myDiff15q005_all$meth.diff>0, "hyper", "hypo"))
blank_theme <- theme_minimal()+ theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
                                      panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),plot.title=element_text(size=14, face="bold"))
percete_data<-data.frame(table(myDiff15q005_all$group))
percete_data$percent<-paste0(round(percete_data$Freq / sum(percete_data$Freq) * 100, 1),"%")

ggplot(percete_data, aes(x = 1, weight = Freq, fill =   Var1)) +
  geom_bar(width = 1, position = "stack")+
  coord_polar("y", start=0)+scale_fill_manual(values=ppCor[3:4])+
  blank_theme + labs(x = "", y = "", title = paste0(compare_name,"_Diff15q005_DMR_percentage"))+
  geom_text(stat="count",aes(label = percent), size=4, position=position_stack(vjust = 0.5))
dev.off()

#rownames(myDiff15q005_all)<-paste(myDiff15q005_all$chr,myDiff15q005_all$start,myDiff15q005_all$end,sep = ".")
myDiff15q005_all_bed<- as.data.frame(str_split_fixed(rownames(myDiff15q005_all), "[.]", 3))
colnames(myDiff15q005_all_bed)<-c("chr","start","end")
myDiff15q005_all<-cbind(myDiff15q005_all_bed,myDiff15q005_all)
myDiff15q005.hyper<-myDiff15q005_all[which(myDiff15q005_all$meth.diff>0),]
myDiff15q005.hypo<-myDiff15q005_all[which(myDiff15q005_all$meth.diff<0),]
write.csv(as.data.frame(myDiff15q005_all),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_all.csv"),quote=F, row.names=T) 
write.csv(as.data.frame(myDiff15q005.hyper),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_hyper.csv"),quote=F, row.names=T) 
write.csv(as.data.frame(myDiff15q005.hypo),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_hypo.csv"),quote=F, row.names=T) 

myDiff15q005.hyper_bed<- as.data.frame(str_split_fixed(rownames(myDiff15q005.hyper), "[.]", 3))
colnames(myDiff15q005.hyper_bed)<-c("chr","start","end")
myDiff15q005.hypo_bed<- as.data.frame(str_split_fixed(rownames(myDiff15q005.hypo), "[.]", 3))
colnames(myDiff15q005.hypo_bed)<-c("chr","start","end")
head(myDiff15q005.hypo_bed)

write.table(myDiff15q005_all_bed,paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_all.bed"),col.names=F,row.names=F,sep="\t")
write.table(myDiff15q005.hyper_bed,paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_hyper.bed"),col.names=F,row.names=F,sep="\t")
write.table(myDiff15q005.hypo_bed,paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_hypo.bed"),col.names=F,row.names=F,sep="\t")

##plot heatmaps for DMRs
compare_name<-"father_AMA_vs_Young"
#heatmap:https://www.jianshu.com/p/4918df41d561
myDiff.all<- as.data.frame(read.csv(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_all.csv"),header = T,sep = ",",row.names = 1,check.names=F))
head(myDiff.all,n=3)
#rownames(myDiff.all)<-as.character(myDiff.all$X)
#myDiff.all<-myDiff.all[,-1]
myDiff.all<-myDiff.all[order(myDiff.all$meth.diff,decreasing = T),]
head(myDiff.all);colnames(myDiff.all)

sample_list<-c("10-2-1D","16-2-4C","17-2-2E","18-2","19-2-2H","20F-4G","6-2","7F","8-2","E16F",
               "11-2","12-2-3H","13-2-1H","14-2","15-2","1F","2-2","3F","4-2","5-2")
heat_plot<-myDiff.all[,sample_list]
huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
rownames(huahua_meta)<-huahua_meta$library_code
huahua_meta2<-huahua_meta[sample_list,]
head(huahua_meta2)
colnames(heat_plot)<-huahua_meta2$analysis_name
as.character(huahua_meta2$analysis_name[order(huahua_meta2$analysis_name,decreasing = F)])
sample_list2 <-c(paste0("AMA_F_",1:10),paste0("YOUNG_F_",1:10))
heat_plot<-heat_plot[,sample_list2]
head(heat_plot)

#设定 annotations
# 生成行 列的注释
annotation_col<-data.frame(Age_group = factor(c(rep("AMA", 10),rep("Young",10))))
rownames(annotation_col) = colnames(heat_plot)

annotation_row = data.frame(Class = factor(rep(c("UP_DMRs","Down_DMRs"), c(length(which(myDiff.all$meth.diff>0)), length(which(myDiff.all$meth.diff<0))))))
rownames(annotation_row) = rownames(myDiff.all)

## 自定义分组的颜色
ann_colors = list(Age_group=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(UP_DMRs=ppCor[3],Down_DMRs=ppCor[4]))
p1<-pheatmap(heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = annotation_col, annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             gaps_row = length(which(myDiff.all$meth.diff>0)),
             gaps_col =10,
             cutree_col = 2,treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","orange","firebrick3"))(50),
             main ="DNA methylation level of DMRs(Def:15% & q<0.05) original value",angle_col ="90")
#scale by row
p2<-pheatmap(heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = annotation_col,annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             gaps_row = length(which(myDiff.all$meth.diff>0)),
             gaps_col =10,cutree_col = 2,
             treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             scale ="row", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(50),
             main ="DNA methylation level of DMRs(Def:15% & q<0.05):scale by row",angle_col ="90")

pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_6X_CpGs_200bp_DMRs_myDiff15q005_heatmap1.pdf"))
print(p1)
dev.off()

pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_6X_CpGs_200bp_DMRs_myDiff15q005_heatmap2.pdf"))
print(p2)
dev.off()
