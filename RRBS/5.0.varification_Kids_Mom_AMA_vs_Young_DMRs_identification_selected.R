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

#############for kids
#selection for DMRs FOR AMA_abortion vs AMA_arrest
kids_200bin_data<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth_PVALUE.csv",row.names = 1,check.names=F)
head(kids_200bin_data);dim(kids_200bin_data)
#602367     16

##test evaluation
str(kids_200bin_data)
myDiffq005_all<-kids_200bin_data[kids_200bin_data$qvalue <0.05,]
dim(myDiffq005_all)#113304     16
myDiffq005_all$group1<-as.factor(ifelse(myDiffq005_all$meth.diff>0, "hyper", "hypo"))
myDiffq005_all$group_mean<-as.factor(ifelse(myDiffq005_all$meth_mean_def>0, "hyper", "hypo"))
head(myDiffq005_all);dim(myDiffq005_all)

pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_kids_AMA_vs_Young_200bin_q005_dif15_distribution.pdf")
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
#    -45.7822  -4.7191  -2.9662  -2.1935  -0.9941  40.3076 

summary(myDiffq005_all[myDiffq005_all$group1=="hyper",]$meth.diff)
#      Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.3191  3.3597  5.1538  6.0517  7.7336 40.3076 
summary(myDiffq005_all[myDiffq005_all$group1=="hypo",]$meth.diff)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-45.7822  -5.2469  -3.5557  -4.2032  -2.3585  -0.2022 

library(psych)
describe(myDiffq005_all$meth.diff)
##非缺失值的数量、平均数、标准差、中位数、截尾均值、绝对中位差、最小值、最大值、值域、偏度、峰度和平均值的标准误
##vars     n    mean    sd median trimmed   mad    min   max  range skew kurtosis   se
## X1     1 113304 -2.19 5.15  -2.97   -2.57 2.73 -45.78 40.31 86.09 0.75     3.73 0.02

describe(myDiffq005_all[myDiffq005_all$group1=="hyper",]$meth.diff)
#vars     n     mean    sd median trimmed  mad  min   max range skew kurtosis   se
##X1    1 22205 6.05 3.86   5.15    5.53 3.05 0.32 40.31 39.99 1.88     6.24 0.03
describe(myDiffq005_all[myDiffq005_all$group1=="hypo",]$meth.diff)
#vars        n   mean   sd median trimmed  mad    min   max range  skew kurtosis   se
##X1       1 91099 -4.2 2.96  -3.56    -3.8 2.04 -45.78 -0.2 45.58 -2.47    12.33 0.01

pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_kids_AMA_vs_Young_Diff_200bp_DMR_pie50_25_20_10_5_def_plot.pdf")
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=5,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=10,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=15,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=20,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=25,]$group1))
dev.off()

#pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/AMA_arrest_vs_AMA_abortion_q005_Diff_200bp_DMR_vocano_plot.pdf")
#https://www.jianshu.com/p/41375176906e

kids_200bin_data$threshold = as.factor(ifelse(kids_200bin_data$qvalue < 0.05 & abs(kids_200bin_data$meth.diff)>= 50, ifelse(kids_200bin_data$meth.diff>= 50 ,'Hyper','Hypo'),'NoSignifi'))
table(kids_200bin_data$threshold)
#NoSignifi 
#602367 
kids_200bin_data$threshold = as.factor(ifelse(kids_200bin_data$qvalue < 0.05 & abs(kids_200bin_data$meth.diff)>= 25, ifelse(kids_200bin_data$meth.diff>= 25 ,'Hyper','Hypo'),'NoSignifi'))
table(kids_200bin_data$threshold)
#Hyper      Hypo NoSignifi 
#74        91    602202 
kids_200bin_data$threshold = as.factor(ifelse(kids_200bin_data$qvalue < 0.05 & abs(kids_200bin_data$meth.diff)>= 20, ifelse(kids_200bin_data$meth.diff>= 20 ,'Hyper','Hypo'),'NoSignifi'))
table(kids_200bin_data$threshold)
# Hyper      Hypo NoSignifi 
#  208       295    601864 
kids_200bin_data$threshold = as.factor(ifelse(kids_200bin_data$qvalue < 0.05 & abs(kids_200bin_data$meth.diff)>= 15, ifelse(kids_200bin_data$meth.diff>= 15 ,'Hyper','Hypo'),'NoSignifi'))
table(kids_200bin_data$threshold)
#  Hyper      Hypo NoSignifi 
#  690       940    600737
kids_200bin_data$threshold = as.factor(ifelse(kids_200bin_data$qvalue < 0.05 & abs(kids_200bin_data$meth.diff)>= 10, ifelse(kids_200bin_data$meth.diff>= 10 ,'Hyper','Hypo'),'NoSignifi'))
table(kids_200bin_data$threshold)
#Hyper      Hypo NoSignifi 
# 2877      3995    595495 

#Final usage:
compare_name<-"varification_kids_AMA_vs_Young"
kids_200bin_data$threshold = as.factor(ifelse(kids_200bin_data$qvalue < 0.05 & abs(kids_200bin_data$meth.diff)>= 25, ifelse(kids_200bin_data$meth.diff>= 25 ,'Hyper','Hypo'),'NoSignifi'))
table(kids_200bin_data$threshold)
kids_200bin_data$threshold<-factor(kids_200bin_data$threshold,level=c("Hyper","Hypo","NoSignifi"))
p_df25<-ggplot(data = kids_200bin_data, aes(x = meth.diff, y = -log10(qvalue), colour=threshold)) +
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

kids_200bin_data$threshold = as.factor(ifelse(kids_200bin_data$qvalue < 0.05 & abs(kids_200bin_data$meth.diff)>= 20, ifelse(kids_200bin_data$meth.diff>= 20 ,'Hyper','Hypo'),'NoSignifi'))
table(kids_200bin_data$threshold)
kids_200bin_data$threshold<-factor(kids_200bin_data$threshold,level=c("Hyper","Hypo","NoSignifi"))
p_df20<-ggplot(data = kids_200bin_data, aes(x = meth.diff, y = -log10(qvalue), colour=threshold)) +
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

kids_200bin_data$threshold = as.factor(ifelse(kids_200bin_data$qvalue < 0.05 & abs(kids_200bin_data$meth.diff)>= 15, ifelse(kids_200bin_data$meth.diff>= 15 ,'Hyper','Hypo'),'NoSignifi'))
table(kids_200bin_data$threshold)
kids_200bin_data$threshold<-factor(kids_200bin_data$threshold,level=c("Hyper","Hypo","NoSignifi"))
p_df15<-ggplot(data = kids_200bin_data, aes(x = meth.diff, y = -log10(qvalue), colour=threshold)) +
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
myDiff15q005_all <- kids_200bin_data[which(kids_200bin_data$qvalue < 0.05 & abs(kids_200bin_data$meth.diff)>= 15), ]

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
#heatmap:https://www.jianshu.com/p/4918df41d561
myDiff.all<- as.data.frame(read.csv(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_all.csv"),header = T,sep = ",",row.names = 1,check.names=F))
head(myDiff.all,n=3)
#rownames(myDiff.all)<-as.character(myDiff.all$X)
#myDiff.all<-myDiff.all[,-1]
myDiff.all<-myDiff.all[order(myDiff.all$meth.diff,decreasing = T),]
head(myDiff.all);colnames(myDiff.all)

sample_list<-c( "Y1Q","Y2Q","Y4Q","Y5Q","Y6Q","A2Q","A3Q","A4Q","A5Q","A6Q")
heat_plot<-myDiff.all[,sample_list]

colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/verification_AMA_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)
colnames(heat_plot)
#colnames(heat_plot)<-colData_whole$sample 

#设定 annotations
# 生成行 列的注释
annotation_col<-data.frame(Age_group = factor(c(rep("Young",5),rep("AMA", 5))))
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

#############for mother
#selection for DMRs FOR AMA_abortion vs AMA_arrest
mother_200bin_data<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth_PVALUE.csv",row.names = 1,check.names=F)
head(mother_200bin_data);dim(mother_200bin_data)
#602367     16

##test evaluation
str(mother_200bin_data)
myDiffq005_all<-mother_200bin_data[mother_200bin_data$qvalue <0.05,]
dim(myDiffq005_all)#83958    16
myDiffq005_all$group1<-as.factor(ifelse(myDiffq005_all$meth.diff>0, "hyper", "hypo"))
myDiffq005_all$group_mean<-as.factor(ifelse(myDiffq005_all$meth_mean_def>0, "hyper", "hypo"))
head(myDiffq005_all);dim(myDiffq005_all)

pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_mother_AMA_vs_Young_200bin_q005_dif15_distribution.pdf")
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
#    -38.799  -4.953  -2.638  -1.678   2.226  39.538 

summary(myDiffq005_all[myDiffq005_all$group1=="hyper",]$meth.diff)
#      Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2973  3.0149  4.6230  5.6845  7.1914 39.5376 
summary(myDiffq005_all[myDiffq005_all$group1=="hypo",]$meth.diff)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-38.7991  -5.9840  -3.8273  -4.6026  -2.3215  -0.1828 

library(psych)
describe(myDiffq005_all$meth.diff)
##非缺失值的数量、平均数、标准差、中位数、截尾均值、绝对中位差、最小值、最大值、值域、偏度、峰度和平均值的标准误
##vars     n    mean    sd median trimmed   mad    min   max  range skew kurtosis   se
## X1     1 83958 -1.68 5.86  -2.64   -1.93 4.11 -38.8 39.54 78.34  0.5     2.29 0.02

describe(myDiffq005_all[myDiffq005_all$group1=="hyper",]$meth.diff)
#vars     n     mean    sd median trimmed  mad  min   max range skew kurtosis   se
##X1    1 23870 5.68 3.96   4.62    5.07 2.82 0.3 39.54 39.24 2.12     7.29 0.03
describe(myDiffq005_all[myDiffq005_all$group1=="hypo",]$meth.diff)
#vars        n   mean   sd median trimmed  mad    min   max range  skew kurtosis   se
##X1      1 60088 -4.6 3.42  -3.83   -4.14 2.58 -38.8 -0.18 38.62 -1.99      7.2 0.01

pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_mother_AMA_vs_Young_Diff_200bp_DMR_pie50_25_20_10_5_def_plot.pdf")
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=5,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=10,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=15,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=20,]$group1))
pie(table(myDiffq005_all[abs(myDiffq005_all$meth.diff)>=25,]$group1))
dev.off()

#pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/AMA_arrest_vs_AMA_abortion_q005_Diff_200bp_DMR_vocano_plot.pdf")
#https://www.jianshu.com/p/41375176906e

mother_200bin_data$threshold = as.factor(ifelse(mother_200bin_data$qvalue < 0.05 & abs(mother_200bin_data$meth.diff)>= 50, ifelse(mother_200bin_data$meth.diff>= 50 ,'Hyper','Hypo'),'NoSignifi'))
table(mother_200bin_data$threshold)
#NoSignifi 
#602367 
mother_200bin_data$threshold = as.factor(ifelse(mother_200bin_data$qvalue < 0.05 & abs(mother_200bin_data$meth.diff)>= 25, ifelse(mother_200bin_data$meth.diff>= 25 ,'Hyper','Hypo'),'NoSignifi'))
table(mother_200bin_data$threshold)
#Hyper      Hypo NoSignifi 
# 85        92    602190 
mother_200bin_data$threshold = as.factor(ifelse(mother_200bin_data$qvalue < 0.05 & abs(mother_200bin_data$meth.diff)>= 20, ifelse(mother_200bin_data$meth.diff>= 20 ,'Hyper','Hypo'),'NoSignifi'))
table(mother_200bin_data$threshold)
# Hyper      Hypo NoSignifi 
# 252       278    601837 
mother_200bin_data$threshold = as.factor(ifelse(mother_200bin_data$qvalue < 0.05 & abs(mother_200bin_data$meth.diff)>= 15, ifelse(mother_200bin_data$meth.diff>= 15 ,'Hyper','Hypo'),'NoSignifi'))
table(mother_200bin_data$threshold)
#  Hyper      Hypo NoSignifi 
# 795       957    600615 
mother_200bin_data$threshold = as.factor(ifelse(mother_200bin_data$qvalue < 0.05 & abs(mother_200bin_data$meth.diff)>= 10, ifelse(mother_200bin_data$meth.diff>= 10 ,'Hyper','Hypo'),'NoSignifi'))
table(mother_200bin_data$threshold)
#Hyper      Hypo NoSignifi 
# 2780      4138    595449 

#Final usage:
compare_name<-"varification_mother_AMA_vs_Young"
mother_200bin_data$threshold = as.factor(ifelse(mother_200bin_data$qvalue < 0.05 & abs(mother_200bin_data$meth.diff)>= 25, ifelse(mother_200bin_data$meth.diff>= 25 ,'Hyper','Hypo'),'NoSignifi'))
table(mother_200bin_data$threshold)
mother_200bin_data$threshold<-factor(mother_200bin_data$threshold,level=c("Hyper","Hypo","NoSignifi"))
p_df25<-ggplot(data = mother_200bin_data, aes(x = meth.diff, y = -log10(qvalue), colour=threshold)) +
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

mother_200bin_data$threshold = as.factor(ifelse(mother_200bin_data$qvalue < 0.05 & abs(mother_200bin_data$meth.diff)>= 20, ifelse(mother_200bin_data$meth.diff>= 20 ,'Hyper','Hypo'),'NoSignifi'))
table(mother_200bin_data$threshold)
mother_200bin_data$threshold<-factor(mother_200bin_data$threshold,level=c("Hyper","Hypo","NoSignifi"))
p_df20<-ggplot(data = mother_200bin_data, aes(x = meth.diff, y = -log10(qvalue), colour=threshold)) +
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

mother_200bin_data$threshold = as.factor(ifelse(mother_200bin_data$qvalue < 0.05 & abs(mother_200bin_data$meth.diff)>= 15, ifelse(mother_200bin_data$meth.diff>= 15 ,'Hyper','Hypo'),'NoSignifi'))
table(mother_200bin_data$threshold)
mother_200bin_data$threshold<-factor(mother_200bin_data$threshold,level=c("Hyper","Hypo","NoSignifi"))
p_df15<-ggplot(data = mother_200bin_data, aes(x = meth.diff, y = -log10(qvalue), colour=threshold)) +
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
myDiff15q005_all <- mother_200bin_data[which(mother_200bin_data$qvalue < 0.05 & abs(mother_200bin_data$meth.diff)>= 15), ]

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
#heatmap:https://www.jianshu.com/p/4918df41d561
myDiff.all<- as.data.frame(read.csv(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_all.csv"),header = T,sep = ",",row.names = 1,check.names=F))
head(myDiff.all,n=3)
#rownames(myDiff.all)<-as.character(myDiff.all$X)
#myDiff.all<-myDiff.all[,-1]
myDiff.all<-myDiff.all[order(myDiff.all$meth.diff,decreasing = T),]
head(myDiff.all);colnames(myDiff.all)

sample_list<-c("Y1M","Y2M","Y4M","Y5M","Y6M","A2M","A3M","A4M","A5M","A6M")
heat_plot<-myDiff.all[,sample_list]

colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/verification_AMA_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)
colnames(heat_plot)
#colnames(heat_plot)<-colData_whole$sample 

#设定 annotations
# 生成行 列的注释
annotation_col<-data.frame(Age_group = factor(c(rep("Young",5),rep("AMA", 5))))
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
