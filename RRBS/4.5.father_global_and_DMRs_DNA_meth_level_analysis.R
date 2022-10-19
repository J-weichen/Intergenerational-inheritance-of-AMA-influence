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
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
head(huahua_meta)
###plot DNA methylation in global level
#for father
father_200bin_data<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv",row.names = 1,check.names = F)
head(father_200bin_data);dim(father_200bin_data)
data_region_raw<-father_200bin_data[,1:22]
head(data_region_raw)

methy_CpG_region_plot<-as.data.frame(t(data_region_raw[,1:20]))
head(methy_CpG_region_plot[1:7,1:6])
methy_CpG_region_plot$library_code<-rownames(methy_CpG_region_plot)
methy_CpG_region_plot[1:6,1:4]
##数据整形为长数据
methy_CpG_region_plot_long <- melt(methy_CpG_region_plot,variable.name="Region",value.name = "Num",id.vars = c("library_code"))
methy_CpG_region_plot_long2<-merge(methy_CpG_region_plot_long,huahua_meta,by="library_code")
head(methy_CpG_region_plot_long2)

methy_CpG_region_plot_long2$Age_group <-factor(methy_CpG_region_plot_long2$Age_group,level=c("YOUNG","AMA"))
#methy_CpG_region_plot_long2$sample <-factor(methy_CpG_region_plot_long2$sample,level=rownames(methy_CpG_region_plot))
methy_CpG_region_plot_long2$analysis_name<-factor(methy_CpG_region_plot_long2$analysis_name,levels=c(paste0("YOUNG_F_",1:10),paste0("AMA_F_",1:10)))
head(methy_CpG_region_plot_long2)
dim(methy_CpG_region_plot_long2)#14038740       10

grid.newpage(); #清空画板，开始画新图
global_meth_plot1<-ggviolin(methy_CpG_region_plot_long2, x="analysis_name", y="Num",fill="Age_group",
                            #color="group", shape="group", 
                            palette = ppCor[2:1],
                            add = "boxplot",add.params = list( width=0.05,fill="white"))+
  xlab("sample_group") +ylab("DNA_methylation")+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(angle=90,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))
#theme(plot.title = element_text(size=10,colour = "blue",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
global_meth_plot2<-ggplot(methy_CpG_region_plot_long2, aes(x =analysis_name,y = Num,fill = Age_group))+
  theme(axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5))+
  # stat_boxplot(geom = "errorbar",width=0.2,aes(color=Group))+ geom_boxplot(alpha = .5,width=0.6)
  #  geom_violin(alpha = .5)+stat_boxplot(geom = "errorbar",width=0.2,aes(color=subtype))+ geom_boxplot(alpha = .9,width=0.2)
  xlab("sample_group") +ylab("DNA_methylation")+scale_fill_manual(values=ppCor[2:1])+
  geom_violin(alpha = .5)+stat_boxplot(geom = "errorbar",width=0.2,color = c("black"))+ 
  geom_boxplot(alpha = .5,fill="white",width=.2)
global_meth_plot1;global_meth_plot2 
ggsave("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample.pdf",global_meth_plot1,width=10, height=6)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample.png",global_meth_plot1,width=10, height=6)

ggsave("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample2.pdf",global_meth_plot2,width=10, height=6)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample2.png",global_meth_plot2,width=10, height=6)

#计算单个样本平均甲基化水平
##使用平均值进行绘图
data_region1<-as.data.frame(colMeans(data_region_raw,na.rm = T))
data_region1$library_code<-rownames(data_region1)
colnames(data_region1)<-c("mean_all_bin","library_code")

##数据整形为长数据
methy_CpG_region_plot_long <- data_region1[1:20,]
methy_CpG_region_plot_long2<-merge(methy_CpG_region_plot_long,huahua_meta,by="library_code")
head(methy_CpG_region_plot_long2)
methy_CpG_region_plot_long2$Age_group <-factor(methy_CpG_region_plot_long2$Age_group,level=c("YOUNG","AMA"))
head(methy_CpG_region_plot_long2);dim(methy_CpG_region_plot_long2)
compare_means(mean_all_bin~Age_group, data=methy_CpG_region_plot_long2,method = "wilcox.test",group_by="group")
#.y.          group1 group2     p p.adj p.format p.signif method  
#1 mean_all_bin YOUNG  AMA    0.315  0.31 0.31     ns       Wilcoxon

global_meth_plot6 <- ggboxplot(methy_CpG_region_plot_long2, x = "Age_group", y = "mean_all_bin",
                               color = "Age_group", palette = ppCor[2:1],
                               add = "jitter")+xlab("sample group") +ylab("DNA methylation")+ ylim(c(70,80))
global_meth_plot6
#  Add p-value
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_mean_all_site.pdf")
global_meth_plot6+stat_compare_means(method = "t.test",comparisons= list(c("YOUNG", "AMA")))
global_meth_plot6+stat_compare_means(method = "wilcox.test",comparisons= list(c("YOUNG", "AMA")))
#p+stat_compare_means(aes(label=..p.signif..), label.x = 1.5, label.y = 1.1)
global_meth_plot6+stat_compare_means(method = "t.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 78)
global_meth_plot6+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 78)
dev.off()

##使用单一位点的平均值进行绘图
data_region1<-as.data.frame(t(data_region_raw[,21:22]))
data_region1$group <-c("Young", "AMA")
data_region1[,1:2]

##数据整形为长数据
methy_CpG_region_plot_long <- melt(data_region1,variable.name="Region",value.name = "Num",id.vars = c("group"))
head(methy_CpG_region_plot_long)
dim(methy_CpG_region_plot_long)
methy_CpG_region_plot_long$group <-factor(methy_CpG_region_plot_long$group,level=c("Young","AMA"))

global_meth_plot7<-ggplot(methy_CpG_region_plot_long, aes(x = group,y = Num,fill = group))+
  # theme(axis.text.x=element_text(hjust=0.5, vjust=0.5))+
  # stat_boxplot(geom = "errorbar",width=0.2,aes(color=Group))+ geom_boxplot(alpha = .5,width=0.6)
  #  geom_violin(alpha = .5)+stat_boxplot(geom = "errorbar",width=0.2,aes(color=subtype))+ geom_boxplot(alpha = .9,width=0.2)
  geom_violin(alpha = .5)+stat_boxplot(geom = "errorbar",width=0.2,color = c("black"))+ geom_boxplot(alpha = .5,fill="white",width=.2)+
  #  stat_compare_means(comparisons = c("A","P")) #不同组间的比较
  xlab("sample_group") +ylab("DNA_methylation")+scale_fill_manual(values=ppCor[2:1])+
  scale_x_discrete(labels=function(x) str_wrap(x, width=20))+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))
#theme(plot.title = element_text(size=10,colour = "blue",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
global_meth_plot7

compare_means(Num~group, data=methy_CpG_region_plot_long,method = "wilcox.test",group_by="group")
compare_means(Num~group, data=methy_CpG_region_plot_long,method = "t.test",group_by="group")
global_meth_plot8<-ggviolin(methy_CpG_region_plot_long, x="group", y="Num",fill="group",
                            #color="group",shape="group", 
                            palette = ppCor[2:1],
                            add = c("boxplot"), add.params = list(width=0.05,fill="white"))+
  xlab("sample_group") +ylab("DNA_methylation")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=20))+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))
#theme(plot.title = element_text(size=10,colour = "blue",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_mean_value_distribution.pdf")
global_meth_plot7;global_meth_plot8
global_meth_plot7+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 100)
global_meth_plot8+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 100)
dev.off()

##使用标准差进行绘图
head(data_region_raw)
sd_Young<-apply(data_region_raw[,1:10], 1, sd,na.rm=T)
sd_AMA<-apply(data_region_raw[,11:20], 1, sd,na.rm=T)
data_region1<-data.frame(sd_Young,sd_AMA)
head(data_region1)
#rownames(data_region1)<-paste(data_region_raw$chrom,"_",data_region_raw$start,"_",data_region_raw$end,sep="")
data_region2<-as.data.frame(t(data_region1))
data_region2[,1:2]
#data_region3<-data_region2*100
data_region2$group <-c("Young","AMA")
dim(data_region2)
##数据整形为长数据
methy_CpG_region_plot_long <- melt(data_region2,variable.name="Region",value.name = "Num",id.vars = c("group"))
head(methy_CpG_region_plot_long)
dim(methy_CpG_region_plot_long)
methy_CpG_region_plot_long$group <-factor(methy_CpG_region_plot_long$group,level=c("Young","AMA"))

#绘制目标基因的目标区域甲基化水平柱状图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_sd1.pdf")
ggviolin(methy_CpG_region_plot_long, x="group", y="Num",fill="group",
         #color="Group",shape="Group", 
         palette = ppCor[2:1],
         add = "boxplot", add.params = list(fill="white"))+
  xlab("sample_group") +ylab("Standard Deviation")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=20))+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))
#theme(plot.title = element_text(size=10,colour = "blue",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
dev.off()

global_meth_plot9<-ggplot(methy_CpG_region_plot_long, aes(x = group,y = Num,fill = group))+
  # theme(axis.text.x=element_text(hjust=0.5, vjust=0.5))+
  # stat_boxplot(geom = "errorbar",width=0.2,aes(color=Group))+ geom_boxplot(alpha = .5,width=0.6)
  #  geom_violin(alpha = .5)+stat_boxplot(geom = "errorbar",width=0.2,aes(color=subtype))+ geom_boxplot(alpha = .9,width=0.2)
  geom_violin(alpha = .5)+stat_boxplot(geom = "errorbar",width=0.2,color = c("black"))+ geom_boxplot(alpha = .5,fill="white",width=.2)+
  xlab("sample_group") +ylab("Standard Deviation")+scale_fill_manual(values=ppCor[2:1])+
  scale_x_discrete(labels=function(x) str_wrap(x, width=20))+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))
#theme(plot.title = element_text(size=10,colour = "blue",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
global_meth_plot9

compare_means(Num~group, data=methy_CpG_region_plot_long,method = "wilcox.test",group_by="Group")
compare_means(Num~group, data=methy_CpG_region_plot_long,method = "t.test",group_by="Group")
global_meth_plot10<-ggviolin(methy_CpG_region_plot_long, x="group", y="Num",fill="group",
                             #color="Group",shape="Group", 
                             palette = ppCor[2:1],
                             add = "boxplot", add.params = list(fill="white"))+
  xlab("sample_group") +ylab("Standard Deviation(%)")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=20))+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_sd2.pdf")
global_meth_plot10+stat_compare_means(method = "t.test",aes(label=paste0(..method..,": "," p  ",..p.format..,"\n",..p.signif..)), label.x = 1.5, label.y = 40)
global_meth_plot10+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p  ",..p.format..,"\n",..p.signif..)), label.x = 1.5, label.y = 40)
global_meth_plot9+stat_compare_means(method = "t.test",aes(label=paste0(..method..,": "," p  ",..p.format..,"\n",..p.signif..)), label.x = 1.5, label.y = 40)
global_meth_plot9+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p  ",..p.format..,"\n",..p.signif..)), label.x = 1.5, label.y = 40)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_sd3.pdf")
global_meth_plot9+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p  ",..p.format..,"\n",..p.signif..)), label.x = 1.5, label.y = 40)
dev.off()
global_meth_plot_sd<-global_meth_plot9+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p  ",..p.format..,"\n",..p.signif..)), label.x = 1.5, label.y = 40)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_sd3.png",global_meth_plot_sd,width = 6, height =6)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Boxplot_father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_sd3_2.png",global_meth_plot9,width = 6, height =6)

##calculation of DNA methylation level of DMRs
#for father 
father_200bin_data<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv",row.names = 1,check.names = F)
father_200bin_data$bin_region <- unlist(lapply(strsplit(as.character(rownames(father_200bin_data)), "[.]"), function(x) paste(x[1],x[2],x[3],sep="_")))
rownames(father_200bin_data)<-father_200bin_data$bin_region
head(father_200bin_data)
dim(father_200bin_data)

#head(data_raw_bin2)
compare_name<-"father_AMA_vs_Young"
DMRs_hyper<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_hyper.bed"),sep="\t")
DMRs_hypo <-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_200bp_DMRs_myDiff15q005_hypo.bed"),sep="\t")
head(DMRs_hyper)
DMRs_hyper_region<-paste(DMRs_hyper$V1,DMRs_hyper$V2,DMRs_hyper$V3,sep="_")
DMRs_hypo_region<-paste(DMRs_hypo$V1,DMRs_hypo$V2,DMRs_hypo$V3,sep="_")
my_comparisons <- list(c("YOUNG","AMA"))
length(DMRs_hyper_region);length(DMRs_hypo_region)#217 #268

##plot DNA methylation in DMRs
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/father_meth_levels_in_6X_DMRs_myDiff15q005_multiple_plot.pdf")

#个体分开绘制与保留所有点
data_raw_bin<-father_200bin_data[,1:20]
DMRs_hyper_data<-data_raw_bin[DMRs_hyper_region,]
DMRs_hypo_data<-data_raw_bin[DMRs_hypo_region,]

meth_dataframe_region<-DMRs_hyper_data
methy_CpG_region_plot_long <- melt(meth_dataframe_region,variable.name="sample",value.name = "Num")
colnames(methy_CpG_region_plot_long) <- c("library_code","mean")
head(methy_CpG_region_plot_long)
methy_CpG_region_plot_long2<-merge(methy_CpG_region_plot_long,huahua_meta,by="library_code")
head(methy_CpG_region_plot_long2)

methy_CpG_region_plot_long2$Age_group <-factor(methy_CpG_region_plot_long2$Age_group,level=c("YOUNG","AMA"))
head(methy_CpG_region_plot_long2);dim(methy_CpG_region_plot_long2)
compare_means(mean~Age_group, data=methy_CpG_region_plot_long2, method = "wilcox.test",group_by="Age_group")
#   .y.  group1 group2        p   p.adj p.format p.signif method  
#1  mean  YOUNG  AMA    7.57e-81 7.6e-81 <2e-16   ****     Wilcoxon

P_mean_site_hyper <- ggboxplot(methy_CpG_region_plot_long2, x = "analysis_name", y = "mean",color = "Age_group", palette = ppCor[2:1],add = "jitter")+
  labs(x = "sample_group", y = "DNA_methylation", title = "DNA_methylation_in hyper_df15q005_DMRs")+
  #xlab("sample_group") +ylab("DNA_methylation")+title("DNA_methylation_in hyper_df50q005_DMRs")+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
P_mean_site_hyper

P_mean_site_hyper1 <- ggboxplot(methy_CpG_region_plot_long2, x = "Age_group", y = "mean",color = "Age_group", palette = ppCor[2:1],add = "jitter")+
  labs(x = "sample_group", y = "DNA_methylation", title = "DNA_methylation_in hyper_df15q005_DMRs")+
  #xlab("sample_group") +ylab("DNA_methylation")+title("DNA_methylation_in hyper_df50q005_DMRs")+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
P_mean_site_hyper1

P_mean_site_hyper2<-ggviolin(methy_CpG_region_plot_long2, x="analysis_name", y="mean",color="Age_group",
                             palette = ppCor[2:1],add = "boxplot", add.params = list(fill="white"))+
  labs(x = "sample_group", y = "DNA_methylation", title = "DNA_methylation_in hypo_df15q005_DMRs")+
  #xlab("sample_group") +ylab("DNA_methylation")+title("DNA_methylation_in_hypo_df50q005_DMRs")+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
P_mean_site_hyper2

meth_dataframe_region<-DMRs_hypo_data
methy_CpG_region_plot_long <- melt(meth_dataframe_region,variable.name="sample",value.name = "Num")
colnames(methy_CpG_region_plot_long) <- c("library_code","mean")
head(methy_CpG_region_plot_long)
methy_CpG_region_plot_long2<-merge(methy_CpG_region_plot_long,huahua_meta,by="library_code")
head(methy_CpG_region_plot_long2)

P_mean_site_hypo <- ggboxplot(methy_CpG_region_plot_long2, x = "analysis_name", y = "mean",color = "Age_group", palette = ppCor[1:2],add = "jitter")+
  labs(x = "sample_group", y = "DNA_methylation", title = "DNA_methylation_in hypo_df15q005_DMRs")+
  #xlab("sample_group") +ylab("DNA_methylation")+title("DNA_methylation_in_hypo_df50q005_DMRs")+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
P_mean_site_hypo

P_mean_site_hypo1 <- ggboxplot(methy_CpG_region_plot_long2, x = "analysis_name", y = "mean",color = "Age_group",palette = ppCor[1:2],add = "jitter")+
  labs(x = "sample_group", y = "DNA_methylation", title = "DNA_methylation_in hypo_df15q005_DMRs")+
  # xlab("sample_group") +ylab("DNA_methylation")+title("DNA_methylation_in_hypo_df50q005_DMRs")+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
P_mean_site_hypo1

P_mean_site_hypo2<-ggviolin(methy_CpG_region_plot_long2, x = "analysis_name", y = "mean",color = "Age_group",
                            palette = ppCor[1:2],add = "boxplot", add.params = list(fill="white"))+
  labs(x = "sample_group", y = "DNA_methylation", title = "DNA_methylation_in hypo_df15q005_DMRs")+
  # xlab("sample_group") +ylab("DNA_methylation")+title("DNA_methylation_in_hypo_df50q005_DMRs")+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
#dev.off()
P_mean_site_hypo2

#method 3
#计算单一样本平均甲基化水平
my_comparisons <- list(c("YOUNG","AMA"))
data_raw_bin3<-father_200bin_data[,1:20]

#多样本平均甲基化水平
DMRs_hyper_data<-data_raw_bin3[DMRs_hyper_region,];DMRs_hypo_data<-data_raw_bin3[DMRs_hypo_region,]
DMRs_hyper_data2<-data.frame(colMeans(DMRs_hyper_data[,1:20],na.rm = T));DMRs_hypo_data2<-data.frame(colMeans(DMRs_hypo_data[,1:20],na.rm = T))
DMRs_hyper_data2$library_code<-rownames(DMRs_hyper_data2);DMRs_hypo_data2$library_code<-rownames(DMRs_hypo_data2)
colnames(DMRs_hyper_data2) <- colnames(DMRs_hypo_data2) <-c("mean","library_code")
DMRs_hyper_data3<-merge(DMRs_hyper_data2,huahua_meta,by="library_code")
DMRs_hypo_data3<-merge(DMRs_hypo_data2,huahua_meta,by="library_code")
head(DMRs_hypo_data3);head(DMRs_hyper_data3)

P_mean_site_hyper6 <- ggboxplot(DMRs_hyper_data3, x = "Age_group", y = "mean",color = "Age_group", palette = ppCor[1:6],add = "jitter")+
  labs(x = "sample_group", y = "DNA_methylation", title = "DNA_methylation_in hyper_df50q005_DMRs")+
  #xlab("sample_group") +ylab("DNA_methylation")+title("DNA_methylation_in_hyper_df50q005_DMRs_mean")+ ylim(0,100)+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "blue",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
P_mean_site_hyper6

P_mean_site_hypo6 <- ggboxplot(DMRs_hypo_data3, x = "Age_group", y = "mean",color = "Age_group", palette = ppCor[1:6],add = "jitter")+
  labs(x = "sample_group", y = "DNA_methylation", title = "DNA_methylation_in hypo_df15q005_DMRs")+
  #xlab("sample_group") +ylab("DNA_methylation")+title("DNA_methylation_in_hypo_df50q005_DMRs_mean")+ ylim(0,100)+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "blue",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
P_mean_site_hypo6

P_mean_site_hyper6 +stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..)), label.x = 1.4, label.y =90)#+labs(title = "wilcox.test")
P_mean_site_hypo6 +stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..)), label.x = 1.4, label.y =90)#+labs(title = "wilcox.test")
dev.off()
