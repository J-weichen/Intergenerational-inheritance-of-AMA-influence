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

grid.newpage(); #清空画板，开始画新图
#setwd("/mnt/data/chenwei/huahua/")
myobj<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_myobj_6X_CpG.rds")
#length(myobj@treatment)
### filter sites with low (<1) or extremly high (>99.9% percentile) sequencing depth to reduce the effect of PCRs
depth<-6
##Mobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
##tiles_200=tileMethylCounts(Mobj,win.size=200,step.size=200,cov.bases =3)

#基因组划分区域进行计算 8 sample per groups
tiles_200 <-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_6X_tiles_200.rds")
meth_200bp=methylKit::unite(tiles_200, destrand=FALSE,min.per.group=8L,mc.cores = 48)
write.table(as.data.frame(meth_200bp), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_200bp_CpG_6X_3site_bismark_8sample_covarage.txt",quote=F, row.names=F, col.names=T) 
#meth_200bp <- read.table(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_200bp_CpG_6X_3site_bismark_8sample_covarage.txt",header=T) 

#将样本聚合在一起
meth<-meth_200bp
head(meth)
meth@sample.ids
unique(meth$chr)
meth2<-meth[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",meth$chr,invert=TRUE),]
unique(meth2$chr)
meth2$chr<-factor(meth2$chr,levels =as.character(unique(meth2$chr)))
dim(meth2); dim(na.omit(meth2))
#[1] 711720    184
#[1] 477743    184

saveRDS(meth2, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth2.rds")

#meth2 <-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth2.rds")
#画出样本关联度的图  #备注：： 绘图速度慢
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_similiraty.pdf")
getCorrelation(meth2,method ="spearman",plot=TRUE)
dev.off()

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth2,rowids=T))
head(data_raw_bin)
#plot count 200bp bin in each samples
huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
head(huahua_meta)
#tiles_200_2[1:6,];meth_200bp_2[1:6,]
rownames(huahua_meta)<-huahua_meta$library_code
#
#filtering and arrange data
evaluation.table<-data_raw_bin
evaluation.table[!(is.na(evaluation.table))]<- 1
evaluation.table[is.na(evaluation.table)]<- 0
bin200_detected_number<-colSums(evaluation.table)
bin200_detected_number<-data.frame(bin200_detected_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin200_detected_number) <- "bin200_number"
rownames(bin200_detected_number);rownames(huahua_meta)
bin200_detected_number2<-merge(bin200_detected_number,huahua_meta,by=0)
head(bin200_detected_number2);dim(bin200_detected_number2)
bin200_detected_number2$group2<-paste(bin200_detected_number2$Age_group,bin200_detected_number2$Sample.types,sep="_")

bin200_detected_number2$Age_group <-factor(bin200_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin200_detected_number2$Sample.types <-factor(bin200_detected_number2$Sample.types,level=c("kids","mother","father"))
bin200_detected_number2$group2 <-factor(bin200_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin200_detected_number2$Gender_RRBS<-factor(bin200_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin200_detected_number2)
bin200_detected_number2$analysis_name<-as.character(bin200_detected_number2$analysis_name)
bin200_detected_number2<-bin200_detected_number2[order(bin200_detected_number2$group2,bin200_detected_number2$bin200_number,decreasing = T),]
bin200_detected_number2$analysis_name<-factor(bin200_detected_number2$analysis_name,levels=as.character(bin200_detected_number2$analysis_name))
head(bin200_detected_number2)
#plot
plot_bin200_detected_number0<-ggplot(data=bin200_detected_number2, mapping=aes(x=analysis_name,y=bin200_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin200_detected_number1 <-plot_bin200_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin200_number), color="black", size=1,position=position_stack(1.01))+
  theme_classic()+labs(x="",y="Number of 200bp bin",title="Number of bin detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin200_detected_number1
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.number_of_200bp_bins_6X_CpG_3L_detected_in_per_group_8_sample.pdf",plot_bin200_detected_number1,width = 16, height = 8)

#plot
plot_bin200_detected_number0<-ggplot(data=bin200_detected_number2, mapping=aes(x=analysis_name,y=bin200_number,fill=type_Batch))+geom_bar(stat="identity",width=0.8)
plot_bin200_detected_number2 <-plot_bin200_detected_number0+scale_fill_manual(values=ppCor[c(8,3)])+
  geom_text(stat="identity",aes(label=bin200_number), color="black", size=1,position=position_stack(1.01))+
  theme_classic()+labs(x="",y="Number of 200bp bin",title="Number of bin detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin200_detected_number2
#plot
plot_bin200_detected_number0<-ggplot(data=bin200_detected_number2, mapping=aes(x=analysis_name,y=bin200_number,fill=seq_Batch))+geom_bar(stat="identity",width=0.8)
plot_bin200_detected_number3 <-plot_bin200_detected_number0+scale_fill_manual(values=ppCor)+
  geom_text(stat="identity",aes(label=bin200_number), color="black", size=1,position=position_stack(1.01))+
  theme_classic()+labs(x="",y="Number of 200bp bin",title="Number of bin detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin200_detected_number3
count_plot<-ggarrange(plot_bin200_detected_number1,plot_bin200_detected_number2,plot_bin200_detected_number3,labels = c("A", "B","C"),ncol = 1, nrow = 3)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.number_of_200bp_bins_6X_CpG_3L_detected_in_per_group_8_sample_three.pdf",count_plot,width = 16, height =20)

#mean DNA methylation level after filters
colnames(data_raw_bin)
data_raw_bin$AMA_Kids_mean<-apply(data_raw_bin[,1:10],1,function(x) mean(x,na.rm=T))
data_raw_bin$YOUNG_Kids_mean<-apply(data_raw_bin[,11:20],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_Mother_mean<-apply(data_raw_bin[,21:30],1,function(x) mean(x,na.rm=T))
data_raw_bin$YOUNG_Mother_mean<-apply(data_raw_bin[,31:40],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_father_mean<-apply(data_raw_bin[,41:50],1,function(x) mean(x,na.rm=T))
data_raw_bin$YOUNG_father_mean<-apply(data_raw_bin[,51:60],1,function(x) mean(x,na.rm=T))
dim(data_raw_bin)
write.table(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.txt",quote=F,row.names=T)
write.csv(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.csv",quote=F, row.names=T)
dim(data_raw_bin[,61:66])#701937      6
head(data_raw_bin[,61:66])#701937      6

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_similarity_logisty.pdf")
p2<-ggplot(data = na.omit(data_raw_bin[,61:62]), mapping = aes(x = AMA_Kids_mean, y = YOUNG_Kids_mean)) + 
  stat_bin2d(bins = 300) + scale_fill_gradient(low = 'steelblue', high = 'darkred', limits = c(0,100), breaks = c(0,25,50,100))
p_plot1 <- p2+ stat_cor(method = "pearson", label.x = 10, label.y = 90)+xlim(0,100)+ xlab("AMA_Kids_mean") + 
  theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5,  hjust = 0.5))+ylab("YOUNG_Kids_mean") + 
  theme(axis.title.y = element_text(size = 16,face = "bold", vjust = 0.5, hjust = 0.5))+theme_bw()
print(p_plot1)

p4<-ggplot(data = na.omit(data_raw_bin[,63:64]), mapping = aes(x = AMA_Mother_mean, y = YOUNG_Mother_mean)) + 
  stat_bin2d(bins = 300) + scale_fill_gradient(low = 'steelblue', high = 'darkred', limits = c(0,100), breaks = c(0,25,50,100))
p_plot3 <- p4+ stat_cor(method = "pearson", label.x = 10, label.y = 90)+xlim(0,100)+ xlab("AMA_Mother_mean") + 
  theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5,  hjust = 0.5))+ylab("YOUNG_Mother_mean") + 
  theme(axis.title.y = element_text(size = 16,face = "bold", vjust = 0.5, hjust = 0.5))+theme_bw()
print(p_plot3)

p3<-ggplot(data = na.omit(data_raw_bin[,65:66]), mapping = aes(x = AMA_father_mean, y = YOUNG_father_mean)) + 
  stat_bin2d(bins = 300) + scale_fill_gradient(low = 'steelblue', high = 'darkred', limits = c(0,100), breaks = c(0,25,50,100))
p_plot4 <- p3+ stat_cor(method = "pearson", label.x = 10, label.y = 90)+xlim(0,100)+ xlab("AMA_father_mean") + 
  theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5,  hjust = 0.5))+ylab("YOUNG_father_mean") + 
  theme(axis.title.y = element_text(size = 16,face = "bold", vjust = 0.5, hjust = 0.5))+theme_bw()
print(p_plot4)

dev.off()

#画出样本聚类的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_clustering.pdf",width=12, height=5)
clusterSamples(meth2, dist="correlation", method="ward", plot=TRUE)
clusterSamples(na.omit(meth2), dist="correlation", method="ward", plot=TRUE)
clusterSamples(meth2, dist="correlation", method="complete", plot=TRUE)
clusterSamples(na.omit(meth2), dist="correlation", method="complete", plot=TRUE)

clusterSamples(meth2, dist="euclidean", method="ward", plot=TRUE)
clusterSamples(na.omit(meth2), dist="euclidean", method="ward", plot=TRUE)
clusterSamples(meth2, dist="euclidean", method="complete", plot=TRUE)
clusterSamples(na.omit(meth2), dist="euclidean", method="complete", plot=TRUE)
dev.off()

#画出样本主成分分析图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot.pdf")
PCASamples(meth2, screeplot=TRUE)
dev.off()

#画出样本主成分分析后的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot2.pdf")
PCASamples(meth2)
dev.off()

##PCA plot
data_raw<-as.data.frame(percMethylation(meth2,rowids=T))
data_raw<-data_raw[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",rownames(data_raw),invert=TRUE),]
data2<- t(data_raw)
set.seed(12)
md <- prep(data2, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 10)
summary(resPPCA) 

##PCA score:    PC1    PC2     PC3     PC4     PC5     PC6    PC7     PC8     PC9    PC10
#R2            0.1662 0.06978 0.05069 0.03624 0.03334 0.02775 0.02288 0.01943 0.01796 0.01718
#Cumulative R2 0.1662 0.23598 0.28666 0.32290 0.35624 0.38399 0.40687 0.42630 0.44426 0.46144
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCA.txt",quote=F, row.names=T, col.names=T)

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
rownames(huahua_meta)<-huahua_meta$library_code
dim(huahua_meta)

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
dim(df)
rownames(df);rownames(huahua_meta)
df1<-merge(df,huahua_meta,by=0)
dim(df1)

variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
df1$group2<-paste(df1$Age_group,df1$Sample.types,sep="_")
df1$group2 <-factor(df1$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))

head(df1)
str(df1)
df1<-df1[order(df1$group2,df1$analysis_name,decreasing = F),]
head(df1);dim(df1)
rownames(df1)<-df1$analysis_name
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot_2D_3D.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(7,10),rep(8,10),rep(2,10),rep(1,10),rep(3,10),rep(5,10))],
                                  pch = c(rep(16,60)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
legend("right",pch=16, legend = levels(df1$group2),col =ppCor[c(7:8,2:1,3,5)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=analysis_name), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(7:8,2:1,3,5)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=analysis_name), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(7:8,2:1,3,5)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()

####evaluation by kids mother and father
meth2<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth2.rds")
meth2@sample.ids
#########################
#for kids
sample_list<-c("10-3-1E","16-3-4D","17-3-2F","18-3","19-3-4E","20Q-4H","6-3","7-3","8-3","E16C",
               "11-3-3F","12-3-4A","13-3-1G","14-3","15-3-2C","1Q","2-3-3C","3-3","4-3","5-3")
meth_kids =reorganize(meth2,sample.ids=sample_list,treatment=c(rep(0,10),rep(1,10)))

#画出样本关联度的图  #备注：： 绘图速度慢
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_similiraty.pdf")
getCorrelation(meth_kids,method ="spearman",plot=TRUE)
dev.off()

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth_kids,rowids=T))
head(data_raw_bin)
#mean DNA methylation level after filters
colnames(data_raw_bin)
data_raw_bin$YOUNG_Kids_mean<-apply(data_raw_bin[,1:10],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_Kids_mean<-apply(data_raw_bin[,11:20],1,function(x) mean(x,na.rm=T))
write.csv(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.csv",quote=F, row.names=T)
dim(data_raw_bin[,21:22])  # 711720 2

#画出样本聚类的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_clustering.pdf")
clusterSamples(meth_kids, dist="correlation", method="ward", plot=TRUE)
clusterSamples(na.omit(meth_kids), dist="correlation", method="ward", plot=TRUE)
clusterSamples(meth_kids, dist="correlation", method="complete", plot=TRUE)
clusterSamples(na.omit(meth_kids), dist="correlation", method="complete", plot=TRUE)

clusterSamples(meth_kids, dist="euclidean", method="ward", plot=TRUE)
clusterSamples(na.omit(meth_kids), dist="euclidean", method="ward", plot=TRUE)
clusterSamples(meth_kids, dist="euclidean", method="complete", plot=TRUE)
clusterSamples(na.omit(meth_kids), dist="euclidean", method="complete", plot=TRUE)
dev.off()

#画出样本主成分分析图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot.pdf")
PCASamples(meth_kids, screeplot=TRUE)
dev.off()

#画出样本主成分分析后的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot2.pdf")
PCASamples(meth_kids)
dev.off()

##PCA plot
data_raw<-as.data.frame(percMethylation(meth_kids,rowids=T))
data_raw<-data_raw[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",rownames(data_raw),invert=TRUE),]
data2<- t(data_raw)
set.seed(12)
md <- prep(data2, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 10)
summary(resPPCA) 

##PCA score:    PC1    PC2     PC3     PC4     PC5     PC6    PC7     PC8     PC9    PC10
#R2            0.1771 0.1020 0.08567 0.07545 0.06839 0.04931 0.04285 0.04127 0.03937 0.03824
#Cumulative R2 0.1771 0.2791 0.36480 0.44025 0.50865 0.55796 0.60081 0.64208 0.68145 0.71969
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCA.txt",quote=F, row.names=T, col.names=T)

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
rownames(huahua_meta)<-huahua_meta$library_code

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rownames(df)
df1<-merge(df,huahua_meta,by=0)
dim(df);dim(df1)
variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
df1$group2<-paste(df1$Age_group,df1$Sample.types,sep="_")
df1$group2 <-factor(df1$group2,level=c("YOUNG_kids","AMA_kids"))

head(df1)
str(df1)
df1<-df1[order(df1$group2,df1$analysis_name,decreasing = F),]
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/Kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot_2D_3D.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(7,10),rep(8,10))],
                                  pch = c(rep(16,20)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend("right",pch=16, legend = levels(df1$group2),col =ppCor[c(7:8)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(7:8)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(7:8)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()
####################
#for mother
sample_list<-c("10-1-in-9","16-1","17-1-2D","18-1","19-1-2G","20M-4F","6-1","7-1","8-1","E16M",
               "11-1","12-1-3G","13-1-1F","14-1","15-1","1M","2-1","3M","4-1","5-1")
meth_mother =reorganize(meth2,sample.ids=sample_list,treatment=c(rep(0,10),rep(1,10)))

#画出样本关联度的图  #备注：： 绘图速度慢
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_similiraty.pdf")
getCorrelation(meth_mother,method ="spearman",plot=TRUE)
dev.off()

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth_mother,rowids=T))
head(data_raw_bin)
#mean DNA methylation level after filters
colnames(data_raw_bin)
data_raw_bin$YOUNG_mother_mean<-apply(data_raw_bin[,1:10],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_mother_mean<-apply(data_raw_bin[,11:20],1,function(x) mean(x,na.rm=T))
write.csv(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.csv",quote=F, row.names=T)
dim(data_raw_bin[,21:22])  # 711720 2

#画出样本聚类的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_clustering.pdf")
clusterSamples(meth_mother, dist="correlation", method="ward", plot=TRUE)
clusterSamples(na.omit(meth_mother), dist="correlation", method="ward", plot=TRUE)
clusterSamples(meth_mother, dist="correlation", method="complete", plot=TRUE)
clusterSamples(na.omit(meth_mother), dist="correlation", method="complete", plot=TRUE)

clusterSamples(meth_mother, dist="euclidean", method="ward", plot=TRUE)
clusterSamples(na.omit(meth_mother), dist="euclidean", method="ward", plot=TRUE)
clusterSamples(meth_mother, dist="euclidean", method="complete", plot=TRUE)
clusterSamples(na.omit(meth_mother), dist="euclidean", method="complete", plot=TRUE)
dev.off()

#画出样本主成分分析图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot.pdf")
PCASamples(meth_mother, screeplot=TRUE)
dev.off()

#画出样本主成分分析后的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot2.pdf")
PCASamples(meth_mother)
dev.off()

##PCA plot
data_raw<-as.data.frame(percMethylation(meth_mother,rowids=T))
data_raw<-data_raw[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",rownames(data_raw),invert=TRUE),]
data2<- t(data_raw)
set.seed(12)
md <- prep(data2, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 10)
summary(resPPCA) 

##PCA score:    PC1    PC2     PC3     PC4     PC5     PC6    PC7     PC8     PC9    PC10
#R2            0.1466 0.1206 0.07473 0.06265 0.05489 0.04985 0.04777 0.04667 0.04452 0.04322
#Cumulative R2 0.1466 0.2672 0.34194 0.40459 0.45949 0.50934 0.55711 0.60378 0.64830 0.69152
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCA.txt",quote=F, row.names=T, col.names=T)

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
rownames(huahua_meta)<-huahua_meta$library_code

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rownames(df)
df1<-merge(df,huahua_meta,by=0)
dim(df);dim(df1)
variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
df1$group2<-paste(df1$Age_group,df1$Sample.types,sep="_")
df1$group2 <-factor(df1$group2,level=c("YOUNG_mother","AMA_mother"))

head(df1)
str(df1)
df1<-df1[order(df1$group2,df1$analysis_name,decreasing = F),]
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot_2D_3D.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(2,10),rep(1,10))],
                                  pch = c(rep(16,20)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend("right",pch=16, legend = levels(df1$group2),col =ppCor[c(7:8)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(2:1)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(2:1)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()
##########################
#for father
sample_list<-c("10-2-1D","16-2-4C","17-2-2E","18-2","19-2-2H","20F-4G","6-2","7F","8-2","E16F",
               "11-2","12-2-3H","13-2-1H","14-2","15-2","1F","2-2","3F","4-2","5-2")
meth_father =reorganize(meth2,sample.ids=sample_list,treatment=c(rep(0,10),rep(1,10)))

#画出样本关联度的图  #备注：： 绘图速度慢
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_similiraty.pdf")
getCorrelation(meth_father,method ="spearman",plot=TRUE)
dev.off()

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth_father,rowids=T))
head(data_raw_bin)
#mean DNA methylation level after filters
colnames(data_raw_bin)
data_raw_bin$YOUNG_father_mean<-apply(data_raw_bin[,1:10],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_father_mean<-apply(data_raw_bin[,11:20],1,function(x) mean(x,na.rm=T))
write.csv(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.csv",quote=F, row.names=T)
dim(data_raw_bin[,21:22])  # 711720 2

#画出样本聚类的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_clustering.pdf")
clusterSamples(meth_father, dist="correlation", method="ward", plot=TRUE)
clusterSamples(na.omit(meth_father), dist="correlation", method="ward", plot=TRUE)
clusterSamples(meth_father, dist="correlation", method="complete", plot=TRUE)
clusterSamples(na.omit(meth_father), dist="correlation", method="complete", plot=TRUE)

clusterSamples(meth_father, dist="euclidean", method="ward", plot=TRUE)
clusterSamples(na.omit(meth_father), dist="euclidean", method="ward", plot=TRUE)
clusterSamples(meth_father, dist="euclidean", method="complete", plot=TRUE)
clusterSamples(na.omit(meth_father), dist="euclidean", method="complete", plot=TRUE)
dev.off()

#画出样本主成分分析图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot.pdf")
PCASamples(meth_father, screeplot=TRUE)
dev.off()

#画出样本主成分分析后的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot2.pdf")
PCASamples(meth_father)
dev.off()

##PCA plot
data_raw<-as.data.frame(percMethylation(meth_father,rowids=T))
data_raw<-data_raw[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",rownames(data_raw),invert=TRUE),]
data2<- t(data_raw)
set.seed(12)
md <- prep(data2, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 10)
summary(resPPCA) 

##PCA score:    PC1    PC2     PC3     PC4     PC5     PC6    PC7     PC8     PC9    PC10
#R2            0.1134 0.09738 0.07566 0.07317 0.06406 0.05889 0.05806 0.05428 0.04566 0.04392
#Cumulative R2 0.1134 0.21080 0.28646 0.35964 0.42370 0.48259 0.54065 0.59494 0.64060 0.68452
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCA.txt",quote=F, row.names=T, col.names=T)

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
rownames(huahua_meta)<-huahua_meta$library_code

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rownames(df)
df1<-merge(df,huahua_meta,by=0)
dim(df);dim(df1)
variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
df1$group2<-paste(df1$Age_group,df1$Sample.types,sep="_")
df1$group2 <-factor(df1$group2,level=c("YOUNG_father","AMA_father"))

head(df1)
str(df1)
df1<-df1[order(df1$group2,df1$analysis_name,decreasing = F),]

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot_2D_3D.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(3,10),rep(5,10))],
                                  pch = c(rep(16,20)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend("right",pch=16, legend = levels(df1$group2),col =ppCor[c(7:8)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(3,5)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(3,5)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()


##evaluation for kids and parents
data_raw_bin<-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.txt", header = T,row.names= 1,check.names=F)
head(data_raw_bin)
dim(data_raw_bin)

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
head(huahua_meta)
rownames(huahua_meta)<-huahua_meta$library_code
#parental_samples<-rownames(huahua_meta[which(!(huahua_meta$Sample.types == "kids")),])
#data_raw_bin2<-data_raw_bin[,parental_samples]
##使用平均值进行绘图
data_region1<-as.data.frame(colMeans(data_raw_bin[,rownames(huahua_meta)],na.rm = T))
rownames(data_region1) 
length(na.omit(rownames(data_region1)));length(na.omit(rownames(huahua_meta)))
all(rownames(data_region1) %in% rownames(huahua_meta))

colnames(data_region1)<-"meth_mean"
data_region<-merge(data_region1,huahua_meta,by=0)
dim(data_region)#60 10
colnames(data_region)
data_region$Age_group <-factor(data_region$Age_group,level=c("YOUNG","AMA"))
data_region$Sample.types <-factor(data_region$Sample.types,level=c("kids","mother","father"))

head(data_region);dim(data_region)
#compare_means(meth_mean~Age_group, data=data_region,method = "t.test",group_by="group_age")
compare_means(meth_mean~Age_group,group.by = "Sample.types", data=data_region,method = "t.test")
#Sample.types .y.       group1 group2     p p.adj p.format p.signif method
#1 mother       meth_mean YOUNG  AMA    0.252  0.5  0.25     ns       T-test
#2 father       meth_mean YOUNG  AMA    0.305  0.5  0.31     ns       T-test
#3 kids         meth_mean YOUNG  AMA    0.144  0.43 0.14     ns       T-test
compare_means(meth_mean~Age_group,group.by = "Sample.types", data=data_region,method = "wilcox.test")
#Sample.types .y.       group1 group2     p p.adj p.format p.signif method
#1 mother       meth_mean YOUNG  AMA    0.436  0.94 0.44     ns       Wilcoxon
#2 father       meth_mean YOUNG  AMA    0.315  0.94 0.31     ns       Wilcoxon
#3 kids         meth_mean YOUNG  AMA    0.436  0.94 0.44     ns       Wilcoxon

global_meth_plot1 <- ggboxplot(data_region, x = "Age_group", y = "meth_mean",facet.by = "Sample.types",
                               color = "Age_group", palette = ppCor[2:1],
                               add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(60,90))
#  Add p-value
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.AMA_group_global_comparison_mean_all_site_by_6X_CpG_200bp_3L_8sample.pdf")
global_meth_plot1+stat_compare_means(method = "t.test",comparisons= list(c("YOUNG","AMA")))
global_meth_plot1+stat_compare_means(method = "wilcox.test",comparisons= list(c("YOUNG","AMA")))
#p+stat_compare_means(aes(label=..p.signif..), label.x = 1.5, label.y = 1.1)
global_meth_plot1+stat_compare_means(method = "t.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
global_meth_plot1+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
dev.off()

#grid.newpage(); #清空画板，开始画新图

#plot only for kids
global_meth_plot <- ggboxplot(data_region[which(data_region$Sample.types == "kids"),], x = "Age_group", y = "meth_mean",
                              color = "Age_group", palette = ppCor[2:1],
                              add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(70,85))+
  stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 80)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.kids_AMA_group_global_comparison_mean_all_site_by_6X_CpG_200bp_3L_8sample.pdf",global_meth_plot,width = 6, height = 6)


##plot for response
data_region2<-data_region[which(data_region$Sample.types == "kids"),]
head(data_region2)
aggregate(data_region2[,c("Age_group","meth_mean")],by=list(data_region2$Age_group),FUN=mean)
# Age_group meth_mean
#  YOUNG  75.30518
#   AMA   74.98727

global_meth_plot <- ggboxplot(data_region2, x = "Age_group", y = "meth_mean",color = "Age_group", palette = ppCor[c(2:1)],
                              add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(70,80))+
  stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 80)+
  stat_summary(fun=mean, geom="point", shape=18, size=4, col="pink", position = position_dodge(0.9))
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.kid_AMA_group_global_comparison_mean_all_site_by_6X_CpG_200bp_3L_8sample_response.pdf",global_meth_plot,width =5, height = 5)

#plot only for parents
data_region2<-data_region[which(data_region$Sample.types != "kids"),]
data_region2$group2<-paste(data_region2$Age_group,data_region2$Sample.types,sep="_")
data_region2$group2 <-factor(data_region2$group2,level=c("YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
my_comparisons <- list(c("YOUNG_mother","AMA_mother"),c("YOUNG_father","AMA_father"))
global_meth_plot <- ggboxplot(data_region2, x = "group2", y = "meth_mean",
                              color = "group2", palette = ppCor[c(2:1,3,5)],
                              add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(60,90))+
  stat_compare_means(method = "kruskal.test",aes(label=paste0(..method..,": "," p = ",..p.format..)), label.x = 1,label.y = 85)+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format") +labs(title = "wilcox.test")
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.parents_AMA_group_global_comparison_mean_all_site_by_6X_CpG_200bp_3L_8sample.pdf",global_meth_plot,width = 8, height = 8)

##plot for response
global_meth_plot <- ggboxplot(data_region2, x = "group2", y = "meth_mean",
                              color = "group2", palette = ppCor[c(2:1,3,5)],
                              add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(70,80))+
  stat_compare_means(method = "kruskal.test",aes(label=paste0(..method..,": "," p = ",..p.format..)), label.x = 1,label.y = 80)+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format") +#labs(title = "wilcox.test")+
  stat_summary(fun=mean, geom="point", shape=18, size=4, col="pink", position = position_dodge(0.9))
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.parents_AMA_group_global_comparison_mean_all_site_by_6X_CpG_200bp_3L_8sample_response.pdf",global_meth_plot,width =5, height = 5)

data_region
data_meth<-data_region[,c("analysis_name","Sample.types","Age_group","meth_mean")]
write.csv(data_meth, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_average_value_of_200bp__binfor_each_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.csv",quote=F, row.names=F) 
head(data_meth)
data_meth$group2<-paste(data_meth$Age_group,data_meth$Sample.types,sep="_")
data_meth$group2 <-factor(data_meth$group2,level=c("YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
aggregate(data_meth$meth_mean, by = list(data_meth$group2), FUN = mean) # 分组统计平均值。

# YOUNG_mother 74.33589
#   AMA_mother 73.97498
# YOUNG_father 74.07423
#   AMA_father 73.68770

#PCA plot for kids
kids_samples<-rownames(huahua_meta[which(huahua_meta$Sample.types == "kids"),])
data_raw_bin2<-data_raw_bin[,kids_samples]
dim(data_raw_bin2)
data2<- t(data_raw_bin2)
set.seed(12)
md <- prep(data2, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 10)
summary(resPPCA) 

##PCA score:    PC1    PC2     PC3     PC4     PC5     PC6    PC7     PC8     PC9    PC10
#R2            0.1771 0.1020 0.08567 0.07545 0.06839 0.04931 0.04285 0.04127 0.03937 0.03824
#Cumulative R2 0.1771 0.2791 0.36480 0.44025 0.50865 0.55796 0.60081 0.64208 0.68145 0.71969
#write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCA.txt",quote=F, row.names=T, col.names=T)

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rownames(df)
df1<-merge(df,huahua_meta,by=0)
dim(df);dim(df1)
variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
df1$group2<-paste(df1$Age_group,df1$Sample.types,sep="_")
df1$group2 <-factor(df1$group2,level=c("YOUNG_kids","AMA_kids"))

head(df1)
dim(df1)
df1<-df1[order(df1$group2,df1$analysis_name,decreasing = F),]

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot_2D_3D_2.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(2,10),rep(1,10))],
                                  pch = c(rep(16,20)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend("top",pch=16, legend = levels(df1$group2),col =ppCor[c(2:1)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(2:1)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(2:1)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()

#PCA plot for parents
parental_samples<-rownames(huahua_meta[which(!(huahua_meta$Sample.types == "kids")),])
data_raw_bin2<-data_raw_bin[,parental_samples]

data2<- t(data_raw_bin2)
set.seed(12)
md <- prep(data2, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 10)
summary(resPPCA) 

##PCA score:    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10
#R2            0.1059 0.07129 0.05405 0.04717 0.03956 0.03843 0.03685 0.02933 0.02777 0.02518
#Cumulative R2 0.1059 0.17720 0.23125 0.27841 0.31797 0.35640 0.39325 0.42258 0.45035 0.47553
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/parental_40_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCA.txt",quote=F, row.names=T, col.names=T)

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rownames(df)
df1<-merge(df,huahua_meta,by=0)
dim(df);dim(df1)
variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
df1$group2<-paste(df1$Age_group,df1$Sample.types,sep="_")
df1$group2 <-factor(df1$group2,level=c("YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))

head(df1)
dim(df1)
df1<-df1[order(df1$group2,df1$analysis_name,decreasing = F),]

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/parental_40_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot_2D_3D.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(2,10),rep(1,10),rep(3,10),rep(5,10))],
                                  pch = c(rep(16,40)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend("top",pch=16, legend = levels(df1$group2),col =ppCor[c(2:1,3,5)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(2:1,3,5)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group2))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(2:1,3,5)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()

###########################
#DMR_called
meth_200bp_three_list<-c(list(meth_kids),list(meth_mother),list(meth_father))
data_raw_bin<-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.txt", header = T,row.names= 1,check.names=F)
dim(data_raw_bin)
colnames(data_raw_bin)
#for kids
myDiff=calculateDiffMeth(meth_kids,num.cores=48)
myDiff@sample.ids
MYDIFF<-data.frame(myDiff)
head(MYDIFF)
rownames(MYDIFF)<-paste(MYDIFF$chr,MYDIFF$start,MYDIFF$end,sep = ".")
str(data_raw_bin[,c(myDiff@sample.ids)])
data_merge<-merge(data_raw_bin[,c(myDiff@sample.ids,"YOUNG_Kids_mean","AMA_Kids_mean")],MYDIFF[,5:7],by=0)
head(data_merge)
data_merge$meth_mean_def<-data_merge$AMA_Kids_mean-data_merge$YOUNG_Kids_mean
write.csv(data_merge, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv",quote=F, row.names=F) 
saveRDS(myDiff, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_DF_Called.rds")

#for mother
myDiff=calculateDiffMeth(meth_mother,num.cores=48)
myDiff@sample.ids
MYDIFF<-data.frame(myDiff)
head(MYDIFF)
rownames(MYDIFF)<-paste(MYDIFF$chr,MYDIFF$start,MYDIFF$end,sep = ".")
str(data_raw_bin[,c(myDiff@sample.ids)])
data_merge<-merge(data_raw_bin[,c(myDiff@sample.ids,"YOUNG_Mother_mean","AMA_Mother_mean")],MYDIFF[,5:7],by=0)
head(data_merge)
data_merge$meth_mean_def<-data_merge$AMA_Mother_mean-data_merge$YOUNG_Mother_mean
write.csv(data_merge, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv",quote=F, row.names=F) 
saveRDS(myDiff, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_DF_Called.rds")

#for father
myDiff=calculateDiffMeth(meth_father,num.cores=48)
myDiff@sample.ids
MYDIFF<-data.frame(myDiff)
head(MYDIFF)
rownames(MYDIFF)<-paste(MYDIFF$chr,MYDIFF$start,MYDIFF$end,sep = ".")
str(data_raw_bin[,c(myDiff@sample.ids)])
data_merge<-merge(data_raw_bin[,c(myDiff@sample.ids,"YOUNG_father_mean","AMA_father_mean")],MYDIFF[,5:7],by=0)
head(data_merge)
data_merge$meth_mean_def<-data_merge$AMA_father_mean-data_merge$YOUNG_father_mean
write.csv(data_merge, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv",quote=F, row.names=F) 
saveRDS(myDiff, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_DF_Called.rds")

###plot boxplot for global methlation level
