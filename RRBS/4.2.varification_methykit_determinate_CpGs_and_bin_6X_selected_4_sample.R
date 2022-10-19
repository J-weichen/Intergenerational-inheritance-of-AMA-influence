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
#myobj<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_myobj_6X_CpG.rds")
#length(myobj@treatment)
### filter sites with low (<1) or extremly high (>99.9% percentile) sequencing depth to reduce the effect of PCRs
depth<-6
#Mobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
#tiles_200=tileMethylCounts(Mobj,win.size=200,step.size=200,cov.bases =3)

#基因组划分区域进行计算 3 sample per groups
tiles_200 <-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_6X_tiles_200.rds") 
meth_200bp=methylKit::unite(tiles_200, destrand=FALSE,min.per.group=4L,mc.cores = 48)
write.table(as.data.frame(meth_200bp), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_200bp_CpG_6X_3site_bismark_4sample_covarage.txt",quote=F, row.names=F, col.names=T) 
#meth_200bp <- read.table(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_200bp_CpG_6X_3site_bismark_8sample_covarage.txt",header=T) 

#将样本聚合在一起
meth<-meth_200bp
head(meth)
meth@sample.ids
unique(meth$chr)
meth2<-meth[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",meth$chr,invert=TRUE),]
unique(meth2$chr)
meth2$chr<-factor(meth2$chr,levels =as.character(unique(meth2$chr)))
dim(meth2); dim(na.omit(meth2))
#[1]659395     64
#[1]514000     64

saveRDS(meth2, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth2.rds")

#meth2 <-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth2.rds")
#画出样本关联度的图  #备注：： 绘图速度慢
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_similiraty.pdf")
getCorrelation(meth2,method ="spearman",plot=TRUE)
dev.off()

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth2,rowids=T))
head(data_raw_bin)

colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/verification_AMA_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)

##huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
#head(huahua_meta)
#tiles_200_2[1:6,];meth_200bp_2[1:6,]
#rownames(huahua_meta)<-huahua_meta$library_code
#filtering and arrange data
evaluation.table<-data_raw_bin
evaluation.table[!(is.na(evaluation.table))]<- 1
evaluation.table[is.na(evaluation.table)]<- 0
bin200_detected_number<-colSums(evaluation.table)
bin200_detected_number<-data.frame(bin200_detected_number)

head(colData_whole);dim(colData_whole)
colnames(bin200_detected_number) <- "bin200_number"
rownames(bin200_detected_number);rownames(colData_whole)
bin200_detected_number2<-merge(bin200_detected_number,colData_whole,by=0)
head(bin200_detected_number2);dim(bin200_detected_number2)

bin200_detected_number2$group <-factor(bin200_detected_number2$group,level=c("Young","AMA"))
bin200_detected_number2$generation <-factor(bin200_detected_number2$generation,level=c("Kid","Mom"))
bin200_detected_number2$group_type <-factor(bin200_detected_number2$group_type,level=c("Young_Kid","AMA_Kid","Young_Mom","AMA_Mom"))
#bin200_detected_number2$Gender_RRBS<-factor(bin200_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin200_detected_number2)
bin200_detected_number2$sample<-as.character(bin200_detected_number2$sample)
bin200_detected_number2<-bin200_detected_number2[order(bin200_detected_number2$group_type,bin200_detected_number2$bin200_number,decreasing = T),]
bin200_detected_number2$sample<-factor(bin200_detected_number2$sample,levels=as.character(bin200_detected_number2$sample))

#plot
plot_bin200_detected_number0<-ggplot(data=bin200_detected_number2, mapping=aes(x=sample,y=bin200_number,fill=group_type))+geom_bar(stat="identity",width=0.8)
plot_bin200_detected_number1 <-plot_bin200_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin200_number), color="black", size=3,position=position_stack(1.01))+
  theme_classic()+labs(x="",y="Number of 200bp bin",title="Number of bin detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin200_detected_number1
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.number_of_varification_20_samples_200bp_bins_6X_CpG_3L_detected_in_per_group_3_sample.pdf",plot_bin200_detected_number1,width = 10, height = 8)

#mean DNA methylation level after filters
colnames(data_raw_bin)
data_raw_bin$AMA_Kids_mean<-apply(data_raw_bin[,1:5],1,function(x) mean(x,na.rm=T))
data_raw_bin$YOUNG_Kids_mean<-apply(data_raw_bin[,6:10],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_Mother_mean<-apply(data_raw_bin[,11:15],1,function(x) mean(x,na.rm=T))
data_raw_bin$YOUNG_Mother_mean<-apply(data_raw_bin[,16:20],1,function(x) mean(x,na.rm=T))

dim(data_raw_bin)#659395     24
write.table(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth.txt",quote=F,row.names=T)
write.csv(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth.csv",quote=F, row.names=T)
dim(data_raw_bin[,21:24])#659395      4
head(data_raw_bin[,21:24])

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_similarity_logisty.pdf")
p2<-ggplot(data = na.omit(data_raw_bin[,21:22]), mapping = aes(x = AMA_Kids_mean, y = YOUNG_Kids_mean)) + 
  stat_bin2d(bins = 300) + scale_fill_gradient(low = 'steelblue', high = 'darkred', limits = c(0,100), breaks = c(0,25,50,100))
p_plot1 <- p2+ stat_cor(method = "pearson", label.x = 10, label.y = 90)+xlim(0,100)+ xlab("AMA_Kids_mean") + 
  theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5,  hjust = 0.5))+ylab("YOUNG_Kids_mean") + 
  theme(axis.title.y = element_text(size = 16,face = "bold", vjust = 0.5, hjust = 0.5))+theme_bw()
print(p_plot1)

p4<-ggplot(data = na.omit(data_raw_bin[,23:24]), mapping = aes(x = AMA_Mother_mean, y = YOUNG_Mother_mean)) + 
  stat_bin2d(bins = 300) + scale_fill_gradient(low = 'steelblue', high = 'darkred', limits = c(0,100), breaks = c(0,25,50,100))
p_plot3 <- p4+ stat_cor(method = "pearson", label.x = 10, label.y = 90)+xlim(0,100)+ xlab("AMA_Mother_mean") + 
  theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5,  hjust = 0.5))+ylab("YOUNG_Mother_mean") + 
  theme(axis.title.y = element_text(size = 16,face = "bold", vjust = 0.5, hjust = 0.5))+theme_bw()
print(p_plot3)
dev.off()

#画出样本聚类的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_clustering.pdf",width=12, height=5)
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
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot.pdf")
PCASamples(meth2, screeplot=TRUE)
dev.off()

#画出样本主成分分析后的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot2.pdf")
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

##                PC1    PC2     PC3     PC4     PC5    PC6     PC7     PC8     PC9    PC10
#R2            0.1808 0.1281 0.05992 0.05685 0.05366 0.0504 0.04858 0.04652 0.04534 0.04038
#Cumulative R2 0.1808 0.3089 0.36880 0.42566 0.47932 0.5297 0.57829 0.62481 0.67014 0.71052
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCA.txt",quote=F, row.names=T, col.names=T)

#huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
#rownames(huahua_meta)<-huahua_meta$library_code
#dim(huahua_meta)
head(colData_whole)
plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
dim(df)
rownames(df);rownames(colData_whole)
df1<-merge(df,colData_whole,by=0)
dim(df1)#20 16

variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
df1$group_type <-factor(df1$group_type,level=c("Young_Kid","AMA_Kid","Young_Mom","AMA_Mom"))

head(df1)
str(df1)
df1<-df1[order(df1$group_type,df1$sample,decreasing = F),]
dim(df1)
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot_2D_3D.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(7,5),rep(8,5),rep(2,5),rep(1,5))],
                                  pch = c(rep(16,20)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
legend("right",pch=10, legend = levels(df1$group_type),col =ppCor[c(7:8,2:1)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group_type))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(7:8,2:1)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group_type))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(7:8,2:1)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()

####evaluation by kids mother and father
meth2<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth2.rds")
meth2@sample.ids
#########################
#for kids
sample_list<-c("Y1Q","Y2Q","Y4Q","Y5Q","Y6Q","A2Q","A3Q","A4Q","A5Q","A6Q")
meth_kids =reorganize(meth2,sample.ids=sample_list,treatment=c(rep(0,5),rep(1,5)))

#画出样本关联度的图  #备注：： 绘图速度慢
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_similiraty.pdf")
getCorrelation(meth_kids,method ="spearman",plot=TRUE)
dev.off()

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth_kids,rowids=T))
head(data_raw_bin)
#mean DNA methylation level after filters
colnames(data_raw_bin)
data_raw_bin$YOUNG_Kids_mean<-apply(data_raw_bin[,1:5],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_Kids_mean<-apply(data_raw_bin[,6:10],1,function(x) mean(x,na.rm=T))
write.csv(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_Kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth.csv",quote=F, row.names=T)
dim(data_raw_bin[,11:12])  # 711720 2

#画出样本聚类的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_Kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_clustering.pdf")
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
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_Kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot.pdf")
PCASamples(meth_kids, screeplot=TRUE)
dev.off()

#画出样本主成分分析后的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_Kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot2.pdf")
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

##PCA score:        PC1    PC2    PC3    PC4     PC5     PC6     PC7    PC8     PC9      PC10
#R2              0.2665 0.1179 0.1044 0.0978 0.08827 0.08471 0.08206 0.0805 0.07784 7.491e-05
#Cumulative R2   0.2665 0.3844 0.4887 0.5866 0.67482 0.75953 0.84158 0.9221 0.99993 1.000e+00
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_Kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCA.txt",quote=F, row.names=T, col.names=T)
head(colData_whole)

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rownames(df)
df1<-merge(df,colData_whole,by=0)
head(df1);dim(df1)

variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
head(df1)
df1<-df1[order(df1$group,df1$sample,decreasing = F),]
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_Kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot_2D_3D.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(7,5),rep(8,5))],
                                  pch = c(rep(16,10)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend("right",pch=16, legend = levels(df1$group2),col =ppCor[c(7:8)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group))+
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
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group))+
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
sample_list<-c("Y1M","Y2M","Y4M","Y5M","Y6M","A2M","A3M","A4M","A5M","A6M")
meth_mother =reorganize(meth2,sample.ids=sample_list,treatment=c(rep(0,5),rep(1,5)))

#画出样本关联度的图  #备注：： 绘图速度慢
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_similiraty.pdf")
getCorrelation(meth_mother,method ="spearman",plot=TRUE)
dev.off()

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth_mother,rowids=T))
head(data_raw_bin)
#mean DNA methylation level after filters
colnames(data_raw_bin)
data_raw_bin$YOUNG_mother_mean<-apply(data_raw_bin[,1:5],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_mother_mean<-apply(data_raw_bin[,6:10],1,function(x) mean(x,na.rm=T))
write.csv(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth.csv",quote=F, row.names=T)
dim(data_raw_bin[,11:12])  # 711720 2

#画出样本聚类的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_clustering.pdf")
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
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot.pdf")
PCASamples(meth_mother, screeplot=TRUE)
dev.off()

#画出样本主成分分析后的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot2.pdf")
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

##PCA score:      PC1    PC2    PC3    PC4    PC5    PC6    PC7     PC8     PC9      PC10
#R2            0.1444 0.1253 0.1184 0.1105 0.1059 0.1016 0.1007 0.09825 0.09483 0.0001387
#Cumulative R2 0.1444 0.2697 0.3881 0.4986 0.6045 0.7061 0.8068 0.90503 0.99986 1.0000000
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCA.txt",quote=F, row.names=T, col.names=T)

colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/verification_AMA_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rownames(df)
df1<-merge(df,colData_whole,by=0)

dim(df);dim(df1)
variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
head(df1)
str(df1)
df1<-df1[order(df1$group,df1$sample,decreasing = F),]
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot_2D_3D.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(2,5),rep(1,5))],
                                  pch = c(rep(16,10)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend("right",pch=16, legend = levels(df1$group2),col =ppCor[c(7:8)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(1:2)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(1:2)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()

##evaluation for kids and parents
data_raw_bin<-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth.txt", header = T,row.names= 1,check.names=F)
head(data_raw_bin)
dim(data_raw_bin)

colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/verification_AMA_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)
#parental_samples<-rownames(huahua_meta[which(!(huahua_meta$Sample.types == "kids")),])
#data_raw_bin2<-data_raw_bin[,parental_samples]
##使用平均值进行绘图
data_region1<-as.data.frame(colMeans(data_raw_bin[,rownames(colData_whole)],na.rm = T))
rownames(data_region1) 
length(na.omit(rownames(data_region1)));length(na.omit(rownames(colData_whole)))
all(rownames(data_region1) %in% rownames(colData_whole))

colnames(data_region1)<-"meth_mean"
data_region<-merge(data_region1,colData_whole,by=0)
dim(data_region)#60 10
colnames(data_region)
data_region$group <-factor(data_region$group,level=c("Young","AMA"))
data_region$generation <-factor(data_region$generation,level=c("Kid","Mom"))

head(data_region);dim(data_region)
#compare_means(meth_mean~Age_group, data=data_region,method = "t.test",group_by="group_age")
compare_means(meth_mean~group,group.by = "generation", data=data_region,method = "t.test")
#  generation .y.       group1 group2      p p.adj p.format p.signif method
#1 Mom        meth_mean AMA    Young  0.0350  0.07 0.035    *        T-test
#2 Kid        meth_mean AMA    Young  0.354   0.35 0.354    ns       T-test
compare_means(meth_mean~group,group.by = "generation", data=data_region,method = "wilcox.test")
#generation .y.       group1 group2      p p.adj p.format p.signif method  
#1 Mom        meth_mean AMA    Young  0.0556  0.11 0.056    ns       Wilcoxon
#2 Kid        meth_mean AMA    Young  0.841   0.84 0.841    ns       Wilcoxon

global_meth_plot1 <- ggboxplot(data_region, x = "group", y = "meth_mean",facet.by = "generation",
                               color = "group", palette = ppCor[1:2],add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(60,90))
#  Add p-value
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.AMA_group_global_comparison_mean_varification_20_samples_all_site_by_6X_CpG_200bp_3L_4sample_.pdf")
global_meth_plot1+stat_compare_means(method = "t.test",comparisons= list(c("Young","AMA")))
global_meth_plot1+stat_compare_means(method = "wilcox.test",comparisons= list(c("Young","AMA")))
#p+stat_compare_means(aes(label=..p.signif..), label.x = 1.5, label.y = 1.1)
global_meth_plot1+stat_compare_means(method = "t.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
global_meth_plot1+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
dev.off()

#grid.newpage(); #清空画板，开始画新图

#plot only for kids
global_meth_plot <- ggboxplot(data_region[which(data_region$generation == "Kid"),], x = "group", y = "meth_mean",
                              color = "group", palette = ppCor[1:2],
                              add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(65,90))+
  stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 80)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.kids_AMA_group_global_comparison_mean_varification_10_samples_all_site_by_6X_CpG_200bp_3L_4sample_.pdf",global_meth_plot,width = 6, height = 6)
#plot only for mom
global_meth_plot <- ggboxplot(data_region[which(data_region$generation == "Mom"),], x = "group", y = "meth_mean",
                              color = "group", palette = ppCor[1:2],
                              add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(65,90))+
  stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 80)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.Mom_AMA_group_global_comparison_mean_varification_10_samples_all_site_by_6X_CpG_200bp_3L_4sample_.pdf",global_meth_plot,width = 6, height = 6)

#PCA plot for kids
kids_samples<-rownames(colData_whole[which(colData_whole$generation == "Kid"),])
data_raw_bin2<-data_raw_bin[,kids_samples]

data2<- t(data_raw_bin2)
set.seed(12)
md <- prep(data2, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 10)
summary(resPPCA) 

##PCA score:   PC1    PC2    PC3    PC4     PC5     PC6     PC7    PC8     PC9      PC10
#R2            0.2665 0.1179 0.1044 0.0978 0.08827 0.08471 0.08206 0.0805 0.07784 7.491e-05
#Cumulative R2 0.2665 0.3844 0.4887 0.5866 0.67482 0.75953 0.84158 0.9221 0.99993 1.000e+00

#write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCA.txt",quote=F, row.names=T, col.names=T)

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rownames(df)

df1<-merge(df,colData_whole,by=0)
dim(df);dim(df1)
variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names

head(df1)
dim(df1)
df1<-df1[order(df1$group,df1$sample,decreasing = F),]

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot_2D_3D_2.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(2,5),rep(1,5))],
                                  pch = c(rep(16,10)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend("top",pch=16, legend = levels(df1$group2),col =ppCor[c(2:1)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(1:2)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(1:2)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()

#PCA plot for parents
mom_samples<-rownames(colData_whole[which(colData_whole$generation == "Mom"),])
data_raw_bin2<-data_raw_bin[,mom_samples]


data2<- t(data_raw_bin2)
set.seed(12)
md <- prep(data2, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 10)
summary(resPPCA) 

##PCA score:       PC1    PC2    PC3    PC4    PC5    PC6     PC7     PC8     PC9      PC10
#R2             0.1463 0.1249 0.1183 0.1120 0.1061 0.1024 0.09942 0.09713 0.09334 9.325e-05
#Cumulative R2  0.1463 0.2712 0.3895 0.5015 0.6076 0.7100 0.80944 0.90656 0.99991 1.000e+00
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mom_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCA.txt",quote=F, row.names=T, col.names=T)

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rownames(df)
df1<-merge(df,colData_whole,by=0)
dim(df);dim(df1)
variation<-data.frame(summary(resPPCA))
rownames(df1)<-df1$Row.names
head(df1);dim(df1)
df1<-df1[order(df1$group,df1$sample,decreasing = F),]

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mom_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_PCAscreeplot_2D_3D.pdf")
plot3d <- with(df1, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(1,5),rep(2,5))],
                                  pch = c(rep(16,10)),
                                  cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend("top",pch=16, legend = levels(df1$group2),col =ppCor[c(2:1,3,5)],horiz =F)
#legend("right",pch=16,legend=levels(df1$Age_group), col =ppCor[2:1],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=df1, mapping=aes(x=PC1,y=PC2,colour = group))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(1:2)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=df1, mapping=aes(x=PC1,y=PC3,colour = group))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=4)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(1:2)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()

###########################
#DMR_called
meth_200bp_two_list<-c(list(meth_kids),list(meth_mother))
data_raw_bin<-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth.txt", header = T,row.names= 1,check.names=F)
dim(data_raw_bin)#602367     24
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
table(data_merge$pvalue<0.05)
#FALSE   TRUE 
#424582 177785 
table(data_merge$qvalue<0.05)
#FALSE   TRUE 
#489063 113304
write.csv(data_merge, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth_PVALUE.csv",quote=F, row.names=F) 
saveRDS(myDiff, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth_DF_Called.rds")

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
table(data_merge$pvalue<0.05)
#FALSE   TRUE 
#447241 155126 
table(data_merge$qvalue<0.05)
#FALSE   TRUE 
#518409  83958 
write.csv(data_merge, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth_PVALUE.csv",quote=F, row.names=F) 
saveRDS(myDiff, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth_DF_Called.rds")
