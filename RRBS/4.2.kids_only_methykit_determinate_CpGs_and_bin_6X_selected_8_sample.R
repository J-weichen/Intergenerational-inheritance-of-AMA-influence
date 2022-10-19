#loading package
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
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

grid.newpage(); #清空画板，开始画新图

AMA_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Kids")
YOUNG_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Kids")

AMA_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Kids/",AMA_Kids_files))
YOUNG_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Kids/",YOUNG_Kids_files))

AMA_Kids_names<- lapply(strsplit(AMA_Kids_files,"_"), function(x) x[1])
YOUNG_Kids_names<- lapply(strsplit(YOUNG_Kids_files,"_"), function(x) x[1])
#covarage 6
myobj=methRead(c(YOUNG_Kids_files2,AMA_Kids_files2),
               sample.id=c(YOUNG_Kids_names,AMA_Kids_names),
               assembly="hg38",pipeline="bismarkCoverage",#default:'amp'
               treatment=c(rep(0,length(YOUNG_Kids_names)),rep(1,length(AMA_Kids_names))),
               context="CpG", resolution = "base", mincov = 6)
dim(as.data.frame(myobj[1]))
saveRDS(myobj, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_single_kids_20_samples_myobj_6X_CpG.rds")

#myobj<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_single_kids_20_samples_myobj_6X_CpG.rds")
#length(myobj@treatment)
### filter sites with low (<1) or extremly high (>99.9% percentile) sequencing depth to reduce the effect of PCRs
depth<-6
Mobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
tiles_200=tileMethylCounts(Mobj,win.size=200,step.size=200,cov.bases =3)
kids_meth_200bp=methylKit::unite(tiles_200, destrand=FALSE,min.per.group=8L,mc.cores = 48)
write.table(as.data.frame(kids_meth_200bp), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_200bp_CpG_6X_3site_bismark_8sample_covarage.txt",quote=F, row.names=F, col.names=T) 

kids_meth_200bp@sample.ids
unique(kids_meth_200bp$chr)
meth_kids<-kids_meth_200bp[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",kids_meth_200bp$chr,invert=TRUE),]
unique(meth_kids$chr)
meth_kids$chr<-factor(meth_kids$chr,levels =as.character(unique(meth2$chr)))
dim(meth_kids); dim(na.omit(meth_kids))#756666 #566219 
meth_kids@sample.ids
#########################

#画出样本关联度的图  #备注：： 绘图速度慢
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_similiraty.pdf")
getCorrelation(meth_kids,method ="spearman",plot=TRUE)
dev.off()

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth_kids,rowids=T))
head(data_raw_bin)
#mean DNA methylation level after filters
colnames(data_raw_bin)
data_raw_bin$YOUNG_Kids_mean<-apply(data_raw_bin[,1:10],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_Kids_mean<-apply(data_raw_bin[,11:20],1,function(x) mean(x,na.rm=T))

write.csv(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.csv",quote=F, row.names=T)
dim(data_raw_bin[,21:22])  # 711720 2

#画出样本聚类的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_clustering.pdf")
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
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot.pdf")
PCASamples(meth_kids, screeplot=TRUE)
dev.off()

#画出样本主成分分析后的图
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot2.pdf")
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
#R2            0.2255 0.1128 0.09288 0.06323 0.05841 0.0535 0.04638 0.03916 0.03597 0.03429
#Cumulative R2 0.2255 0.3383 0.43118 0.49441 0.55282 0.6063 0.65271 0.69186 0.72784 0.76213
write.table((as.data.frame(summary(resPPCA))), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCA.txt",quote=F, row.names=T, col.names=T)

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
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_PCAscreeplot_2D_3D.pdf")
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

#######################
#DMR_called
data_raw_bin<-read.csv(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.csv", header = T,row.names= 1,check.names=F)
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
write.csv(data_merge, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv",quote=F, row.names=F) 
saveRDS(myDiff, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/single_kids_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_DF_Called.rds")
