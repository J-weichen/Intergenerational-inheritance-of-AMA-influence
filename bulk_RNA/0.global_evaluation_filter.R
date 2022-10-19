##second time work directly goto 100 line

rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

library(scales)
library(ggsci)
library(dplyr)
library(ggpubr)
library(DESeq2) 
library(RColorBrewer)
library(amap)
library(scatterplot3d)
library(pcaMethods)
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)
pal2<-pal_jama("default",alpha = 1)(7)
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)
pal5 <- pal_npg("nrc", alpha=0.5)(10)
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor_all)
#step 0 making  count data
files_AMA_mom<- list.files("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts/AMA_mom")
files_Young_mom<- list.files("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts/Young_mom")
files_AMA_kid<- list.files("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts/AMA_kid")
files_Young_kid<- list.files("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts/Young_kid")

files_AMA_mom2<-paste0("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts/AMA_mom/",files_AMA_mom)
files_Young_mom2<- paste0("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts/Young_mom/",files_Young_mom)
files_AMA_kid2<- paste0("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts/AMA_kid/",files_AMA_kid)
files_Young_kid2<- paste0("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts/Young_kid/",files_Young_kid)

sample_name<-c(files_AMA_mom,files_Young_mom,files_AMA_kid,files_Young_kid)
files<-c(files_AMA_mom2,files_Young_mom2,files_AMA_kid2,files_Young_kid2)

read.delim(paste0("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts/AMA_mom/",files_AMA_mom[1]),header = F, nrow=2,row.names = 1,sep = "\t")
sample1 <- read.delim(files[1],header = F,row.names = 1,sep = "\t",col.names = c("gene_id",substring(sample_name[1], 1, nchar(files_AMA_mom[1])-19)))
head(sample1)
#read.delim 是 read.table 的变形，read_delim 比 read.table要快
count.table <- data.frame(sample1)
head(count.table)
#options(stringsAsFactors = FALSE)
for ( sample_order in 2:length(files)){
  #sample_order <-2
  f_name<- files[sample_order]
  names<-sample_name[sample_order]
  print(f_name)
  column <- read.delim(f_name,header=F,row.names = 1,sep = "\t",col.names = c("gene_id",substring(names, 1, nchar(names)-19)))
  head(column)
  count.table <- cbind(count.table,column)
}
head(count.table) 
count.table<-count.table[which(rowSums(count.table)>0),]
str(count.table);head(count.table);dim(count.table)#20985    26
write.table(count.table,file="/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/file1_All_filter_Verification_AMA_versus_Young_rawcount.xls",sep="\t", quote=F, row.names=T,col.names=T)

#step one:read data
count.table2<-read.table(file="/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/file1_All_filter_Verification_AMA_versus_Young_rawcount.xls",sep="\t", row.names=1,header =T)
head(count.table2)

#担心性别影响因此去除性染染色体上基因
data_chromX<-read.table("/mnt/data/chenwei/huahua/0.hua_script/0.RNA_script/ChromX_genelist.txt",head = F,sep = "\t"	)
data_chromY<-read.table("/mnt/data/chenwei/huahua/0.hua_script/0.RNA_script/ChromY_genelist.txt",head = F,sep = "\t")
head(data_chromX)
geneX<-unlist(as.character(unique(data_chromX$V13)))
geneY<-unlist(as.character(unique(data_chromY$V13)))
geneXY<-c(geneX,geneY);length(geneXY)#1246
count.table_XY_removed<-count.table2[-which(rownames(count.table2) %in% geneXY),]
head(count.table_XY_removed)
dim(count.table_XY_removed);dim(count.table2)#20230 20985
write.table(count.table_XY_removed,file="/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/file2_All_filter_Verification_AMA_versus_Young_rawcount_XY_removed.xls",sep="\t", quote=F, row.names=T,col.names=T)

#colData可以直接提供也可以自己构建
samplenames2<-colnames(count.table2)#将要分析的亚集列名字赋予此
group<-ifelse(grepl("A",samplenames2),"AMA","Young")
group<-as.factor(group)
generation<-ifelse(grepl("M",samplenames2),"Mom","Kid")
generation<-as.factor(generation)
sample<-as.factor(samplenames2)
colData_whole <- data.frame(sample = sample, group=group,generation=generation)
#最后构建数据:
rownames(colData_whole)<-colData_whole$sample
colData_whole$group_type<-paste0(colData_whole$group,"_",colData_whole$generation)
colData_whole$type<-colData_whole$generation
colData_whole$type<-ifelse(colData_whole$sample %in% c("A2Q","A6Q","Y1Q","A7Q"),"Boy",ifelse(grepl("M",colData_whole$sample),"Mom","Girl"))
colData_whole$time<-colData_whole$type
colData_whole$time<-ifelse(colData_whole$sample %in% c("Y7Q","Y8Q","A7Q","Y7M","Y8M","A7M"),"Second","First")
head(colData_whole)
write.table(colData_whole,file="/mnt/data/chenwei/huahua/file3_All_filter_verification_AMA_analysis_metadata.csv",sep=",", quote=F, row.names=T,col.names=T)
colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/file3_All_filter_verification_AMA_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)

#step2 reading count and colData and select used files
count.table_XY_removed<-read.table(file="/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/file2_All_filter_Verification_AMA_versus_Young_rawcount_XY_removed.xls",sep="\t", row.names=1,header =T)
colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/file3_All_filter_verification_AMA_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole);head(count.table_XY_removed)
colData_whole<-colData_whole[which(!(colData_whole$sample %in% c("Y5Q","Y6Q","A2Q","Y5M","Y6M","A2M"))),]
count.table_XY_removed<-count.table_XY_removed[,which(!(colnames(count.table_XY_removed) %in% c("Y5Q","Y6Q","A2Q","Y5M","Y6M","A2M")))]
head(colData_whole);head(count.table_XY_removed)

#绘图
tag_name<-"All_sample" #test line
colData<-colData_whole
expression_matrix<-count.table_XY_removed[,as.character(colData$sample)]

#step two : evaluation gene detective number
evaluation.table<-expression_matrix
evaluation.table[evaluation.table !=0]<- 1
gene_detected_number<-colSums(evaluation.table)
gene_detected_number<-data.frame(gene_detected_number)
colnames(gene_detected_number) <- "gene_number"
gene_detected_number2<-merge(gene_detected_number,colData,by=0)
gene_detected_number2$group_type<-factor(gene_detected_number2$group_type,levels = c("AMA_Mom","Young_Mom","AMA_Kid","Young_Kid"))
gene_detected_number2$Row.names<-as.character(gene_detected_number2$Row.names)
gene_detected_number2<-gene_detected_number2[order(gene_detected_number2$group_type,decreasing = T),]
gene_detected_number2$Row.names<-factor(gene_detected_number2$Row.names,levels=as.character(gene_detected_number2$Row.names))
head(gene_detected_number2)

plot_gene_detected_number<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number,fill=Row.names))+geom_bar(stat="identity",width=0.8)

plot_gene_detected_number2<-plot_gene_detected_number+scale_fill_manual(values=ppCor_all)+
  geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
  theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_gene_detected_number2

plot_gene_detected_number0<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number,fill=group))+geom_bar(stat="identity",width=0.8)
plot_gene_detected_number3 <-plot_gene_detected_number0+scale_fill_manual(values=ppCor)+
  geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
  theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_gene_detected_number3

plot_gene_detected_number12<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number,fill=group_type))+geom_bar(stat="identity",width=0.8)
plot_gene_detected_number4<-plot_gene_detected_number12+scale_fill_manual(values=ppCor[c(10:7)])+
  geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
  theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_gene_detected_number4

plot_gene_detected_number13<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number,fill=type))+geom_bar(stat="identity",width=0.8)
plot_gene_detected_number5<-plot_gene_detected_number13+scale_fill_manual(values=ppCor[c(10:7)])+
  geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
  theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_gene_detected_number5

count_plot<-ggarrange(plot_gene_detected_number2, plot_gene_detected_number3,plot_gene_detected_number4,plot_gene_detected_number5,
                      labels = c("A", "B", "C", "D"),ncol = 1, nrow = 4)
ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_gene_detected_in_AMA_Young_RNA_library.pdf"),count_plot,width = 200, height =240, units = "mm")
#}
#evaluation global effect
dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix, colData = colData, design = ~ group)
#data pre_cleaning
dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
dds2 <- estimateSizeFactors(dds1)
#data transform:sample number < 30 -> rlog
#show transform result
rld <- rlog(dds2, blind = FALSE);vsd <- vst(dds2, blind = FALSE)
colnames(rld)
df <- bind_rows(
  as_data_frame(log2(counts(dds2, normalized=TRUE)[, 2:3]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 2:3]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 2:3]) %>% mutate(transformation = "vst"))
colnames(df)[1:2] <- c("sample2","sample3")  
trans_plot<-ggplot(df, aes(x = sample2, y = sample3)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)
ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_trans_three_type_in_AMA_Young_RNA_library.pdf"),trans_plot,width = 620, height =200, units = "mm")

#样本距离： RNA-Seq分析第一步通常是评估样本间的总体相似度
anno1 <- as.data.frame(colData(dds2)[, c("group","generation")])
#count not normalization
sampleDists <- dist(t(log2(counts(dds2,normalized=F)+1)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_distance_count_not_normal_in_AMA_Young_RNA_library.pdf"),width = 10, height =10)
p1<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "Log2_trans_no_normalization")
print(p1)
dev.off()
#count  normalization
sampleDists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_distance_count_normal_in_AMA_Young_RNA_library.pdf"),width = 10, height =10)
p2<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "Log2_trans_normal")
print(p2)
dev.off()
#count rld transform
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_distance_count_rld_in_AMA_Young_RNA_library.pdf"),width = 10, height =10)
p3<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "rld_trans")
print(p3)
dev.off()
#count vsd transform
sampleDists <- dist(t(assay(vsd)))
sampleDists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_distance_count_vsd_in_AMA_Young_RNA_library.pdf"),width = 10, height =10)
p4<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "vsd_trans")
print(p4)
dev.off()

#cluster
dists <- dist(t(log2(counts(dds2, normalized=F)+1)),method = "euclidean") 
hc <- hcluster(dists, method="pearson")
row_dend <- as.dendrogram(hc)
hc$order
lable_cor<-rep(ppCor[2],20)
lable_cor[which(hc$order %in% c(1:5,11:15))] <- ppCor[1]
labels_colors(row_dend) <- lable_cor
pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_cluster_not_normal_in_AMA_Young_RNA_library.pdf"),width = 10, height =6)
p5<-plot(row_dend, main = "samples similarity \n count not normalization (pearson)")
print(p5)
dev.off()
#count normalization
dists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)),method = "euclidean") 
hc <- hcluster(dists, method="pearson")
#hc <- hcluster(dists, method="spearman")
row_dend <- as.dendrogram(hc)
lable_cor<-rep(ppCor[2],26)
lable_cor[which(hc$order %in% c(1:6,14:19))] <- ppCor[1]
labels_colors(row_dend) <- lable_cor
pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_cluster_normal_in_AMA_Young_RNA_library.pdf"),width = 10, height =6)
p6<-plot(row_dend, main = "samples similarity \n count normalization (pearson)")
print(p6)
dev.off()
#rld
dists <- dist(t(assay(rld)),method = "euclidean") 
hc <- hcluster(dists, method="pearson")
row_dend <- as.dendrogram(hc)
lable_cor<-rep(ppCor[2],26)
lable_cor[which(hc$order %in% c(1:6,14:19))] <- ppCor[1]
labels_colors(row_dend) <- lable_cor
pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_cluster_rld_in_AMA_Young_RNA_library.pdf"),width = 10, height =6)
p7<-plot(row_dend, main = "samples similarity \n rlg (pearson)")
print(p7)
dev.off()
#vsd
dists <- dist(t(assay(vsd)),method = "euclidean") 
hc <- hcluster(dists, method="pearson")
row_dend <- as.dendrogram(hc)
lable_cor<-rep(ppCor[2],26)
lable_cor[which(hc$order %in% c(1:6,14:19))] <- ppCor[1]
labels_colors(row_dend) <- lable_cor
pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_cluster_vsd_in_AMA_Young_RNA_library.pdf"),width = 10, height =6)
p8<-plot(row_dend, main = "samples similarity \n vsd (pearson)")
print(p8)
dev.off()

#另一个可视化样本-样本距离的方法就是主成分分析。
#DESeq2提供了专门的方法用于作图
plotPCA(rld,intgroup=c("group"))
#plotPCA(rld,intgroup=c("Gender_RRBS"))

rld_PCA<-plotPCA(rld,intgroup=c("group","generation"),returnData = TRUE)
percentVar<-round(100*attr(rld_PCA,"percentVar"),1)
rld_PCA$group<-factor(rld_PCA$group,levels=c("AMA:Kid","AMA:Mom","Young:Kid","Young:Mom"))

plot_rld_PCA<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour = group,shape= group))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  labs(title ="Distribution of samples \n rld Deseq2 inner method")+scale_color_manual(values=c(ppCor))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_rld_PCA
ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_rld_PCA_in_AMA_Young_RNA_library.pdf"),plot_rld_PCA,width = 10, height =10)


# perform PCA analysis
rlog_Mat <- assay(rld)
PCA_data1 <- t(rlog_Mat)
set.seed(19921010)
md <- prep(PCA_data1, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
summary(resPPCA) 
#Importance of component(s):
#                 PC1    PC2    PC3    PC4     PC5
#R2            0.5674 0.1009 0.06117 0.05009 0.03397
#Cumulative R2 0.5674 0.6683 0.72947 0.77957 0.81354

variation<-data.frame(summary(resPPCA))

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rld_PCA<-merge(colData,df,by=0)
rld_PCA$group_type<-factor(rld_PCA$group_type,levels=c("AMA_Kid","AMA_Mom","Young_Kid","Young_Mom"))

#rld_PCA$Gender_RRBS<-factor(rld_PCA$Gender_RRBS,levels=c("Female","Male"))
plot_rld_PCA2<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour =group_type,shape= group_type))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples \n rld")+scale_color_manual(values=c(ppCor[1:4]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_rld_PCA2
ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_rld_PCA_self_cal_in_AMA_Young_RNA_library.pdf"),plot_rld_PCA2,width = 10, height =10)

#plot3d <- with(rld_PCA,scatterplot3d(PC1, PC2, PC3, color =c(rep(ppCor[2],5), rep(ppCor[1],5)), pch = c(rep(8,5), rep(16,5)),cex.symbols = 1.2, font.lab = 2, font.axis = 2))
#legend(plot3d$xyz.convert(-15,0,10),pch=16,legend=levels(df$group), col =ppCor[1:2],cex = rel(1.1), bty = 'n',yjust=0.5, xjust = 0, horiz = F)

# perform PCA analysis for count not normalization
count_Mat <- log2(counts(dds2, normalized=F)+1)
PCA_data1 <- t(count_Mat)
set.seed(19921010)
md <- prep(PCA_data1, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
summary(resPPCA) 
#Importance of component(s):
#                 PC1    PC2    PC3     PC4     PC5
#R2             0.6428 0.06457 0.05142 0.03214 0.02434
#Cumulative R2  0.6428 0.70735 0.75877 0.79091 0.81525

variation<-data.frame(summary(resPPCA))
plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
count_PCA<-merge(colData,df,by=0)
count_PCA$group_type<-factor(count_PCA$group_type,levels=c("AMA_Kid","AMA_Mom","Young_Kid","Young_Mom"))
#count_PCA$Gender_RRBS<-factor(count_PCA$Gender_RRBS,levels=c("Female","Male"))
plot_not_normal_PCA<-ggplot(data=count_PCA, mapping=aes(x=PC1,y=PC2,colour = group_type,shape= group_type))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples \n count not normalization")+scale_color_manual(values=c(ppCor[1:4]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_not_normal_PCA
ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_not_normal_PCA_in_AMA_Young_RNA_library.pdf"),plot_not_normal_PCA,width = 10, height =10)

# perform PCA analysis for count  normalization
count_Mat <- log2(counts(dds2, normalized=T)+1)
PCA_data1 <- t(count_Mat)
set.seed(19921010)
md <- prep(PCA_data1, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
summary(resPPCA) 
#Importance of component(s):
#                 PC1    PC2    PC3    PC4     PC5
#R2            0.4772 0.09482 0.06238 0.05015 0.04215
#Cumulative R2 0.4772 0.57207 0.63444 0.68459 0.72674

variation<-data.frame(summary(resPPCA))
#variation$PC1[1]
#variation$PC2[1]
plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
count_PCA<-merge(colData,df,by=0)
count_PCA$group_type<-factor(count_PCA$group_type,levels=c("AMA_Kid","AMA_Mom","Young_Kid","Young_Mom"))
#count_PCA$Gender_RRBS<-factor(count_PCA$Gender_RRBS,levels=c("Female","Male"))
plot_normal_PCA<-ggplot(data=count_PCA, mapping=aes(x=PC1,y=PC2,colour = group_type,shape= group_type))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples \n count normalization")+scale_color_manual(values=c(ppCor[1:4]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_normal_PCA
ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_normal_PCA_in_AMA_Young_RNA_library.pdf"),plot_normal_PCA,width = 10, height =10)



#母子分开
type<-c("Mom","Kid")
for (tag_name in type){
  # tag_name<-"Kid" #test line
  
  colData<-colData_whole[which(colData_whole$generation==tag_name),]
  expression_matrix<-count.table_XY_removed[,as.character(colData$sample)]
  
  #step two : evaluation gene detective number
  evaluation.table<-expression_matrix
  evaluation.table[evaluation.table !=0]<- 1
  gene_detected_number<-colSums(evaluation.table)
  gene_detected_number<-data.frame(gene_detected_number)
  colnames(gene_detected_number) <- "gene_number"
  gene_detected_number2<-merge(gene_detected_number,colData,by=0)
  gene_detected_number2$group_type<-factor(gene_detected_number2$group_type,levels = paste0(c("AMA_","Young_"),tag_name))
  gene_detected_number2$Row.names<-as.character(gene_detected_number2$Row.names)
  gene_detected_number2<-gene_detected_number2[order(gene_detected_number2$group_type,decreasing = T),]
  gene_detected_number2$Row.names<-factor(gene_detected_number2$Row.names,levels=as.character(gene_detected_number2$Row.names))
  head(gene_detected_number2)
  
  plot_gene_detected_number<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number,fill=Row.names))+geom_bar(stat="identity",width=0.8)
  plot_gene_detected_number2<-plot_gene_detected_number+scale_fill_manual(values=ppCor_all)+
    geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
    theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
    theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
          axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
  plot_gene_detected_number2
  
  plot_gene_detected_number0<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number,fill=group))+geom_bar(stat="identity",width=0.8)
  plot_gene_detected_number3 <-plot_gene_detected_number0+scale_fill_manual(values=ppCor)+
    geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
    theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
    theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
          axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
  plot_gene_detected_number3
  
  plot_gene_detected_number12<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number,fill=group_type))+geom_bar(stat="identity",width=0.8)
  plot_gene_detected_number4<-plot_gene_detected_number12+scale_fill_manual(values=ppCor[c(10:7)])+
    geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
    theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
    theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
          axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
  plot_gene_detected_number4
  
  count_plot<-ggarrange(plot_gene_detected_number2, plot_gene_detected_number3,plot_gene_detected_number4,
                        labels = c("A", "B", "C"),ncol = 1, nrow = 3)
  ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_gene_detected_in_AMA_Young_RNA_library.pdf"),count_plot,width = 200, height =200, units = "mm")
  #}
  #evaluation global effect
  dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix, colData = colData, design = ~ group)
  #data pre_cleaning
  dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
  dds2 <- estimateSizeFactors(dds1)
  #data transform:sample number < 30 -> rlog
  #show transform result
  rld <- rlog(dds2, blind = FALSE);vsd <- vst(dds2, blind = FALSE)
  colnames(rld)
  df <- bind_rows(
    as_data_frame(log2(counts(dds2, normalized=TRUE)[, 2:3]+1)) %>% mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(rld)[, 2:3]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 2:3]) %>% mutate(transformation = "vst"))
  colnames(df)[1:2] <- c("sample2","sample3")  
  trans_plot<-ggplot(df, aes(x = sample2, y = sample3)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)
  ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_trans_three_type_in_AMA_Young_RNA_library.pdf"),trans_plot,width = 620, height =200, units = "mm")
  
  #样本距离： RNA-Seq分析第一步通常是评估样本间的总体相似度
  anno1 <- as.data.frame(colData(dds2)[, c("group","generation")])
  #count not normalization
  sampleDists <- dist(t(log2(counts(dds2,normalized=F)+1)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
  pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_distance_count_not_normal_in_AMA_Young_RNA_library.pdf"),width = 10, height =10)
  p1<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "Log2_trans_no_normalization")
  print(p1)
  dev.off()
  #count  normalization
  sampleDists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
  pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_distance_count_normal_in_AMA_Young_RNA_library.pdf"),width = 10, height =10)
  p2<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "Log2_trans")
  print(p2)
  dev.off()
  #count rld transform
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
  pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_distance_count_rld_in_AMA_Young_RNA_library.pdf"),width = 10, height =10)
  p3<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "rld_trans")
  print(p3)
  dev.off()
  #count vsd transform
  sampleDists <- dist(t(assay(vsd)))
  sampleDists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
  pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_distance_count_vsd_in_AMA_Young_RNA_library.pdf"),width = 10, height =10)
  p4<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "vsd_trans")
  print(p4)
  dev.off()
  
  #cluster
  dists <- dist(t(log2(counts(dds2, normalized=F)+1)),method = "euclidean") 
  hc <- hcluster(dists, method="pearson")
  row_dend <- as.dendrogram(hc)
  hc$order
  lable_cor<-rep(ppCor[2],13)
  lable_cor[which(hc$order %in% 1:6)] <- ppCor[1]
  labels_colors(row_dend) <- lable_cor
  pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_cluster_not_normal_in_AMA_Young_RNA_library.pdf"),width = 10, height =6)
  p5<-plot(row_dend, main = "samples similarity \n count not normalization (pearson)")
  print(p5)
  dev.off()
  #count normalization
  dists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)),method = "euclidean") 
  hc <- hcluster(dists, method="pearson")
  row_dend <- as.dendrogram(hc)
  lable_cor<-rep(ppCor[2],13)
  lable_cor[which(hc$order %in% 1:6)] <- ppCor[1]
  labels_colors(row_dend) <- lable_cor
  pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_cluster_normal_in_AMA_Young_RNA_library.pdf"),width = 10, height =6)
  p6<-plot(row_dend, main = "samples similarity \n count normalization (pearson)")
  print(p6)
  dev.off()
  #rld
  dists <- dist(t(assay(rld)),method = "euclidean") 
  hc <- hcluster(dists, method="pearson")
  row_dend <- as.dendrogram(hc)
  lable_cor<-rep(ppCor[2],13)
  lable_cor[which(hc$order %in% 1:6)] <- ppCor[1]
  labels_colors(row_dend) <- lable_cor
  pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_cluster_rld_in_AMA_Young_RNA_library.pdf"),width = 10, height =6)
  p7<-plot(row_dend, main = "samples similarity \n rlg (pearson)")
  print(p7)
  dev.off()
  #vsd
  dists <- dist(t(assay(vsd)),method = "euclidean") 
  hc <- hcluster(dists, method="pearson")
  row_dend <- as.dendrogram(hc)
  lable_cor<-rep(ppCor[2],13)
  lable_cor[which(hc$order %in% 1:6)] <- ppCor[1]
  labels_colors(row_dend) <- lable_cor
  pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_cluster_vsd_in_AMA_Young_RNA_library.pdf"),width = 10, height =6)
  p8<-plot(row_dend, main = "samples similarity \n vsd (pearson)")
  print(p8)
  dev.off()
  
  #另一个可视化样本-样本距离的方法就是主成分分析。
  #DESeq2提供了专门的方法用于作图
  plotPCA(rld,intgroup=c("group"))
  #plotPCA(rld,intgroup=c("Gender_RRBS"))
  
  rld_PCA<-plotPCA(rld,intgroup=c("group"),returnData = TRUE)
  percentVar<-round(100*attr(rld_PCA,"percentVar"),1)
  rld_PCA$group<-factor(rld_PCA$group,levels=c("AMA","Young"))
  #rld_PCA$Batch_RNA<-factor(rld_PCA$Batch_RNA,levels=c("First","Second"))
  
  plot_rld_PCA<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour = group,shape= group))+
    geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
    xlab(paste0("PC1: ",percentVar[1],"% variance"))+
    ylab(paste0("PC2: ",percentVar[2],"% variance"))+
    labs(title ="Distribution of samples \n rld Deseq2 inner method")+scale_color_manual(values=c(ppCor[1:2]))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                     panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                     axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                     axis.text.y = element_text(size=15),legend.position = "right")
  ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_rld_PCA_in_AMA_Young_RNA_library.pdf"),plot_rld_PCA,width = 10, height =10)
  
  
  # perform PCA analysis
  rlog_Mat <- assay(rld)
  PCA_data1 <- t(rlog_Mat)
  set.seed(19921010)
  md <- prep(PCA_data1, scale="none", center=TRUE)
  resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
  summary(resPPCA) 
  #Importance of component(s):
  #                 PC1    PC2    PC3    PC4     PC5
  #R2                 0.3498 0.1651 0.1240 0.09127 0.05611
  #Cumulative R2 0.3498 0.5149 0.6389 0.73017 0.78628
  
  variation<-data.frame(summary(resPPCA))
  
  plot(sDev(resPPCA))
  df <-as.data.frame(scores(resPPCA))
  rld_PCA<-merge(colData,df,by=0)
  rld_PCA$group<-factor(rld_PCA$group,levels=c("AMA","Young"))
  #rld_PCA$Gender_RRBS<-factor(rld_PCA$Gender_RRBS,levels=c("Female","Male"))
  plot_rld_PCA2<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour =group,shape= group))+
    geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
    xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
    ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
    labs(title ="Distribution of samples \n rld")+scale_color_manual(values=c(ppCor[1:4]))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                     panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                     axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                     axis.text.y = element_text(size=15),legend.position = "right")
  ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_rld_PCA_self_cal_in_AMA_Young_RNA_library.pdf"),plot_rld_PCA2,width = 10, height =10)
  
  #plot3d <- with(rld_PCA,scatterplot3d(PC1, PC2, PC3, color =c(rep(ppCor[2],5), rep(ppCor[1],5)), pch = c(rep(8,5), rep(16,5)),cex.symbols = 1.2, font.lab = 2, font.axis = 2))
  #legend(plot3d$xyz.convert(-15,0,10),pch=16,legend=levels(df$group), col =ppCor[1:2],cex = rel(1.1), bty = 'n',yjust=0.5, xjust = 0, horiz = F)
  
  # perform PCA analysis for count not normalization
  count_Mat <- log2(counts(dds2, normalized=F)+1)
  PCA_data1 <- t(count_Mat)
  set.seed(19921010)
  md <- prep(PCA_data1, scale="none", center=TRUE)
  resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
  summary(resPPCA) 
  #Importance of component(s):
  #                 PC1    PC2    PC3     PC4     PC5
  #R2            0.3193 0.1609 0.09265 0.07602 0.0565
  #Cumulative R2 0.3193 0.4802 0.57285 0.64887 0.7054
  
  variation<-data.frame(summary(resPPCA))
  #variation$PC1[1]
  #variation$PC2[1]
  plot(sDev(resPPCA))
  df <-as.data.frame(scores(resPPCA))
  count_PCA<-merge(colData,df,by=0)
  count_PCA$group<-factor(count_PCA$group,levels=c("AMA","Young"))
  #count_PCA$Gender_RRBS<-factor(count_PCA$Gender_RRBS,levels=c("Female","Male"))
  plot_not_normal_PCA<-ggplot(data=count_PCA, mapping=aes(x=PC1,y=PC2,colour = group,shape= group))+
    geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
    xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
    ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
    labs(title ="Distribution of samples \n count not normalization")+scale_color_manual(values=c(ppCor[1:4]))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                     panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                     axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                     axis.text.y = element_text(size=15),legend.position = "right")
  ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/All_filter_",tag_name,"_not_normal_PCA_in_AMA_Young_RNA_library.pdf"),plot_not_normal_PCA,width = 10, height =10)
  
  # perform PCA analysis for count  normalization
  count_Mat <- log2(counts(dds2, normalized=T)+1)
  PCA_data1 <- t(count_Mat)
  set.seed(19921010)
  md <- prep(PCA_data1, scale="none", center=TRUE)
  resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
  summary(resPPCA) 
  #Importance of component(s):
  #                 PC1    PC2    PC3    PC4     PC5
  #R2           0.2401 0.1557 0.1019 0.08844 0.06583
  # Cumulative R2 0.2401 0.3958 0.4978 0.58623 0.65206
  
  variation<-data.frame(summary(resPPCA))
  #variation$PC1[1]
  #variation$PC2[1]
  plot(sDev(resPPCA))
  df <-as.data.frame(scores(resPPCA))
  count_PCA<-merge(colData,df,by=0)
  count_PCA$group<-factor(count_PCA$group,levels=c("AMA","Young"))
  #count_PCA$Gender_RRBS<-factor(count_PCA$Gender_RRBS,levels=c("Female","Male"))
  plot_normal_PCA<-ggplot(data=count_PCA, mapping=aes(x=PC1,y=PC2,colour = group,shape= group))+
    geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
    xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
    ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
    labs(title ="Distribution of samples \n count normalization")+scale_color_manual(values=c(ppCor[1:4]))+
    theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                     panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                     axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                     axis.text.y = element_text(size=15),legend.position = "right")
  ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_normal_PCA_in_AMA_Young_RNA_library.pdf"),plot_normal_PCA,width = 10, height =10)
}
