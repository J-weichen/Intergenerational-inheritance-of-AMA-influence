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
library(dendextend)
library(pheatmap)

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
files_MCDAN_Large<- list.files("/mnt/data/chenwei/huangnana/2.mapping_result/bulk_RNA_count/4.count/feactureCounts/MCDAN_Large")
files_MCDAN_Small<- list.files("/mnt/data/chenwei/huangnana/2.mapping_result/bulk_RNA_count/4.count/feactureCounts/MCDAN_Small")
files_sIUGR_Large<- list.files("/mnt/data/chenwei/huangnana/2.mapping_result/bulk_RNA_count/4.count/feactureCounts/sIUGR_Large")
files_sIUGR_Small<- list.files("/mnt/data/chenwei/huangnana/2.mapping_result/bulk_RNA_count/4.count/feactureCounts/sIUGR_Small")

files_MCDAN_Large2<-paste0("/mnt/data/chenwei/huangnana/2.mapping_result/bulk_RNA_count/4.count/feactureCounts/MCDAN_Large/",files_MCDAN_Large)
files_MCDAN_Small2<- paste0("/mnt/data/chenwei/huangnana/2.mapping_result/bulk_RNA_count/4.count/feactureCounts/MCDAN_Small/",files_MCDAN_Small)
files_sIUGR_Large2<- paste0("/mnt/data/chenwei/huangnana/2.mapping_result/bulk_RNA_count/4.count/feactureCounts/sIUGR_Large/",files_sIUGR_Large)
files_sIUGR_Small2<- paste0("/mnt/data/chenwei/huangnana/2.mapping_result/bulk_RNA_count/4.count/feactureCounts/sIUGR_Small/",files_sIUGR_Small)

sample_name<-c(files_MCDAN_Large,files_MCDAN_Small,files_sIUGR_Large,files_sIUGR_Small)
files<-c(files_MCDAN_Large2,files_MCDAN_Small2,files_sIUGR_Large2,files_sIUGR_Small2)

read.delim(paste0("/mnt/data/chenwei/huangnana/2.mapping_result/bulk_RNA_count/4.count/feactureCounts/MCDAN_Large/",files_MCDAN_Large[1]),header = F, nrow=2,row.names = 1,sep = "\t")
sample1 <- read.delim(files[1],header = F,row.names = 1,sep = "\t",col.names = c("gene_id",substring(sample_name[1], 1, nchar(files_MCDAN_Large[1])-19)))
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
str(count.table);head(count.table);dim(count.table)#23077    24
write.table(count.table,file="/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file1_all_no_filter_Verification_MCDAN_sIUGR_rawcount.xls",sep="\t", quote=F, row.names=T,col.names=T)

#step one:read data
count.table2<-read.table(file="/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file1_all_no_filter_Verification_MCDAN_sIUGR_rawcount.xls",sep="\t", row.names=1,header =T)
head(count.table2)

#担心性别影响因此去除性染染色体上基因
data_chromX<-read.table("/mnt/data/chenwei/huangnana/0.script/bulk_RNA/ChromX_genelist.txt",head = F,sep = "\t"	)
data_chromY<-read.table("/mnt/data/chenwei/huangnana/0.script/bulk_RNA/ChromY_genelist.txt",head = F,sep = "\t")
head(data_chromX)
geneX<-unlist(as.character(unique(data_chromX$V13)))
geneY<-unlist(as.character(unique(data_chromY$V13)))
geneXY<-c(geneX,geneY);length(geneXY)#1246
count.table_XY_removed<-count.table2[-which(rownames(count.table2) %in% geneXY),]
head(count.table_XY_removed)
dim(count.table_XY_removed);dim(count.table2)#20230 20985
write.table(count.table_XY_removed,file="/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file2_all_no_filter_MCDAN_sIUGR_rawcount_XY_removed.xls",sep="\t", quote=F, row.names=T,col.names=T)
#count.table_XY_removed<-read.table(file="/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file2_all_no_filter_MCDAN_sIUGR_rawcount_XY_removed.xls",sep="\t", row.names=1,header =T)

#colData可以直接提供也可以自己构建
samplenames2<-colnames(count.table2)#将要分析的亚集列名字赋予此
group<-ifelse(grepl("MCDAN",samplenames2),"MCDAN","sIUGR")
group2<-ifelse(grepl("MCDAN",samplenames2),"MCDAN",ifelse(grepl("sIUGR_L",samplenames2),"sIUGR_L","sIUGR_S"))
group3<-ifelse(grepl("MCDAN_L",samplenames2),"MCDAN_L",ifelse(grepl("MCDAN_S",samplenames2),"MCDAN_S",ifelse(grepl("sIUGR_L",samplenames2),"sIUGR_L","sIUGR_S")))
group4<-ifelse(grepl("_L",samplenames2),"Large","Small")

colData_whole <- data.frame(sample = samplenames2, group=group,group2=group2,group3=group3,group4=group4)
rownames(colData_whole)<-colData_whole$sample

colData_whole$group<-factor(colData_whole$group,levels = c("MCDAN","sIUGR"))
colData_whole$group2<-factor(colData_whole$group2,levels = c("MCDAN","sIUGR_L","sIUGR_S"))
colData_whole$group3<-factor(colData_whole$group3,levels = c("MCDAN_L","MCDAN_S","sIUGR_L","sIUGR_S"))
colData_whole$group4<-factor(colData_whole$group4,levels = c("Large","Small"))

#colData_whole$type<-ifelse(colData_whole$sample %in% c("A2Q","A6Q","Y1Q","A7Q"),"Boy",ifelse(grepl("M",colData_whole$sample),"Mom","Girl"))
head(colData_whole)
write.table(colData_whole,file="/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file3_all_no_filter_MCDAN_sIUGR_analysis_metadata.csv",sep=",", quote=F, row.names=T,col.names=T)

#count.table_XY_removed<-read.table(file="/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file2_all_no_filter_MCDAN_sIUGR_rawcount_XY_removed.xls",sep="\t", row.names=1,header =T)
#colData_whole<-read.csv(file="/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file3_all_no_filter_MCDAN_sIUGR_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)
sample_notuse<-c("sIUGR_S1","sIUGR_S5","sIUGR_L5","MCDAN_S2")

#绘图
tag_name<-"All_sample_remove_notused" #test line
colData<-colData_whole[!(as.character(colData_whole$sample) %in% sample_notuse),]
expression_matrix<-count.table_XY_removed[,as.character(colData$sample)]

#step two : evaluation gene detective number
evaluation.table<-expression_matrix
evaluation.table[evaluation.table !=0]<- 1
gene_detected_number<-colSums(evaluation.table)
gene_detected_number<-data.frame(gene_detected_number)
colnames(gene_detected_number) <- "gene_number"
gene_detected_number2<-merge(gene_detected_number,colData,by=0)

gene_detected_number2$group3<-factor(gene_detected_number2$group3,levels = c("MCDAN_L","MCDAN_S","sIUGR_L","sIUGR_S"))
gene_detected_number2$Row.names<-as.character(gene_detected_number2$Row.names)
gene_detected_number2<-gene_detected_number2[order(gene_detected_number2$group3,decreasing = F),]
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

plot_gene_detected_number12<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number,fill=group3))+geom_bar(stat="identity",width=0.8)
plot_gene_detected_number4<-plot_gene_detected_number12+scale_fill_manual(values=ppCor[c(10:7)])+
  geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
  theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_gene_detected_number4

plot_gene_detected_number13<-ggplot(data=gene_detected_number2, mapping=aes(x=Row.names,y=gene_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_gene_detected_number5 <-plot_gene_detected_number13+scale_fill_manual(values=ppCor)+
  geom_text(stat="identity",aes(label=gene_number), color="black", size=2,position=position_stack(1.05))+
  theme_classic()+labs(x="",y="Number of genes",title="Number of gene detected")+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_gene_detected_number5

count_plot<-ggarrange(plot_gene_detected_number2, plot_gene_detected_number3,plot_gene_detected_number4,plot_gene_detected_number5,
                      labels = c("A", "B", "C", "D"),ncol = 1, nrow = 4)
ggsave(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_gene_detected_in_MCDAN_sIUGR_RNA_library.pdf"),count_plot,width = 200, height =400, units = "mm")
#}
#evaluation global effect
dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix, colData = colData, design = ~ group3)
#data pre_cleaning
dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
dds2 <- estimateSizeFactors(dds1)
dds3 <- DESeq(dds2)
#data transform:sample number < 30 -> rlog
#show transform result
rld <- rlog(dds2, blind = FALSE);vsd <- vst(dds2, blind = FALSE)
colnames(rld)
df <- bind_rows(
  as_data_frame(log2(counts(dds2, normalized=TRUE)[, 2:3]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 2:3]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 2:3]) %>% mutate(transformation = "vst"))
colnames(df)[1:2] <- c("sample2","sample3")  #[1] "MCDAN_L2"       "MCDAN_L3"       "transformation"
trans_plot<-ggplot(df, aes(x = sample2, y = sample3)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)
trans_plot
ggsave(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_trans_three_type_in_MCDAN_sIUGR_RNA_library.pdf"),trans_plot,width = 620, height =200, units = "mm")

#样本距离： RNA-Seq分析第一步通常是评估样本间的总体相似度
anno1 <- as.data.frame(colData(dds2)[, c("group","group3")])
#count not normalization
sampleDists <- dist(t(log2(counts(dds2,normalized=F)+1)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
p1<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "Log2_trans_no_normalization")

pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_distance_count_not_normal_in_MCDAN_sIUGR_RNA_library.pdf"),width = 10, height =10)
print(p1)
dev.off()
#count  normalization
sampleDists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
p2<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "Log2_trans_normal")

pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_distance_count_normal_in_MCDAN_sIUGR_RNA_library.pdf"),width = 10, height =10)
print(p2)
dev.off()
#count rld transform
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
p3<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "rld_trans")

pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_distance_count_rld_in_MCDAN_sIUGR_RNA_library.pdf"),width = 10, height =10)
print(p3)
dev.off()
#count vsd transform
sampleDists <- dist(t(assay(vsd)))
sampleDists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
p4<-pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "vsd_trans")

pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_distance_count_vsd_in_MCDAN_sIUGR_RNA_library.pdf"),width = 10, height =10)
print(p4)
dev.off()

#cluster
dists <- dist(t(log2(counts(dds2, normalized=F)+1)),method = "euclidean") 
hc <- hcluster(dists, method="pearson")
row_dend <- as.dendrogram(hc)
hc$order
table(colData$group3)

lable_cor<-rep(ppCor[2],20)
lable_cor[which(hc$order %in% c(1:6))] <- ppCor[1]
lable_cor[which(hc$order %in% c(7:11))] <- ppCor[3]
lable_cor[which(hc$order %in% c(17:20))] <- ppCor[6]
labels_colors(row_dend) <- lable_cor
pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_cluster_not_normal_in_MCDAN_sIUGR_RNA_library.pdf"),width = 10, height =6)
plot(row_dend, main = "samples similarity \n count not normalization (pearson)")
dev.off()
#count normalization
dists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)),method = "euclidean") 
hc <- hcluster(dists, method="pearson")
row_dend <- as.dendrogram(hc)
lable_cor<-rep(ppCor[2],20)
lable_cor[which(hc$order %in% c(1:6))] <- ppCor[1]
lable_cor[which(hc$order %in% c(7:11))] <- ppCor[3]
lable_cor[which(hc$order %in% c(17:20))] <- ppCor[6]
labels_colors(row_dend) <- lable_cor
pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_cluster_normal_in_MCDAN_sIUGR_RNA_library.pdf"),width = 10, height =6)
plot(row_dend, main = "samples similarity \n count normalization (pearson)")
dev.off()
#rld
dists <- dist(t(assay(rld)),method = "euclidean") 
hc <- hcluster(dists, method="pearson")
row_dend <- as.dendrogram(hc)
lable_cor<-rep(ppCor[2],20)
lable_cor[which(hc$order %in% c(1:6))] <- ppCor[1]
lable_cor[which(hc$order %in% c(7:11))] <- ppCor[3]
lable_cor[which(hc$order %in% c(17:20))] <- ppCor[6]
labels_colors(row_dend) <- lable_cor
pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_cluster_rld_in_MCDAN_sIUGR_RNA_library.pdf"),width = 10, height =6)
plot(row_dend, main = "samples similarity \n rlg (pearson)")
dev.off()
#vsd
dists <- dist(t(assay(vsd)),method = "euclidean") 
hc <- hcluster(dists, method="pearson")
row_dend <- as.dendrogram(hc)
lable_cor<-rep(ppCor[2],20)
lable_cor[which(hc$order %in% c(1:6))] <- ppCor[1]
lable_cor[which(hc$order %in% c(7:11))] <- ppCor[3]
lable_cor[which(hc$order %in% c(17:20))] <- ppCor[6]
labels_colors(row_dend) <- lable_cor
pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_cluster_vsd_in_MCDAN_sIUGR_RNA_library.pdf"),width = 10, height =6)
plot(row_dend, main = "samples similarity \n vsd (pearson)")
dev.off()

#另一个可视化样本-样本距离的方法就是主成分分析。
#DESeq2提供了专门的方法用于作图
plotPCA(rld,intgroup=c("group3"))
#plotPCA(rld,intgroup=c("Gender_RRBS"))

rld_PCA<-plotPCA(rld,intgroup=c("group","group3"),returnData = TRUE)
percentVar<-round(100*attr(rld_PCA,"percentVar"),1)
rld_PCA$group3<-factor(rld_PCA$group3,levels= c("MCDAN_L","MCDAN_S","sIUGR_L","sIUGR_S"))

plot_rld_PCA<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour = group3,shape= group3))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  labs(title ="Distribution of samples \n rld Deseq2 inner method")+scale_color_manual(values=c(ppCor))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
ggsave(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_rld_PCA_in_MCDAN_sIUGR_RNA_library.pdf"),plot_rld_PCA,width = 10, height =10)


# perform PCA analysis
rlog_Mat <- assay(rld)
PCA_data1 <- t(rlog_Mat)
set.seed(19921010)
md <- prep(PCA_data1, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
summary(resPPCA) 
#Importance of component(s):
#                 PC1    PC2    PC3    PC4     PC5
#R2             0.2516 0.1538 0.1269 0.07308 0.06482
#Cumulative R2 0.2516 0.4054 0.5323 0.60540 0.67022

variation<-data.frame(summary(resPPCA))

plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
rld_PCA<-merge(colData,df,by=0)
rld_PCA$group3<-factor(rld_PCA$group3,levels=c("MCDAN_L","MCDAN_S","sIUGR_L","sIUGR_S"))

#rld_PCA$Gender_RRBS<-factor(rld_PCA$Gender_RRBS,levels=c("Female","Male"))
plot_rld_PCA2<-ggplot(data=rld_PCA, mapping=aes(x=PC1,y=PC2,colour =group3,shape= group3))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples \n rld")+scale_color_manual(values=c(ppCor[1:4]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_rld_PCA2
ggsave(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_rld_PCA_self_cal_in_MCDAN_sIUGR_RNA_library.pdf"),plot_rld_PCA2,width = 10, height =10)

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
#R2              0.3186 0.1237 0.07158 0.05229 0.04833
#Cumulative R2 0.3186 0.4423 0.51390 0.56619 0.61452
 

variation<-data.frame(summary(resPPCA))
plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
count_PCA<-merge(colData,df,by=0)
count_PCA$group3<-factor(count_PCA$group3,levels= c("MCDAN_L","MCDAN_S","sIUGR_L","sIUGR_S"))
#count_PCA$Gender_RRBS<-factor(count_PCA$Gender_RRBS,levels=c("Female","Male"))
plot_not_normal_PCA<-ggplot(data=count_PCA, mapping=aes(x=PC1,y=PC2,colour = group3,shape= group3))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples \n count not normalization")+scale_color_manual(values=c(ppCor[1:4]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_not_normal_PCA
ggsave(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/All_no_filter_",tag_name,"_not_normal_PCA_in_MCDAN_sIUGR_RNA_library.pdf"),plot_not_normal_PCA,width = 10, height =10)

# perform PCA analysis for count  normalization
count_Mat <- log2(counts(dds2, normalized=T)+1)
PCA_data1 <- t(count_Mat)
set.seed(19921010)
md <- prep(PCA_data1, scale="none", center=TRUE)
resPPCA <- pcaMethods::pca(md, method="ppca", center=FALSE,nPcs = 5)
summary(resPPCA) 
#Importance of component(s):
#                 PC1    PC2    PC3    PC4     PC5
#R2            0.1847 0.1309 0.09017 0.06489 0.05903
#Cumulative R2 0.1847 0.3156 0.40575 0.47064 0.52967

variation<-data.frame(summary(resPPCA))
#variation$PC1[1]
#variation$PC2[1]
plot(sDev(resPPCA))
df <-as.data.frame(scores(resPPCA))
count_PCA<-merge(colData,df,by=0)
count_PCA$group3<-factor(count_PCA$group3,levels= c("MCDAN_L","MCDAN_S","sIUGR_L","sIUGR_S"))
#count_PCA$Gender_RRBS<-factor(count_PCA$Gender_RRBS,levels=c("Female","Male"))
plot_normal_PCA<-ggplot(data=count_PCA, mapping=aes(x=PC1,y=PC2,colour = group3,shape= group3))+
  geom_point(stat= "identity",size=4,alpha=0.5,show.legend = TRUE)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples \n count normalization")+scale_color_manual(values=c(ppCor[1:4]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
plot_normal_PCA
ggsave(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/",tag_name,"_normal_PCA_in_MCDAN_sIUGR_RNA_library.pdf"),plot_normal_PCA,width = 10, height =10)

#3D and lable names
pdf(paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/",tag_name,"_normal_PCA_in_MCDAN_sIUGR_RNA_library_PCAscreeplot_2D_3D.pdf"))
plot3d <- with(count_PCA, scatterplot3d(PC1, PC2, PC3,id=T,color=ppCor[c(rep(1,6),rep(2,5),rep(3,5),rep(5,4))],
                                        pch = c(rep(16,20)),
                                        cex.symbols = 1.2, font.lab = 2, font.axis = 2))
legend("right",pch=16, legend = levels(count_PCA$group3),col =ppCor[c(1:3,5)],horiz =F)
#legend("right",pch=16,legend=levels(count_PCA$group3), col =ppCor[c(1:3,5)],cex = rel(1.1), bty = 'n',yjust=0, xjust = 0.5, horiz = F)

PCA_plot<-ggplot(data=count_PCA, mapping=aes(x=PC1,y=PC2,colour = group3))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=2)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC2: ",round(variation$PC2[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(1:3,5)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot
PCA_plot2<-ggplot(data=count_PCA, mapping=aes(x=PC1,y=PC3,colour = group3))+
  geom_point(stat= "identity",size=4,alpha=0.7,show.legend = TRUE)+
  geom_text(stat="identity",aes(label=Row.names), color="black", size=2)+
  xlab(paste0("PC1: ",round(variation$PC1[1]*100,1),"% variance"))+
  ylab(paste0("PC3: ",round(variation$PC3[1]*100,1),"% variance"))+
  labs(title ="Distribution of samples")+scale_color_manual(values=c(ppCor[c(1:3,5)]))+
  theme_bw()+theme(panel.border = element_rect(colour="grey",fill=NA),
                   panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
                   axis.title = element_text(size=15),axis.text.x = element_text(size=15), 
                   axis.text.y = element_text(size=15),legend.position = "right")
PCA_plot2
dev.off()

#step seven: count matrix generation
##the normalized counts (normalized for library size, which is the total number of gene counts per sample, while accounting for library composition)
normalized_counts <- counts(dds3, normalized=TRUE)
#根据基因在不同的样本中表达变化的差异程度mad值对数据排序，差异越大的基因排位越前。
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.table(normalized_counts,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file4_",tag_name,"_MCDAN_sIUGR_deseq2_normalized_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
#log2
normalized_counts_log2 <- log2(counts(dds3, normalized=T)+1)
normalized_counts_mad2 <- apply(normalized_counts_log2, 1, mad)
normalized_counts_log2 <- normalized_counts_log2[order(normalized_counts_mad2, decreasing=T), ]
write.table(normalized_counts_log2,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file4_",tag_name,"_MCDAN_sIUGR_deseq2_deseq2_log2_normalized_count_order.txt"),sep="\t", quote=F, row.names=T,col.names=T)
# other log transform
rlog_dds3_Mat <- assay(rld)
rlog_dds3_Mat <- rlog_dds3_Mat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlog_dds3_Mat,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file4_",tag_name,"_MCDAN_sIUGR_deseq2_deseq2_rlg_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
# 只在Linux下能运行，目的是填补表格左上角内容
#system(paste("sed -i '1 s/^/ID\t/'", "ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.xls"))

#call DEGs
#step2 reading count and colData
count.table_XY_removed<-read.table(file="/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file2_all_no_filter_MCDAN_sIUGR_rawcount_XY_removed.xls",sep="\t", row.names=1,header =T)
colData_whole<-read.csv(file="/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file3_all_no_filter_MCDAN_sIUGR_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)
sample_notuse<-c("sIUGR_S1","sIUGR_S5","sIUGR_L5","MCDAN_S2")
colData_whole<-colData_whole[!(as.character(colData_whole$sample) %in% sample_notuse),]
head(colData_whole)
colData_whole$group3


#step3 selected specific compare_group: Disease compared with CTRL 
type<-c("sIUGR_L","sIUGR_S")
for (tag_name in type){
  # tag_name<-"sIUGR_S" #test line
  print(tag_name)
  sample_type<-c("MCDAN_L","MCDAN_S",tag_name)
  colData<-colData_whole[which(as.character(colData_whole$group3) %in% sample_type),]
  expression_matrix<-count.table_XY_removed[,as.character(colData$sample)]
  dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix, colData = colData, design = ~ group)
  #step 4 deseq2 normalization and DEGs called
  dds1 <- dds0[rowSums(counts(dds0)) > 1, ]
  dds2 <- estimateSizeFactors(dds1)
  dds3 <- DESeq(dds2)
  rld <- rlog(dds2, blind = FALSE);vsd <- vst(dds2, blind = FALSE)
  
  ##the normalized counts (normalized for library size, which is the total number of gene counts per sample, while accounting for library composition)
  normalized_counts <- counts(dds3, normalized=TRUE)
  normalized_counts_mad <- apply(normalized_counts, 1, mad)
  normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
  write.table(normalized_counts,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file5_",tag_name,"_MCDAN_sIUGR_deseq2_normalized_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
  #log2
  normalized_counts_log2 <- log2(counts(dds3, normalized=T)+1)
  normalized_counts_mad2 <- apply(normalized_counts_log2, 1, mad)
  normalized_counts_log2 <- normalized_counts_log2[order(normalized_counts_mad2, decreasing=T), ]
  write.table(normalized_counts_log2,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file5_",tag_name,"_MCDAN_sIUGR_deseq2_deseq2_log2_normalized_count_order.txt"),sep="\t", quote=F, row.names=T,col.names=T)
  # other log transform
  rlog_dds3_Mat <- assay(rld)
  rlog_dds3_Mat <- rlog_dds3_Mat[order(normalized_counts_mad, decreasing=T), ]
  write.table(rlog_dds3_Mat,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file5_",tag_name,"_MCDAN_sIUGR_deseq2_deseq2_rlg_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
  
  #step eight: DEGs output
  # 定义变量：如果想比较不同的组，只需在这修改即可？能否直接写一个两两比较的循环进行自动处理？
  Conditions<-"group";sampleA <-"MCDAN";sampleB <-"sIUGR" 
  contrastV <- c(Conditions, sampleA, sampleB)
  res1 <- results(dds3, contrast=contrastV) 
  summary(res1)
  
  ##给DESeq2的原始输出结果增加样品平均表达信息，使得结果更容易理解和解析。
  # 获得第一组数据均值
  baseA <- counts(dds3, normalized=TRUE)[, colData(dds3)$group == sampleA]
  if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
  } else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
  }
  colnames(baseMeanA) <-paste(Conditions,"-",sampleA,sep="")
  #head(baseMeanA) 
  
  # 获得第二组数据均值
  baseB <- counts(dds3, normalized=TRUE)[, colData(dds3)$group == sampleB]
  if (is.vector(baseB)){
    baseMeanB <- as.data.frame(baseB)
  } else {
    baseMeanB <- as.data.frame(rowMeans(baseB))
  }
  colnames(baseMeanB) <- paste(Conditions,"-",sampleB,sep="")
  #head(baseMeanB) 
  
  # 结果组合
  res2 <- cbind(baseMeanA, baseMeanB, as.data.frame(res1));head(res2) 
  res3 <- cbind(ID=rownames(res2), as.data.frame(res2));head(res3)
  # 校正后p-value为NA的复制为1
  res3$padj[is.na(res3$padj)] <- 1;res3$pvalue[is.na(res3$pvalue)] <- 1
  # 按pvalue排序, 把差异大的基因放前面
  res3 <- res3[order(res3$pvalue,decreasing = F),]
  head(res3);tail(res3) ;dim(res3)
  #whole RNA_data combined
  head(as.data.frame(res3));dim(as.data.frame(res3))
  head(normalized_counts_log2);dim(normalized_counts_log2)
  
  merge_data <- merge(normalized_counts_log2,as.data.frame(res3),by="row.names",sort=FALSE)
  merge_data <- merge_data[order(merge_data$pvalue,decreasing = F),]
  file <- paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file6_",tag_name,"_special_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
  write.table(merge_data, file=file, sep="\t", quote=F, row.names=F,col.names=T)
  
  table(merge_data$pvalue<0.05)  
  table(merge_data$pvalue<0.01) 
  table(merge_data$pvalue<0.05 & abs(merge_data$log2FoldChange)>=0.585)   
  table(merge_data$padj<0.1)  
  table(merge_data$padj<0.05)   
  table(merge_data$padj<0.05 & abs(merge_data$log2FoldChange)>=0.585)
  
  #差异基因的进一步筛选
  # p_value<0.05 & Foldchange >=1.5
  res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
  res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
  nrow(res_de_up);nrow(res_de_dw)
  DEG_information<-rbind(res_de_up,res_de_dw)
  file <- paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file7_",tag_name,"_special_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
  write.table(as.data.frame(DEG_information), file=file, sep="\t", quote=F, row.names=T, col.names=T)
  
  #plot 
  head(merge_data) 
  merge_data$threshold = as.factor(ifelse(merge_data$pvalue < 0.05 & abs(merge_data$log2FoldChange)>= 0.585, ifelse(merge_data$log2FoldChange>= 0.585 ,'Up','Down'),'NoSignifi'))
  table(merge_data$threshold)
  merge_data$threshold<-factor(merge_data$threshold,level=c("Up","Down","NoSignifi"))
  range(merge_data$log2FoldChange)#-5.598469  4.432278
  range(-log10(merge_data$pvalue))#0.00000 20.19249
  
  plot_vocano<-ggplot(data = merge_data, aes(x = log2FoldChange, y = -log10(pvalue), colour=threshold)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("red","blue","grey")) +
    scale_x_continuous(limits=c(-ceiling(max(abs(merge_data$log2FoldChange))),ceiling(max(abs(merge_data$log2FoldChange)))),
                       breaks=seq(-ceiling(max(abs(merge_data$log2FoldChange))),ceiling(max(abs(merge_data$log2FoldChange))), 1))+
    #xlim(-ceiling(max(abs(merge_data$log2FoldChange)))-1,ceiling(max(abs(merge_data$log2FoldChange)))+1)+
    ylim(0,ceiling(max(-log10(merge_data$pvalue))))+
    geom_vline(xintercept=c(-0.585,0.585),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(Foldchange(AMA/Young)) in gene expression level",y="-log10 (p value)",
         title=paste(tag_name,":",sampleA,"_vs_",sampleB,":Up_DEGs: ",length(which(merge_data$threshold == "Up"))," & Down_DEG: ",length(which(merge_data$threshold == "Down"))," [p value<0.05 & FoldChange > 1.5]",sep="")) +
    theme_bw()+
    theme(panel.border = element_rect(colour="black",fill=NA),panel.grid.major =element_blank())+
    theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5), legend.position="right", legend.title = element_blank(),axis.line = element_line(colour="black"))+
    theme(axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
  plot_vocano
  ggsave(plot_vocano,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/",tag_name,"_",sampleA,"_vs_",sampleB,"_vocano_DEGs.pdf"), width = 200, height = 200, units = "mm")
  
  
  #绘制各组特异高表达基因热图
  rownames(merge_data)<- merge_data$Row.names
  merge_data<-merge_data[,-1]
  res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585 )
  res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585 )
  de_id_whole = rbind(res_de_up, res_de_dw)
  nrow(de_id_whole);head(de_id_whole)
  colnames(de_id_whole)
 # mat<-de_id_whole[,1:16]
  mat<-de_id_whole[,1:15]
  anno1 <- data.frame(sample=colnames(mat))
  anno1$group<-ifelse(grepl("MCDAN",anno1$sample),"MCDAN","sIUGR")
  rownames(anno1)<-anno1$sample
  ann_colors = list(group=c(MCDAN=ppCor[2],sIUGR=ppCor[1]))
  range(mat)# 0.00000 18.32706
  pheatmap(mat, cluster_rows=F, cluster_cols =F,scale="row", annotation_col=anno1[2],show_rownames = F,main="scale_Genes_pvalue<0.05 & FC >1.5")
  pheatmap(mat,scale =  "row",
           color = colorRampPalette(colors = c("blue","white","red"))(100),
           annotation_col = anno1[2],main="Genes_pvalue <0.05&log2FoldChange>=0.585")
  
  bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01)) 
  # 做热图： 
  heat_plot<-pheatmap(mat, scale = "row", 
                      # color = colorRampPalette(colors = c("blue","white","red"))(100),
                      color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), 
                      cluster_row =FALSE,cluster_col =FALSE,                    
                      legend_breaks=seq(-3,3,1),  breaks=bk,
                      annotation_colors = ann_colors, 
                      gaps_row = nrow(res_de_up),gaps_col =11,cutree_col = 2,
                      show_rownames = F,
                      annotation_col = anno1[2],
                      main=paste0(tag_name,":DEGs_pvalue <0.05& FoldChange>=1.5")
  )
  
  pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/heatmap_All_filter_",tag_name,"_",sampleA,"_vs_",sampleB,"_DEGs_p005_FC1.5.pdf"),width = 7, height =7)
  print(heat_plot)
  dev.off()
}

#step3 selected specific compare_group: Small compared with Large 
type<-c("MCDAN","sIUGR")
for (tag_name in type){
  # tag_name<-"sIUGR" #test line
  print(tag_name)
  sample_type<-c(paste0(tag_name,"_L"),paste0(tag_name,"_S"))
  colData<-colData_whole[which(as.character(colData_whole$group3) %in% sample_type),]
  expression_matrix<-count.table_XY_removed[,as.character(colData$sample)]
  dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix, colData = colData, design = ~ group4)
  #step 4 deseq2 normalization and DEGs called
  dds1 <- dds0[rowSums(counts(dds0)) > 1, ]
  dds2 <- estimateSizeFactors(dds1)
  dds3 <- DESeq(dds2)
  rld <- rlog(dds2, blind = FALSE);vsd <- vst(dds2, blind = FALSE)
  
  ##the normalized counts (normalized for library size, which is the total number of gene counts per sample, while accounting for library composition)
  normalized_counts <- counts(dds3, normalized=TRUE)
  normalized_counts_mad <- apply(normalized_counts, 1, mad)
  normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
  write.table(normalized_counts,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file5_",tag_name,"_special_lager_small_deseq2_normalized_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
  #log2
  normalized_counts_log2 <- log2(counts(dds3, normalized=T)+1)
  normalized_counts_mad2 <- apply(normalized_counts_log2, 1, mad)
  normalized_counts_log2 <- normalized_counts_log2[order(normalized_counts_mad2, decreasing=T), ]
  write.table(normalized_counts_log2,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file5_",tag_name,"_special_lager_small_deseq2_deseq2_log2_normalized_count_order.txt"),sep="\t", quote=F, row.names=T,col.names=T)
  # other log transform
  rlog_dds3_Mat <- assay(rld)
  rlog_dds3_Mat <- rlog_dds3_Mat[order(normalized_counts_mad, decreasing=T), ]
  write.table(rlog_dds3_Mat,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file5_",tag_name,"_special_lager_small_deseq2_deseq2_rlg_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
  
  #step eight: DEGs output
  # 定义变量：如果想比较不同的组，只需在这修改即可？能否直接写一个两两比较的循环进行自动处理？
  Conditions<-"group4";sampleA <-"Large";sampleB <-"Small" 
  contrastV <- c(Conditions, sampleA, sampleB)
  res1 <- results(dds3, contrast=contrastV) 
  summary(res1)
  
  ##给DESeq2的原始输出结果增加样品平均表达信息，使得结果更容易理解和解析。
  # 获得第一组数据均值
  baseA <- counts(dds3, normalized=TRUE)[, colData(dds3)$group4 == sampleA]
  if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
  } else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
  }
  colnames(baseMeanA) <-paste(Conditions,"-",sampleA,sep="")
  head(baseMeanA) 
  
  # 获得第二组数据均值
  baseB <- counts(dds3, normalized=TRUE)[, colData(dds3)$group4 == sampleB]
  if (is.vector(baseB)){
    baseMeanB <- as.data.frame(baseB)
  } else {
    baseMeanB <- as.data.frame(rowMeans(baseB))
  }
  colnames(baseMeanB) <- paste(Conditions,"-",sampleB,sep="")
  head(baseMeanB) 
  
  # 结果组合
  res2 <- cbind(baseMeanA, baseMeanB, as.data.frame(res1));head(res2) 
  res3 <- cbind(ID=rownames(res2), as.data.frame(res2));head(res3)
  # 校正后p-value为NA的复制为1
  res3$padj[is.na(res3$padj)] <- 1;res3$pvalue[is.na(res3$pvalue)] <- 1
  # 按pvalue排序, 把差异大的基因放前面
  res3 <- res3[order(res3$pvalue,decreasing = F),]
  head(res3);tail(res3) ;dim(res3)
  #whole RNA_data combined
  head(as.data.frame(res3));dim(as.data.frame(res3))
  head(normalized_counts_log2);dim(normalized_counts_log2)
  
  merge_data <- merge(normalized_counts_log2,as.data.frame(res3),by="row.names",sort=FALSE)
  merge_data <- merge_data[order(merge_data$pvalue,decreasing = F),]
  file <- paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file6_",tag_name,"_special_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
  write.table(merge_data, file=file, sep="\t", quote=F, row.names=F,col.names=T)
  
  table(merge_data$pvalue<0.05)  
  table(merge_data$pvalue<0.01) 
  table(merge_data$pvalue<0.05 & abs(merge_data$log2FoldChange)>=0.585)   
  table(merge_data$pvalue<0.05 & abs(merge_data$log2FoldChange)>=1)   
  table(merge_data$padj<0.1)  
  table(merge_data$padj<0.05)   
  table(merge_data$padj<0.05 & abs(merge_data$log2FoldChange)>=0.585)
  
  #差异基因的进一步筛选
  # p_value<0.05 & Foldchange >=1.5
  res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
  res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
  nrow(res_de_up);nrow(res_de_dw)
  DEG_information<-rbind(res_de_up,res_de_dw)
  file <- paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file7_",tag_name,"_special_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
  write.table(as.data.frame(DEG_information), file=file, sep="\t", quote=F, row.names=T, col.names=T)
  
  #plot 
  head(merge_data) 
  merge_data$threshold = as.factor(ifelse(merge_data$pvalue < 0.05 & abs(merge_data$log2FoldChange)>= 0.585, ifelse(merge_data$log2FoldChange>= 0.585 ,'Up','Down'),'NoSignifi'))
  table(merge_data$threshold)
  merge_data$threshold<-factor(merge_data$threshold,level=c("Up","Down","NoSignifi"))
  range(merge_data$log2FoldChange)#-5.598469  4.432278
  range(-log10(merge_data$pvalue))#0.00000 20.19249
  
  plot_vocano<-ggplot(data = merge_data, aes(x = log2FoldChange, y = -log10(pvalue), colour=threshold)) +
    geom_point(alpha=0.4, size=3.5) +
    scale_color_manual(values=c("red","blue","grey")) +
    scale_x_continuous(limits=c(-ceiling(max(abs(merge_data$log2FoldChange))),ceiling(max(abs(merge_data$log2FoldChange)))),
                       breaks=seq(-ceiling(max(abs(merge_data$log2FoldChange))),ceiling(max(abs(merge_data$log2FoldChange))), 1))+
    #xlim(-ceiling(max(abs(merge_data$log2FoldChange)))-1,ceiling(max(abs(merge_data$log2FoldChange)))+1)+
    ylim(0,ceiling(max(-log10(merge_data$pvalue))))+
    geom_vline(xintercept=c(-0.585,0.585),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(Foldchange(AMA/Young)) in gene expression level",y="-log10 (p value)",
         title=paste(tag_name,":",sampleA,"_vs_",sampleB,":Up_DEGs: ",length(which(merge_data$threshold == "Up"))," & Down_DEG: ",length(which(merge_data$threshold == "Down"))," [p value<0.05 & FoldChange > 1.5]",sep="")) +
    theme_bw()+
    theme(panel.border = element_rect(colour="black",fill=NA),panel.grid.major =element_blank())+
    theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5), legend.position="right", legend.title = element_blank(),axis.line = element_line(colour="black"))+
    theme(axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
  plot_vocano
  ggsave(plot_vocano,file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/",tag_name,"_",sampleA,"_vs_",sampleB,"_vocano_DEGs.pdf"), width = 200, height = 200, units = "mm")
  
  
  #绘制各组特异高表达基因热图
  rownames(merge_data)<- merge_data$Row.names
  merge_data<-merge_data[,-1]
  res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585 )
  res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585 )
  de_id_whole = rbind(res_de_up, res_de_dw)
  nrow(de_id_whole);head(de_id_whole)
  colnames(de_id_whole)
  mat<-de_id_whole[,1:11]
 # mat<-de_id_whole[,1:9]
  anno1 <- data.frame(sample=colnames(mat))
  anno1$group<-ifelse(grepl("_L",anno1$sample),"Lager","Small")
  rownames(anno1)<-anno1$sample
  ann_colors = list(group=c(Lager=ppCor[2],Small=ppCor[1]))
  range(mat)# 0.00000 18.32706
  pheatmap(mat, cluster_rows=F, cluster_cols =F,scale="row", annotation_col=anno1[2],show_rownames = F,main="scale_Genes_pvalue<0.05 & FC >1.5")
  pheatmap(mat,scale =  "row",
           color = colorRampPalette(colors = c("blue","white","red"))(100),
           annotation_col = anno1[2],main="Genes_pvalue <0.05&log2FoldChange>=0.585")
  
  bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01)) 
  # 做热图： 
  heat_plot<-pheatmap(mat, scale = "row", 
                      # color = colorRampPalette(colors = c("blue","white","red"))(100),
                      color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), 
                      cluster_row =FALSE,cluster_col =FALSE,                    
                      legend_breaks=seq(-3,3,1),  breaks=bk,
                      annotation_colors = ann_colors, 
                      gaps_row = nrow(res_de_up),gaps_col =5,cutree_col = 2,
                      show_rownames = F,
                      annotation_col = anno1[2],
                      main=paste0(tag_name,":DEGs_pvalue <0.05& FoldChange>=1.5")
  )
  
  pdf(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/heatmap_All_filter_",tag_name,"_",sampleA,"_vs_",sampleB,"_DEGs_p005_FC1.5.pdf"),width = 7, height =7)
  print(heat_plot)
  dev.off()
}

#以下待修整
#绘制配对表达量变化
normalized_counts<-read.table(file=paste0("/mnt/data/chenwei/huangnana/4.bulk_RNA_result/file5_",tag_name,"_special_lager_small_deseq2_normalized_count.txt"),sep="\t", row.names=1,header =T)
head(normalized_counts);range(normalized_counts)
normalized_counts2<-normalized_counts+0.01
fold_change_data<-as.data.frame(t(apply(normalized_counts2,1,function(x) c(x[1]/x[7],x[2]/x[8],x[3]/x[9],x[4]/x[10],x[5]/x[11],x[6]/x[12]))))
fold_change_data$mad_value <- apply(fold_change_data, 1, mad)
fold_change_data$FC_median<-apply(fold_change_data[,1:6], 1, median)
fold_change_data$FC_sd<-apply(fold_change_data[,1:6], 1, sd)
head(fold_change_data);dim(fold_change_data)

fold_change_data <- fold_change_data[order(fold_change_data$mad_value, decreasing=F), ]
#fold_change_data[fold_change_data =="Inf"]<-200
range(fold_change_data$mad_value)
range(fold_change_data$FC_sd)
plot_data<-fold_change_data[which((fold_change_data$FC_median>2 |fold_change_data$FC_median<0.5) & fold_change_data$mad_value <20),]
dim(plot_data)# 2120    9

heat_plot3 <-pheatmap(log2(plot_data[,1:6]), cluster_row =T,cluster_col =FALSE,na_col = "black",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      #annotation_col = annotation_col,annotation_row=annotation_row,
                      #annotation_colors = ann_colors, 
                      # gaps_row =c(4,8,32),gaps_col =c(10),cutree_col = 2,
                      #treeheight_col = 20, 
                      #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      border=FALSE,scale ="none",
                      color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), 
                      # filename = "D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/DMR_methlevel_genes_over_DEGs_single_heatmap.pdf",
                      main ="Kids DNA methylation level of meth_trans_inter genes original value",angle_col ="90")
