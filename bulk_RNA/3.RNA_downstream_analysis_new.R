#转录组总的分析流程：https://www.jianshu.com/p/fb15249200d7
##另外可以参考：https://www.cnblogs.com/weipeng-loaded/p/10601285.html
#待学习Limma差异分析、火山图、功能富集： https://www.jianshu.com/p/0a8841715804
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("bladderbatch")
#参考：http://www.gene123.com/n/2885/
#数据读取
#如果已经有了整合好的readscount matrixs
rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
#setwd("D:/3.秦萌高龄甲基化/1.RNA_result/merge_count/feactureCounts")
library("dendextend")
library("pcaMethods")
library("scatterplot3d")
library("DESeq2") 
library("dplyr")
library("stringr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("limma")
library("amap")
#library("PoiClaClu")
library("gplots")
library("enrichplot")
library("ggpubr")
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)#org.Hs.eg.db 的数据类型以及使用简介：可用于数据类型转换
keytypes(org.Hs.eg.db)
grid.newpage()
#set colors
pal <- pal_npg("nrc", alpha=1)(9)
show_col(pal)
pal
#[1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF"
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
show_col(nejm)
nejm
#[1] "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF" "#EE4C97FF"
#decide color
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

#n=5;barplot(rep(1,n), col=colorRampPalette(colors = c('red', 'white'))( n ))
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
show_col(Cells_col_raw)
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
show_col(pal1)
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
show_col(pal2)
pal3<- pal_aaas("default",alpha=1)(10)
show_col(pal3)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
show_col(pal4)
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
show_col(pal5)
pal5
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor_all)
ppCor_all2<-ppCor_all[c(c(1:8),c(11:14),17,19,23,24,c(26:32),33,34,c(38:44))]
Cells_col<-colorRampPalette(colors = ppCor_all2[c(24,2:9,1,12,14:16,19:23,13)])(80)
show_col(Cells_col)

#step2 reading count and colData
count.table_XY_removed<-read.table(file="/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Verification_AMA_versus_Young_rawcount_XY_removed.xls",sep="\t", row.names=1,header =T)

colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/verification_AMA_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)

#step3 selected specific compare_group
type<-c("Mom","Kid")

for (tag_name in type){
  # tag_name<-"Mom" #test line
  colData<-colData_whole[which(colData_whole$generation==tag_name),]
  expression_matrix<-count.table_XY_removed[,as.character(colData$sample)]
  dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix, colData = colData, design = ~ group )
  #26485    10
  #22083    10
  #step 4 deseq2 normalization and DEGs called
  dds1 <- dds0[rowSums(counts(dds0)) > 1, ]
  dds2 <- estimateSizeFactors(dds1)
  dds3 <- DESeq(dds2)
  AMA_deseq2_log2_normalized_count<-log2(counts(dds3, normalized=TRUE)+1)
  write.table(AMA_deseq2_log2_normalized_count,file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_AMA_deseq2_log2_normalized_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
 #step four : show data transformation
  rld <- rlog(dds3, blind = FALSE);head(assay(rld), 3)
  vsd <- vst(dds3, blind = FALSE);head(assay(vsd), 3)
  colnames(rld)#"P22" "P24" "P25" "P38" "P31" "P32" "P33" "P34" "P36" "P41"
  df <- bind_rows(
  as_data_frame(log2(counts(dds3, normalized=TRUE)[, 2:3]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 2:3]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 2:3]) %>% mutate(transformation = "vst"))
  #colnames(df)[1:2] <- c("A3M","A4M")  
  ggplot(df) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)

  #标准化后进行样本评估：样品层级聚类分析，判断样品的相似性和组间组内差异
  nor_Mat<-counts(dds3, normalized=TRUE)
  No_nor_Mat<-counts(dds3, normalized=F)
  rlog_dds3_Mat <- assay(rld)
  vlog_dds3_Mat <- assay(vsd)
# 生成颜色
hmcol <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
# 计算相关性pearson correlation:rlog_dds3_Mat
hc <- hcluster(t(rlog_dds3_Mat), method="pearson")
plot(hc, cex = 2) 
#pearson_cor <- as.matrix(cor(nor_Mat, method="pearson"))
pearson_cor <- as.matrix(cor(rlog_dds3_Mat, method="pearson"))
#pearson_cor <- as.matrix(cor(vlog_dds3_Mat, method="pearson"))
heatmap.2(pearson_cor, 
          Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, 
          margins=c(11,11), 
          main="Pearson correlation(rlg)")

hc <- hcluster(t(rlog_dds3_Mat), method="spearman")
plot(hc, cex = 2) 
spearman_cor <- as.matrix(cor(rlog_dds3_Mat, method="spearman"))
heatmap.2(spearman_cor, 
          Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, 
          margins=c(11,11), 
          main="Spearman  correlation(rlg)")

# 计算相关性pearson correlation:nor_Mat
hc <- hcluster(t(nor_Mat), method="pearson")
plot(hc, cex = 2) 
pearson_cor <- as.matrix(cor(nor_Mat, method="pearson"))
heatmap.2(pearson_cor,Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11),  main="Pearson correlation \n (count normalization)")

hc <- hcluster(t(nor_Mat), method="spearman")
plot(hc, cex = 2) 
spearman_cor <- as.matrix(cor(nor_Mat, method="spearman"))
heatmap.2(spearman_cor,Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11),  main="Spearman  correlation \n (count normalization)")

# 计算相关性pearson correlation:No_nor_Mat
hc <- hcluster(t(No_nor_Mat), method="pearson")
plot(hc, cex = 2) 
pearson_cor <- as.matrix(cor(No_nor_Mat, method="pearson"))
heatmap.2(pearson_cor,Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11),  main="Pearson correlation \n (count no normalization)")

hc <- hcluster(t(No_nor_Mat), method="spearman")
plot(hc, cex = 2) 
spearman_cor <- as.matrix(cor(No_nor_Mat, method="spearman"))
heatmap.2(spearman_cor,Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11),  main="Spearman  correlation \n (count no normalization)")

#pdf("ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.pearson.pdf", pointsize=10)
#dev.off()

#step seven: count matrix generation
#获取标准化后的数据
##the normalized counts (normalized for library size, which is the total number of gene counts per sample, while accounting for library composition)
normalized_counts <- counts(dds3, normalized=TRUE)
#根据基因在不同的样本中表达变化的差异程度mad值对数据排序，差异越大的基因排位越前。
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.table(normalized_counts,file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_AMA_deseq2_normalized_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
#log2
normalized_counts_log2 <- log2(counts(dds3, normalized=T)+1)
normalized_counts_mad2 <- apply(normalized_counts_log2, 1, mad)
normalized_counts_log2 <- normalized_counts_log2[order(normalized_counts_mad2, decreasing=T), ]
write.table(normalized_counts_log2,file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_AMA_deseq2_log2_normalized_count_order.txt"),sep="\t", quote=F, row.names=T,col.names=T)
# other log transform
rlog_dds3_Mat <- assay(rld)
rlog_dds3_Mat <- rlog_dds3_Mat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlog_dds3_Mat,file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_AMA_deseq2_rlg_count.txt"),sep="\t", quote=F, row.names=T,col.names=T)
# 只在Linux下能运行，目的是填补表格左上角内容
#system(paste("sed -i '1 s/^/ID\t/'", "ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.xls"))

#step eight: DEGs output
# 定义变量：如果想比较不同的组，只需在这修改即可？能否直接写一个两两比较的循环进行自动处理？
Conditions<-"group"
sampleA <-"AMA"
sampleB <-"Young" 
contrastV <- c(Conditions, sampleA, sampleB)
res1 <- results(dds3, contrast=contrastV) 
#results函数：contrast <- c("conditions", sampleA, sampleB)  提取差异基因分析结果，包含log2 fold changes, p values和adjusted p values. 
#输出的log2FoldChange为log2(SampleA/SampleB)。
#如果直接调用results，只会提取出design matrix最后两个变量的 log2 fold changes and p values等结果。 
#results(dds)
#可以指定比较对象,并可以用mcols查看结果存储的元数据，了解列名的含义
#res1 <- results(dds3) 
#一些总结性的内容
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
res2 <- cbind(baseMeanA, baseMeanB, as.data.frame(res1))
head(res2) 

# 增加ID信息
res3 <- cbind(ID=rownames(res2), as.data.frame(res2))
#res3$baseMean2 <- rowMeans(cbind(baseA, baseB))
head(res3)
# 校正后p-value为NA的复制为1
res3$padj[is.na(res3$padj)] <- 1
res3$pvalue[is.na(res3$pvalue)] <- 1
# 按pvalue排序, 把差异大的基因放前面
res3 <- res3[order(res3$pvalue,decreasing = T),]
head(res3);tail(res3) ;dim(res3)
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_",sampleA,"_vs_",sampleB,".gene_information_unselect.txt")
write.table(as.data.frame(res3), file=file, sep="\t", quote=F, row.names=T,col.names=T)

#whole RNA_data combined
head(as.data.frame(res3));dim(as.data.frame(res3))
head(normalized_counts_log2);dim(normalized_counts_log2)
merge_data <- merge(normalized_counts_log2,as.data.frame(res3),by="row.names",sort=FALSE)
merge_data <- merge_data[order(merge_data$pvalue,decreasing = T),]
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
write.table(merge_data, file=file, sep="\t", quote=F, row.names=F,col.names=T)

table(merge_data$pvalue<0.05)#400 
table(merge_data$pvalue<0.01)# 97 
table(merge_data$pvalue<0.05 & abs(merge_data$log2FoldChange)>=0.585)#364  
table(merge_data$padj<0.1)# 3 

#差异基因的进一步筛选
# padj<0.05 Foldchange >=1.5
res_de <- subset(merge_data, merge_data$padj<0.05, select=c('ID', 'log2FoldChange',"pvalue",'padj'))
res_de<-res_de[order(res_de$padj, decreasing = F), ]
res_de_up <- subset(res_de, res_de$log2FoldChange>=0.585)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*0.585)
DEG_information<-rbind(res_de_up,res_de_dw)
file <-  paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_padj005_FC1.5.txt")
write.table(as.data.frame(DEG_information), file=file, sep="\t", quote=F, row.names=T, col.names=T)

# padj<0.1 & Foldchange >=1.5
res_de <- subset(merge_data, merge_data$padj<0.1, select=c('ID', 'log2FoldChange',"pvalue", 'padj'))
res_de<-res_de[order(res_de$padj, decreasing = F), ]
res_de_up <- subset(res_de, res_de$log2FoldChange>=0.585)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*0.585)
DEG_information<-rbind(res_de_up,res_de_dw);dim(DEG_information)
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_padj01_FC1.5.txt")
write.table(as.data.frame(DEG_information), file=file, sep="\t", quote=F, row.names=T, col.names=T)

# p_value<0.01 & Foldchange >=1.5
res_de <- subset(merge_data, merge_data$pvalue <0.01, select=c('ID','log2FoldChange',"pvalue", 'padj'))
res_de<-res_de[order(res_de$padj, decreasing = F), ]
res_de_up <- subset(res_de, res_de$log2FoldChange>=0.585)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*0.585)
DEG_information<-rbind(res_de_up,res_de_dw)
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_pvalue001_FC1.5.txt")
write.table(as.data.frame(DEG_information), file=file, sep="\t", quote=F, row.names=T, col.names=T)

# p_value<0.05 & Foldchange >=1.5
res_de <- subset(merge_data, merge_data$pvalue <0.05, select=c('ID','log2FoldChange',"pvalue",'padj'))
res_de<-res_de[order(res_de$padj, decreasing = F), ]
res_de_up <- subset(res_de, res_de$log2FoldChange>=0.585)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*0.585)
DEG_information<-rbind(res_de_up,res_de_dw)
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
write.table(as.data.frame(DEG_information), file=file, sep="\t", quote=F, row.names=T, col.names=T)

#step nine: plot for DEGs
#raw plot for evaluation
#Counts plot
#最显著基因
topGene <- as.character(merge_data[which.min(merge_data$padj),]$ID) 
plotCounts(dds3, gene = topGene, intgroup=c("group"))
plotCounts(dds3, gene = topGene, intgroup=c("group"), returnData=TRUE) %>% ggplot(aes(group, count)) + 
  geom_boxplot(aes(fill=group)) +
  geom_jitter(width = 0.2)+scale_y_log10() + ggtitle(paste0(topGene,": (min padj)"))
#MA-plot（也叫mean-difference plot,Bland-Altman plot）了解模型（如所有基因在不同处理比较的结果）的估计系数的分布
#plotMA(res1, ylim = c(-100, 100),text(paste("up_genes:",nrow(res_de_up_id),sep="")))
plotMA(res1)

topGene <- rownames(res1)[which.min(res1$padj)]
with(res1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

res3.shrink <- lfcShrink(dds3, contrast = contrastV, res=res1)
plotMA(res3.shrink, ylim = c(-5,5))
topGene <- rownames(res3.shrink)[which.min(res3.shrink$padj)]
with(res3.shrink[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#Volcano
logCounts <- log2(merge_data$baseMean+1)
logFC <- merge_data$log2FoldChange
FDR <- merge_data$padj
plot(logFC, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", "black"),xlab="logFC", ylab="-1*log1o(FDR)", main="Volcano plot", pch= 20)

#formal Volcano
#png(filename=paste(comp,"Volcano.png", sep=".")) 

res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
nrow(res_de_up)#204
nrow(res_de_dw)#160

filename<-paste(sampleA,"_versus_",sampleB, sep="")
#for pvalue
cols <- rep("#000000", nrow(merge_data))
cols[merge_data$log2FoldChange >= 0.585 & merge_data$pvalue < 0.05] <- "#CC0000"
cols[merge_data$log2FoldChange <= -0.585 & merge_data$pvalue < 0.05] <- "#0000CC"
x_lim<-max(merge_data$log2FoldChange,-merge_data$log2FoldChange)
y_lim<-max(-log10(na.omit(merge_data$pvalue)),log10(na.omit(merge_data$pvalue)))
plot(merge_data$log2FoldChange, -log10(merge_data$pvalue),main=filename, pch=16, cex=0.75,col=cols,las=1,xlim=c(-x_lim,x_lim), ylim=c(0,y_lim+2))
mtext(paste("up_gene: ",nrow(res_de_up)," & down_gene: ",nrow(res_de_dw)," [pvalue <0.05 & FoldChange > 1.5]",sep="") , side=3 , line = 0.3 , outer = FALSE,cex=0.8,col="black") 
#添加辅助线
abline(h=-1*log10(0.05),lwd=3,lty=3,col="#4C5B61")
abline(v=log2(1.5) ,lwd=3,lty=3,col="#4C5B61")
abline(v=log2(2/3) ,lwd=3,lty=3,col="#4C5B61")

#for padj<0.05
res_de_up <- subset(merge_data, padj<0.05 & log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
res_de_dw <- subset(merge_data, padj<0.05 & log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
nrow(res_de_up)
nrow(res_de_dw)
cols <- rep("#000000", nrow(merge_data))
cols[merge_data$log2FoldChange >= 0.585 & merge_data$padj < 0.05] <- "#CC0000"
cols[merge_data$log2FoldChange <= -0.585 & merge_data$padj < 0.05] <- "#0000CC"
x_lim<-max(merge_data$log2FoldChange,-merge_data$log2FoldChange)
y_lim<-max(-log10(na.omit(merge_data$padj)),log10(na.omit(merge_data$padj)))
plot(merge_data$log2FoldChange, -log10(merge_data$padj),main=filename, pch=16, cex=0.75,col=cols,las=1,xlim=c(-x_lim,x_lim), ylim=c(0,y_lim+2))
mtext(paste("up_gene: ",nrow(res_de_up)," & down_gene: ",nrow(res_de_dw),"  [padj<0.05 & FoldChange > 1.5]",sep="") , side=3 , line = 0.3 , outer = FALSE,cex=1,col="black") 
#添加辅助线
abline(h=-1*log10(0.05),lwd=3,lty=3,col="#4C5B61")
abline(v=log2(3/2) ,lwd=3,lty=3,col="#4C5B61")
abline(v=log2(2/3) ,lwd=3,lty=3,col="#4C5B61")

#添加统计值
sig <- merge_data[which(merge_data$padj < 0.05),]
topGene <- as.character(merge_data[which.min(merge_data$padj),]$ID)
with(sig[which(sig$ID == topGene), ], {
  points(log2FoldChange,-log10(padj),col="grey",cex=2, lwd=2)
  text(log2FoldChange,-log10(padj),paste(topGene,sep=""), pos=2, col="black")
  text(log2FoldChange, -log10(padj),"most significant", pos=1,cex=0.7,col="black")
})

topGene <- as.character(merge_data[which.max(merge_data$log2FoldChange),]$ID)
with(merge_data[which(merge_data$ID == topGene), ], {
  points(log2FoldChange, -log10(padj), col="grey", cex=2, lwd=2)
  text(log2FoldChange, -log10(padj), topGene, pos=2, col="green")
  text(log2FoldChange, -log10(padj),"most UP_change", pos=1,  cex=0.7,col="green")
})

topGene <- as.character(merge_data[which.min(merge_data$log2FoldChange),]$ID)
with(merge_data[which(merge_data$ID == topGene), ], {
  points(log2FoldChange, -log10(padj), col="grey", cex=2, lwd=0.5)
  text(log2FoldChange, -log10(padj), topGene, pos=4, col="orange")
  text(log2FoldChange, -log10(padj),"most down_change", pos=1,  cex=0.7,col="orange")
})

##基因聚类的热图展示。
#前20个样本件差异比较大，然后看他们在不同样本间的表达情况。
library("genefilter")
rownames(merge_data)<-merge_data$Row.names
merge_data<-merge_data[,-1]
merge_data[1:6,1:6]
draw_count<-merge_data[,1:10]
#sample_list<-c("P22","P23","P24","P25","P37","P38","P31","P32","P33","P34","P36","P41")
#draw_count<-draw_count[,sample_list]

topVarGenes <- head(order(rowVars(draw_count), decreasing = TRUE), 20)
#rowVars(x, ...): Row variance and standard deviation of a numeric array
#order返回的这个向量刚好是xx向量中各数的原始的位置
mat  <- draw_count[topVarGenes, ]
#mat  <- mat - rowMeans(mat)
anno1 <- as.data.frame(colData(dds2)[, c("group","group_type")])
pheatmap(mat,annotation_col = anno1,cluster_rows=T, cluster_cols =F,)
pheatmap(mat, cluster_rows=T, cluster_cols =F,scale="row", annotation_col=anno1)

#找样本间全部差异显著基因在不同样本间的表达情况--类似计算
resSig_RT <- subset(merge_data, pvalue< 0.05)
nrow(resSig_RT)
resSig_RT<-resSig_RT[order(resSig_RT$pvalue, decreasing = F), ]
namm<-rownames(draw_count)
mat  <- draw_count[which(namm %in% rownames(resSig_RT)), ]
#mat  <- mat - rowMeans(mat)
#pheatmap(mat, cluster_row=T, scale="none", annotation_col=anno1,main="scale_Genes_padjust<0.05")
pheatmap(mat,cluster_row=T,cluster_cols =F, scale= "row", annotation_col=anno1,main="scale_Genes_pvalue<0.05")

#绘制各组特异高表达基因热图
res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
nrow(res_de_up )
nrow(res_de_dw )

res_de_up_id = data.frame(ID=res_de_up$ID,type=paste(sampleA,"_higherThan_", sampleB, sep=""))
res_de_dw_id = data.frame(ID=res_de_dw$ID,type=paste(sampleA,"_lowerThan_", sampleB, sep=""))
de_id_whole = rbind(res_de_up_id, res_de_dw_id)
nrow(de_id_whole)
#by normal count
DEG_expr <- draw_count[which(rownames(draw_count) %in% unique(de_id_whole$ID)),]
pheatmap(DEG_expr, cluster_rows=T, cluster_cols =F,scale="row", annotation_col=anno1)
pheatmap(DEG_expr, cluster_rows=T, cluster_cols =F,scale="row", annotation_col=anno1,show_rownames = F,main="scale_Genes_pvalue<0.05 & FC >1.5")

#上调与下调最明显的基因分别选取20个，然后看他们在不同样本间的表达情况--normalized_counts
res_de_up <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
res_de_dw <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
nrow(res_de_up )
nrow(res_de_dw )
res_de_up_top20_id <-  rownames(res_de_up[head(order(res_de_up$padj,decreasing = F),20),])
res_de_dw_top20_id <-  rownames(res_de_dw[head(order(res_de_dw$padj,decreasing = F),20),])
red_de_top20 <- c(res_de_up_top20_id, res_de_dw_top20_id)
red_de_top20 
red_de_top20_expr <- draw_count[rownames(draw_count) %in% red_de_top20,]
pheatmap(red_de_top20_expr, cluster_row=T,cluster_cols =T, scale="row", annotation_col=anno1,main="Top 20 UP/Down DEGs(pj0.05&FC1.5)")


#热图修饰
#查看数据范围：
mat<-DEG_expr #padjust <0.01&log2FoldChange>=0.585
range(mat)# 0.00000 18.32706
#使用上面的矩阵做热图，要求低值为蓝色，高值为红色，中间值为白色：
pheatmap(mat,
         scale =  "row",
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         annotation_col = anno1[2],main="Genes_pvalue <0.05&log2FoldChange>=0.585"
)
#数据data1的数值范围是-4到8.6，所以上图的色条范围也是-4到8.6。
#但由于0的位置比较特殊，现在要求0的位置为白色，并且色条范围为-9到9，
#这里使用breaks参数重新定义色条范围并根据break范围划分颜色范围，代码如下：
#breaks 
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01)) 
# 做热图： 
heat_plot<-pheatmap(mat, 
         scale = "row", 
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), 
         #         color = c(colorRampPalette(colors = c("blue","black"))(length(bk)/3),colorRampPalette(colors = c("black","red"))(length(bk)/3)), 
         cluster_cols=F,
         legend_breaks=seq(-2,2,1), 
         breaks=bk,
         show_rownames = F,
         annotation_col = anno1,main="Genes_pvalue <0.05& FoldChange>=1.5"
)

pdf(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_heatmap_DEGs_p005_FC1.5_AMA_Young.pdf"),width = 7, height =7)
print(heat_plot)
dev.off()
#进一步美化：https://shengxin.ren/article/107

#step ten: perform enrichment_AMA_nobatch 
#file <- paste0("D:/3.秦萌高龄甲基化/1.RNA_result/merge_count/",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
#merge_data<-read.table(file=file, sep="\t",header = T,row.names=1)
#for padj<0.05& FoldChange>=1.5
res_de_up_padj005_FC1.5 <- subset(merge_data, padj<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
res_de_dw_padj005_FC1.5 <- subset(merge_data, padj<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))

#for padj<0.1& FoldChange>=1.5
res_de_up_padj01_FC1.5 <- subset(merge_data, padj<0.1&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
res_de_dw_padj01_FC1.5 <- subset(merge_data, padj<0.1&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))

#for pvalue<0.01& FoldChange>=1.5
res_de_up_pvalue001_FC1.5 <- subset(merge_data, pvalue<0.01&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'pvalue'))
res_de_dw_pvalue001_FC1.5 <- subset(merge_data, pvalue<0.01&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'pvalue'))
#for pvalue<0.05& FoldChange>=1.5
res_de_up_pvalue005_FC1.5 <- subset(merge_data, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'pvalue'))
res_de_dw_pvalue005_FC1.5 <- subset(merge_data, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'pvalue'))

gene_list<-c(list(res_de_up_padj005_FC1.5),list(res_de_dw_padj005_FC1.5),
             list(res_de_up_padj01_FC1.5),list(res_de_dw_padj01_FC1.5),
             list(res_de_up_pvalue001_FC1.5),list(res_de_dw_pvalue001_FC1.5),
             list( res_de_up_pvalue005_FC1.5),list(res_de_dw_pvalue005_FC1.5))
list_name<-c("up_padj005_FC1.5","dw_padj005_FC1.5",
             "up_padj01_FC1.5","dw_padj01_FC1.5",
             "up_pvalue001_FC1.5","dw_pvalue001_FC1.5",
             "up_pvalue005_FC1.5","dw_pvalue005_FC1.5")
gene_list_used<-gene_list[c(7,8)]
list_name_used<-list_name[c(7,8)]
for ( i in 1:length(list_name_used)){
  #i=2
  gene<-as.character(gene_list_used[[i]]$ID)
  print(as.character(list_name_used[[i]]))
  gene.df <- bitr(gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  target<-unique(gene.df$ENTREZID)
  if(length(target) ==0 ){next;}
  gene_GO_RPEA_erichment_results<-list()
  BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T)  
  MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
  CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
  head(as.data.frame(BP@result))
  CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
  BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
  MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
  head(as.data.frame(BP_simp@result))
  ## Reactome pathway enrichment_AMA_nobatch analysis
  rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
  head(as.data.frame(rpea@result))
  
  write.table(as.data.frame(CC@result), file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_",as.character(list_name_used[[i]]),"_CC_GO.txt"), row.names=T, col.names=T) 
  write.table(as.data.frame(BP@result), file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_",as.character(list_name_used[[i]]),"_BP_GO.txt"),row.names=T, col.names=T) 
  write.table(as.data.frame(MF@result), file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_",as.character(list_name_used[[i]]),"_MF_GO.txt"),row.names=T, col.names=T) 
  write.table(as.data.frame(CC_simp@result), file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_",as.character(list_name_used[[i]]),"_BP_GO_simp.txt"),row.names=T, col.names=T) 
  write.table(as.data.frame(MF_simp@result), file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_",as.character(list_name_used[[i]]),"_MF_GO_simp.txt"), row.names=T, col.names=T) 
  head(as.data.frame(BP_simp))
  write.table((as.data.frame(rpea@result)), file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_",as.character(list_name_used[[i]]),"_RPEA.txt"),row.names=T, col.names=T) 
  
  kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
  if(is.null(kk)){gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,rpea)
  }else{
    kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
    head(as.data.frame(kk@result))
    write.table((as.data.frame(kk@result)), file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_",as.character(list_name_used[[i]]),"_KEGG.txt"),row.names=T, col.names=T) 
    gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
  }
  saveRDS(gene_GO_RPEA_erichment_results, file = paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_",as.character(list_name_used[[i]]),"_GO_RPEA_result.rds"))
}

#step eleven: plot target  enrichment_AMA_nobatch result
# https://www.jianshu.com/p/e133ab3169fa 
gene_erichment_results_Up<-readRDS(file = paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_up_pvalue005_FC1.5_GO_RPEA_result.rds"))
gene_erichment_results_Down<-readRDS(file = paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_dw_pvalue005_FC1.5_GO_RPEA_result.rds"))

gene_erichment_results<-gene_erichment_results_Up #gene_erichment_results_hyper_unselect
CC_plot<-gene_erichment_results[[1]]
BP_plot<-gene_erichment_results[[2]]
MF_plot<-gene_erichment_results[[3]]

#结果展示
merge_three<-as.data.frame(rbind(head(BP_plot@result,n=10),head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
merge_three$qvalue_log10<-c(-log(merge_three$qvalue,10))
merge_three$pvalue_log10<-c(-log(merge_three$pvalue,10))

merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))

pdf(paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_Up_BP_CC_MF.pdf"))
ggbarplot(merge_three, x="Description", y="qvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=0)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55)) #按组排序 x.text.angle=90
#修进：http://blog.sciencenet.cn/blog-3334560-1091714.html
ggbarplot(merge_three, x="Description", y="pvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          #          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=60)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55)) #按组排序 x.text.angle=90
#修进：http://blog.sciencenet.cn/blog-3334560-1091714.html
dev.off()

#For hypO
gene_erichment_results<-gene_erichment_results_Down#gene_erichment_results_hypo_unselct
CC_plot<-gene_erichment_results[[1]]
BP_plot<-gene_erichment_results[[2]]
MF_plot<-gene_erichment_results[[3]]

#结果展示
merge_three<-as.data.frame(rbind(head(BP_plot@result,n=10),head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
dim(merge_three)
#merge_three$qlogtran<-c(-log(merge_three$qvalue,10))
#merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))
#merge_two<-as.data.frame(rbind(head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
merge_three$qvalue_log10<-c(-log(merge_three$qvalue,10))
merge_three$pvalue_log10<-c(-log(merge_three$pvalue,10))
merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))

pdf(paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_Down_BP_CC_MF.pdf"))
ggbarplot(merge_three, x="Description", y="qvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=0)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55)) #按组排序 x.text.angle=90
#修进：http://blog.sciencenet.cn/blog-3334560-1091714.html
ggbarplot(merge_three, x="Description", y="pvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          #          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=60)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55))
dev.off()

KEGG_plot_Up<-gene_erichment_results_Up[[7]] #gene_erichment_results_hyper_unselect
KEGG_plot_Down<-gene_erichment_results_Down[[7]] #gene_erichment_results_hypo_unselect

dim(KEGG_plot_Up@result)
dim(KEGG_plot_Down@result)

KEGG_plot<-as.data.frame(rbind(head(KEGG_plot_Up@result,n=10),head(KEGG_plot_Down@result,n=10)))
KEGG_plot$group<-c(rep("Up",10),rep("Down",10))
str(KEGG_plot)
KEGG_plot$qvalue_log10<-c(-log(KEGG_plot$qvalue,10))
KEGG_plot$pvalue_log10<-c(-log(KEGG_plot$pvalue,10))

#BP_plot$group<- factor(x =BP_plot$group,levels =c("Trobt_Vims","EVTs","VCTs","STBs") )
KEGG_plot$group<- factor(x =KEGG_plot$group,levels =c("Up","Down") )
levels(KEGG_plot$group)
KEGG_plot<-KEGG_plot[,c("Description","group","qvalue_log10","Count")]
KEGG_plot$Description<-as.character(KEGG_plot$Description)
KEGG_plot<-KEGG_plot[order(KEGG_plot$group,KEGG_plot$qvalue_log10),]
str(KEGG_plot)
head(KEGG_plot)

PP_plot<-KEGG_plot
head(PP_plot)
pdf(paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_Up_Down_gene_KEGG.pdf"))

p <- ggplot(PP_plot,aes(qvalue_log10 ,reorder(Description,qvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(p.adjust)",y="Description",title="KEGG enrichment_AMA")+
  scale_size_continuous(range=c(5,14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=55))+
  geom_text(aes(label=sprintf("%.3f",qvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")
p
KEGG_plot<-as.data.frame(rbind(head(KEGG_plot_Up@result,n=10),head(KEGG_plot_Down@result,n=10)))
KEGG_plot$group<-c(rep("Up",10),rep("Down",10))
str(KEGG_plot)
KEGG_plot$qvalue_log10<-c(-log(KEGG_plot$qvalue,10))
KEGG_plot$pvalue_log10<-c(-log(KEGG_plot$pvalue,10))

#BP_plot$group<- factor(x =BP_plot$group,levels =c("Trobt_Vims","EVTs","VCTs","STBs") )
KEGG_plot$group<- factor(x =KEGG_plot$group,levels =c("Up","Down") )
levels(KEGG_plot$group)

KEGG_plot<-KEGG_plot[,c("Description","group","pvalue_log10","Count")]
KEGG_plot$Description<-as.character(KEGG_plot$Description)
KEGG_plot<-KEGG_plot[order(KEGG_plot$group,KEGG_plot$pvalue_log10),]
str(KEGG_plot);head(KEGG_plot)

PP_plot<-KEGG_plot
head(PP_plot)
ggplot(PP_plot,aes(pvalue_log10 ,reorder(Description,pvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(pvalue)",y="Description",title="KEGG enrichment_AMA")+
  scale_size_continuous(range=c(5,14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=55))+
  geom_text(aes(label=sprintf("%.3f",pvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")
dev.off()

#for no selected GO terms of DEGs
GO_BP_up<-gene_erichment_results_Up[[2]]
GO_BP_down<-gene_erichment_results_Down[[2]]

Enrich_all_plot<-as.data.frame(rbind(GO_BP_up@result,GO_BP_down@result))
Enrich_all_plot$group<-c(rep("DEGs_Up",nrow(GO_BP_up@result)),rep("DEGs_Down",nrow(GO_BP_down@result)))
Enrich_all_plot<-Enrich_all_plot[which(Enrich_all_plot$pvalue<0.05 & Enrich_all_plot$Count >1),]
dim(Enrich_all_plot)# 433  10
Enrich_all_plot$Log10_padjust<- c(-log10(Enrich_all_plot$p.adjust))
Enrich_all_plot$Log10_pvalue<- c(-log10(Enrich_all_plot$pvalue))
range(Enrich_all_plot$Log10_pvalue)#1.301309 4.679950
range(Enrich_all_plot$Log10_padjust)#0.304503 1.404408
range(Enrich_all_plot$Count)#  2 11

Enrich_all_number <- data.frame(table(as.character(Enrich_all_plot$Description)))
#Enrich_all_number2 <- Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot2<-merge(Enrich_all_plot,Enrich_all_number)
write.table(Enrich_all_plot2, file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_GO_BP_merge.txt"),row.names=T, col.names=T) 

Enrich_all_plot_noselect<-read.table(paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_GO_BP_merge.txt"),header =T) 
head(Enrich_all_plot_noselect)
length(unique(Enrich_all_plot_noselect$Description))#417 
length(unique(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>=2),]$Description))#255 
Enrich_all_plot_noselect$group<-factor(Enrich_all_plot_noselect$group,levels = c("DEGs_Down","DEGs_Up"))
Enrich_all_plot_noselect <- Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$group),]
Enrich_all_plot_noselect<-Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$freq,decreasing = F),]
Enrich_all_plot_noselect$Description<-factor(Enrich_all_plot_noselect$Description,levels = c(unique(as.character(Enrich_all_plot_noselect$Description))))

range(Enrich_all_plot_noselect$pvalue) #2.089536e-05 4.996786e-02
range(Enrich_all_plot_noselect$Log10_pvalue) #1.301309 4.679950
range(Enrich_all_plot_noselect$Count) #3 51

GO_BP_plot<-ggplot(Enrich_all_plot_noselect,aes(x=group,y=Description,size=Count,colour=Log10_pvalue))+
  geom_point(alpha =0.8,na.rm = TRUE)+
  scale_size(breaks = c(0,2,4,6,8,10),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(1,2,3,4,5),name='-log10(p_value)')+  
  scale_y_discrete(labels=function(x) str_wrap(x, width=150))+
  theme_classic()+labs(x="",y="GO terms",title="AMA related DEGs BP enrichment")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),
        legend.position ="right",legend.direction = "vertical")
write.table(Enrich_all_plot_noselect, file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_wait_selected_GO_BP_merge.txt"),row.names=T, col.names=T) 
ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_wait_selected_GO_BP.pdf"),GO_BP_plot,width =15,height =45)

#for no selected kegg pathway of DEGs
KEGG_up<-gene_erichment_results_Up[[7]]
KEGG_down<-gene_erichment_results_Down[[7]]

Enrich_all_plot<-as.data.frame(rbind(KEGG_up@result,KEGG_down@result))
Enrich_all_plot$group<-c(rep("DEGs_Up",nrow(KEGG_up@result)),rep("DEGs_Down",nrow(KEGG_down@result)))
dim(Enrich_all_plot)
Enrich_all_plot$Log10_pvalue<- c(-log10(Enrich_all_plot$pvalue))
range(Enrich_all_plot$Log10_pvalue)#0.02578149 3.33675597
Enrich_all_plot<-Enrich_all_plot[which(Enrich_all_plot$pvalue<0.05 & Enrich_all_plot$Count >1),]
range(Enrich_all_plot$Count)#2 7

Enrich_all_number <- data.frame(table(as.character(Enrich_all_plot$Description)))
#Enrich_all_number2 <- Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot2<-merge(Enrich_all_plot,Enrich_all_number)
write.table(Enrich_all_plot2, file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_KEGG_merge.txt"),row.names=T, col.names=T) 

Enrich_all_plot_noselect<-read.table(paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_KEGG_merge.txt"),header =T) 
head(Enrich_all_plot_noselect)
length(unique(Enrich_all_plot_noselect$Description))#15 
length(unique(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>1),]$Description))#104 
Enrich_all_plot_noselect$group<-factor(Enrich_all_plot_noselect$group,levels = c("DEGs_Down","DEGs_Up"))
Enrich_all_plot_noselect <- Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$group),]
Enrich_all_plot_noselect<-Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$freq,decreasing = F),]
Enrich_all_plot_noselect$Description<-factor(Enrich_all_plot_noselect$Description,levels = c(unique(as.character(Enrich_all_plot_noselect$Description))))

range(Enrich_all_plot_noselect$Log10_pvalue) #1.337537 3.336756
range(Enrich_all_plot_noselect$Count) # 1 23
KEGG_plot<-ggplot(Enrich_all_plot_noselect,aes(x=group,y=Description,size=Count,colour=Log10_pvalue))+
  geom_point(alpha =0.8,na.rm = TRUE)+scale_size(breaks = c(0,2,4,6,8),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,1,2,3,4),name='-log10(p_value)')+  
  scale_y_discrete(labels=function(x) str_wrap(x, width=150))+
  theme_classic()+labs(x="",y="GO terms",title="AMA related DEGs KEGG enrichment")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),
        legend.position ="right",legend.direction = "vertical")
write.table(Enrich_all_plot_noselect, file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_wait_selected_KEGG_merge.txt"),row.names=T, col.names=T) 
ggsave(file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/",tag_name,"_enrichment_",sampleA,"_vs_",sampleB,"_wait_selected_KEGG.pdf"),KEGG_plot,width =15,height =20)

}

###=================待看
merge_three<-as.data.frame(rbind(head(BP@result,n=10),head(CC@result,n=10),head(MF@result,n=10)))
merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))
ggbarplot(merge_three, x="Description", y="pvalue", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序 
          sort.by.groups=TRUE,x.text.angle=60) #按组排序 x.text.angle=90

##具体绘图
GO_plot_data<-BP
#泡泡图
dotplot(GO_plot_data,showCategory = 14, color="pvalue", font.size=18) 
#条形图
barplot(GO_plot_data, showCategory=15,title="enrichment_AMA_nobatchGO_ALL")#条状图，按p从小到大排的#（这个图...好丑)

#绘制富集网格图 enrichment_AMA_nobatch Map 可以将富集条目和重叠的基因集整合为一个网络图，相互重叠的基因集则趋向于成簇，从而易于分辨功能模型。
#对于富集到的GO terms之间的基因重叠关系进行展示，如果两个GO terms系的差异基因存在重叠，说明这两个节点存在overlap关系，在图中用线条连接起来
emapplot(GO_plot_data, showCategory = 10)
#函数 goplot() 接受 enrichGO() 的输出，并将其可视化。
goplot(GO_plot_data)
#还可以绘制GO的网络关系图，但是值得注意的是这里的数据只能是富集一个GO通路（BP、CC或MF）的数据
plotGOgraph(GO_plot_data)

#函数 cnetplot() 可以将基因与生物学概念 (e.g.* GO terms or KEGG pathways) 的关系绘制成网状图。
#绘制各过程中各个元素的表达变化量
cnetplot(GO_plot_data, foldChange=gene$log2FoldChange)
cnetplot(GO_plot_data, categorySize="pvalue", foldChange=gene$log2FoldChange)  ## categorySize 可以是 "pvalue" 或 "geneNum" 
cnetplot(GO_plot_data, foldChange=gene$log2FoldChange, circular = TRUE, colorEdge = TRUE)  ## 圆形布局，给线条上色

#Heatmap 热图能够简化结果，更容易分辨表达模式 (expression patterns) 。
heatplot(GO_plot_data, gene$log2FoldChange)
#UpSet Plot 函数 upsetplot() 是 cneplot() 的一种替代方案，用于可视化基因与基因集间的复杂关联，而 upsetplot() 更着重于不同基因集间基因的重叠情况。
upsetplot(GO_plot_data)

#绘制文本所属节点图
plotGOgraph(GO_plot_data)

#GOplot包是用于GO富集结果可视化的R包，尤其其中的弦图
#https://bioinfohome.com/index.php/2019/07/08/goplot-01/
#GOplot输入文件，包含三个，GO富集结果、基因列表（包含logFC）、GO条目。
#如果，我们处理的数据不包含logFC，或者不想使用logFC。我们可以自己生产随机数代替。
#最后的pdf文件，借助AI绘制软件，将logFC的图标删除即可得到漂亮的弦图。
library(GOplot)
library(stringr)

#构建GOplot的输入文件。
#GOplot的输入文件包含terms，genes，process。
#terms： 是一个数据框，包含的字段：Category,ID,term,Genes,adj_pval；
#genes：基因名加logFC（ logFC 随机生成）；
#process：GO的条目描述；
ego1<-as.data.frame(ego@result)
GO=ego1[1:12,c(2,1,2,8,6)]
##因为GO富集结果非常多，不可能也不需要全部显示，一般选取前几行，这里取12，也可以设定为其他值
GO$geneID=str_replace_all(GO$geneID,"/",",")
#这里替换富集条目中基因一栏的分隔符为逗号，这是GOplot默认的。
names(GO)=c("Category","ID","term","Genes","adj_pval")#重新命名数据框的字段名，这些都是固定的，不可更改
#提供带有基因名加logFC的数据框，这点和GSEA的要求有点像，此处自己构建
gene=data.frame(ID=gene.df$SYMBOL,logFC=rnorm(length(gene.df$SYMBOL),mean=0,sd=2))
#生成输入数据框 #随机生成相同数量的正负小数，充作logFC

#第四步，使用GOplot作出弦图，并导出pdf，使用AI，微调
circ <- circle_dat(GO,gene)
chord <- chord_dat(data = circ, genes = gene, process = GO$term)

#GO$term就是富集的GO条目描述
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
# Excluding process with less than 5 assigned genes
GOChord(chord, limit = c(0,5))
# Creating the chord plot genes ordered by logFC and a different logFC color scale
GOChord(chord,space=0.02,gene.order='logFC',lfc.col=c('red','black','cyan'))

##KEGG富集

kk <- enrichKEGG(gene = gene$ID,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 1,use_internal_data =FALSE)

file <- paste("H:/ART_Epigenome_Papers/ART项目数据/wholeBlood_RNAseq_20180916/DEseq2_20190809_by_CW/",sampleA,"_vs_",sampleB,sep="")
write.table(as.data.frame(CC_simp@result), file=paste(file,"hyper_KEGG.txt",sep="."),quote=F, row.names=F, col.names=T) 

dotplot(kk,showCategory = 14,colorBy="pvalue",font.size=18)

## Reactome pathway enrichment_AMA_nobatch analysis
require(ReactomePA)
rpea <- enrichPathway(gene=gene$ID,pvalueCutoff=0.05, readable=T,organism = "human")
head(as.data.frame(rpea))

file <- paste("H:/ART_Epigenome_Papers/ART项目数据/wholeBlood_RNAseq_20180916/DEseq2_20190809_by_CW/",sampleA,"_vs_",sampleB,sep="")
write.table(as.data.frame(CC_simp@result), file=paste(file,"hyper_RPEA.txt",sep="."),quote=F, row.names=F, col.names=T) 


#visulization for pathway
pathway_data<-kk
barplot(pathway_data,showCategory=10)
dotplot(pathway_data,showCategory=10)
emapplot(pathway_data)  #enrichment_AMA_nobatch map
cnetplot(pathway_data, categorySize="pvalue", foldChange=gene) 
heatplot(pathway_data, gene$log2FoldChange)

#browseKEGG
#函数 browseKEGG 可以帮你打开浏览器，嗯。
browseKEGG(pathway_data, 'hsa04110')

###多个基因集时（eg，WGCNA得到的不同module的基因集），用compareCluster。
#参考：https://www.jianshu.com/p/c01b4cc1b98a




##GSEA分析 Gene Set enrichment_AMA_nobatch Analysis(GSEA) 参考：https://www.jianshu.com/p/feaefcbdf986
#GSEA(geneList, exponent = 1, nPerm = 1000, minGSSize = 10,maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE,TERM2NAME = NA, verbose = TRUE, seed = FALSE, by = "fgsea")

#1 获取按照log2FC大小来排序的基因列表
sig.gene<-res_de #res_de
##排序--建议先转为ENTREZID
#数据构建
Tans_genes <- bitr(sig.gene$ID, fromType = "REFSEQ", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
rownames(Tans_genes)<-Tans_genes$REFSEQ
sig_gene<-merge(sig.gene,Tans_genes,by= "row.names")

#排序
genelist_frame<-sig_gene[,c("ENTREZID","log2FoldChange")]
genelist_frame <- arrange(genelist_frame, desc(log2FoldChange))
gene<-geneList.sort$ENTREZID
gene[duplicated(gene)]

#2 这里使用的是broad GSEA提供的gene sets 来提供TERM2GENE：
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
head(c5)
##选择GSEA数据库用的，http://software.broadinstitute.org/gsea/msigdb/collections.jsp，
#前面俩行大概就是来选择GSEA的数据库用的，c5.cc.v5.0.entrez.gmt这里他事例用的是GO的cc数据库，这个可以调整的，参考GSEA的网站来选择

##3 先使用基于超几何分布的富集分析
enrich <- enricher(gene, TERM2GENE=c5); head(enrich)

#用于GSEA的list构建
genelist <- genelist_frame$log2FoldChange
names(genelist) <- as.character(genelist_frame$ENTREZID)
genelist <- sort(genelist, decreasing = TRUE); head(genelist)

#4输入文件准备好了尽可以进行GSEA富集分析了：
gsea <- GSEA(genelist, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.1); head(gsea)

# clusterProfiler的GSEA分析
#gseGO(geneList, ont = "BP", OrgDb, keyType = "ENTREZID", exponent = 1,
#      nPerm = 1000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
#      pAdjustMethod = "BH", verbose = TRUE, seed = FALSE, by = "fgsea")
#主要参数：
#geneList： 排序数据， 可以根据log2foldchange, 也可以是pvalues
#nPerm： 重抽取次数
#minGSSize： 每个基因集的最小数目maxGSSize： 用于测试的基因注释最大数目

#gsemf <- gseGO(genelist,OrgDb = org.Hs.eg.db,pvalueCutoff = 0.5,ont="BP") ; head(gsemf)
gsemf <- gseGO(genelist,OrgDb = org.Hs.eg.db,pvalueCutoff = 0.5) ; head(gsemf)
# 画出GSEA图
gseaplot(gsemf, geneSetID="GO:0008150")

#GO GSEA
BP_GSEA<- gseGO(geneList, OrgDb= org.Hs.eg.db, ont= "BP",
                nPerm = 1000,minGSSize= 100,maxGSSize = 500,pvalueCutoff = 0.05,verbose = FALSE)  ## 不输出结果
#KEGG GSEA
kk_GSEA <- gseKEGG(geneList, organism  = 'hsa',
                   nPerm = 1000,minGSSize= 120,pvalueCutoff = 0.05,verbose = FALSE)  ## 不输出结果
head(kk_GSEA,n=3)

# 画出GSEA图
gseaplot(BP_GSEA, geneSetID="GO:0001819")
gseaplot(kk_GSEA, geneSetID="GO:0001819")

##GSEA结果的表达分布叠嶂图 (ridgeline plot for expression distribution of GSEA result)
#更直观地展示上调/下调的通路。
ridgeplot(gsemf)

#利用上面处理好的glist进行一下KEGG富集到的某一条通路的可视化：
library(pathview)
pathview(gene.data = glist, pathway.id = 'hsa04658',species="hsa", limit=list(gene=max(abs(genelist)), cpd=1))
sessionInfo()
