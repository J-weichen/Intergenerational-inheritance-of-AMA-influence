#转录组总的分析流程：https://www.jianshu.com/p/fb15249200d7
##另外可以参考：https://www.cnblogs.com/weipeng-loaded/p/10601285.html
#待学习Limma差异分析、火山图、功能富集： https://www.jianshu.com/p/0a8841715804

#参考：http://www.gene123.com/n/2885/
#数据读取
#如果已经有了整合好的readscount matrixs
rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

setwd("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts")

library("DESeq2") 
library("dplyr")
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

files<- list.files("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/feactureCounts")
#list.files命令将input文件夹下所有文件名输入
#查看下文件的数据格式,分别记录了每个基因的EntrezID，基因长度和数量
read.delim(files[1],header = F, nrow=2,row.names = 1,sep = "\t")
#         V2
#DDX11L1   0
#WASH7P  442
#之后我们需要用DESeqDataSetFromMatrix为DESeq2提供数据,命令形式如下：

#DESeqDataSetFromMatrix(countData, colData, design, tidy = FALSE,ignoreRank = FALSE, ...)
#其中countData存放counts数据，colData存放样本信息的数据，design就是实验设计。
#首先从各个count matrix文件中读取count，基因长度部分可以舍弃，
#因为DESeq2只需要为标准化的count数据，不需要提供基因长度信息。
#逻辑就是分别读取每一个文件的count列，然后赋予样本名。
sample1 <- read.delim(files[1],header = F,row.names = 1,sep = "\t",col.names = c("gene_id",substring(files[1], 1, nchar(files[1])-19)))
head(sample1)
#read.delim 是 read.table 的变形，read_delim 比 read.table要快
count.table <- data.frame(sample1)
head(count.table)
#options(stringsAsFactors = FALSE)

for ( f in files[2:length(files)]){
  #f<- files[2]
  column <- read.delim(f,header=F,row.names = 1,sep = "\t",col.names = c("gene_id",substring(f, 1, nchar(f)-19)))
  head(column)
  #row.names是一个可选参数，用以指定一个或多个表示行标识符的变量
  count.table <- cbind(count.table,column[rownames(column),])
}
head(count.table) 
samplenames <- substring(files, 1, nchar(files)-19)
colnames(count.table) <- samplenames
str(count.table);head(count.table);dim(count.table)#26485    20
count.table<-count.table[which(rowSums(count.table)>0),]
str(count.table);head(count.table);dim(count.table)#20567    20

write.table(count.table,file="/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Verification_AMA_versus_Young_rawcount.xls",sep="\t", quote=F, row.names=T,col.names=T)

count.table2<-read.table(file="/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Verification_AMA_versus_Young_rawcount.xls",sep="\t", row.names=1,header =T)
head(count.table2)

#然后要提供colData,其中colData存放行名要和countData的列名相同。
#colData可以直接提供也可以自己构建
samplenames2<-colnames(count.table2)#将要分析的亚集列名字赋予此
group<-ifelse(grepl("A",samplenames2),"AMA","Young")
group<-as.factor(group)
generation<-ifelse(grepl("M",samplenames2),"Mom","Kid")
generation<-as.factor(generation)
sample<-as.factor(samplenames2)
colData_whole <- data.frame(sample = sample, group=group,generation=generation)
#最后构建数据:
colnames(colData)
#for mom
colData<-colData_whole[which(colData_whole$generation=="Mom"),]
expression_matrix<-count.table2[,as.character(colData$sample)]
dds0 <- DESeqDataSetFromMatrix(countData = expression_matrix, colData = colData, design = ~ group )

#数据预处理
#预过滤#虽然DESeq2会自动屏蔽那些低count的基因，但是剔除那些几乎不表达基因的部分能够提高运行速度。
dds1 <- dds0[ rowSums(counts(dds0)) > 1, ]
dds2 <- estimateSizeFactors(dds1)

#不同的数据转换效果展示【结果就是转换后更加集中了】：简单粗暴的方法就是对count matrix取log后加1；数据集小于30 -> rlog，大数据集 -> VST。
#还有这个处理过程不是用于差异检验的，在DESeq分析中会自动选择最合适的
rld <- rlog(dds1, blind = FALSE)
head(assay(rld), 3)
vsd <- vst(dds1, blind = FALSE)
head(assay(vsd), 3)

df <- bind_rows(
  as_data_frame(log2(counts(dds2, normalized=TRUE)[, 1:2]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)


#样本距离： RNA-Seq分析第一步通常是评估样本间的总体相似度。
#那些样本最接近
#那些样本差异较大
#这与实验设计预期符合么
#这里使用R内置的 dist 计算 rlog变换数据的Euclidean distance，然后用pheatmap根据结果画热图

anno1 <- as.data.frame(colData(dds2)[, c("sample","group")])

sampleDists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "Log2_trans")

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "rld_trans")

sampleDists <- dist(t(assay(vsd)))
sampleDists <- dist(t(log2(counts(dds2, normalized=TRUE)+1)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation_col = anno1, main = "vsd_trans")
#同样的可以用Poisson Distance (Witten 2011)计算距离，计算方式如下：
#library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds2)))
#??怎么画图呢

#另一个可视化样本-样本距离的方法就是主成分分析。
#DESeq2提供了专门的方法用于作图
plotPCA(rld,intgroup=c("group"))
plotPCA(vsd,intgroup=c("group"))


#差异表达基因分析（DEA）
dds3 <- DESeq(dds2)

#标准化后进行样本评估：样品层级聚类分析，判断样品的相似性和组间组内差异
rld3 <- rlog(dds3)
vsd3 <- vst(dds3) 

nor_Mat<-counts(dds3, normalized=TRUE)
rlog_dds3_Mat <- assay(rld3)
vlog_dds3_Mat <- assay(vsd3)

# 层级聚类
hc <- hcluster(t(nor_Mat), method="pearson")
plot(hc, cex = 2) 
hc <- hcluster(t(rlog_dds3_Mat), method="pearson")
plot(hc, cex = 2) 
hc <- hcluster(t(vlog_dds3_Mat), method="pearson")
plot(hc, cex = 2) 
hc <- hcluster(t(nor_Mat), method="spearman")
plot(hc, cex = 2) 
hc <- hcluster(t(rlog_dds3_Mat), method="spearman")
plot(hc, cex = 2) 
hc <- hcluster(t(vlog_dds3_Mat), method="spearman")
plot(hc, cex = 2) 

# 生成颜色
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# 计算相关性pearson correlation
hc <- hcluster(t(rlog_dds3_Mat), method="pearson")
plot(hc, cex = 2) 
#pearson_cor <- as.matrix(cor(nor_Mat, method="pearson"))
pearson_cor <- as.matrix(cor(rlog_dds3_Mat, method="pearson"))
#pearson_cor <- as.matrix(cor(vlog_dds3_Mat, method="pearson"))
heatmap.2(pearson_cor, 
          Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, 
          margins=c(11,11), 
          main="The pearson correlation of each sample")

hc <- hcluster(t(rlog_dds3_Mat), method="spearman")
plot(hc, cex = 2) 
spearman_cor <- as.matrix(cor(rlog_dds3_Mat, method="spearman"))
heatmap.2(spearman_cor, 
          Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, 
          margins=c(11,11), 
          main="The spearman correlation of each sample")

#pdf("ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.pearson.pdf", pointsize=10)
#dev.off()
plotPCA(rld3,intgroup=c("group"))


#获取标准化后的数据
# ?counts查看此函数功能# normalized=T, 返回标准化的数据
##the normalized counts (normalized for library size, which is the total number of gene counts per sample, while accounting for library composition)
normalized_counts <- nor_Mat
head(normalized_counts) 
#根据基因在不同的样本中表达变化的差异程度mad值对数据排序，差异越大的基因排位越前。
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
# 标准化后的数据输出
write.table(normalized_counts, file="/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Verification_AMA_versus_Young_trans.Count.DESeq2.normalized.raw.xls",quote=F, sep="\t", row.names=T, col.names=T)
# log转换后的结果
rlog_dds3_Mat <- rlog_dds3_Mat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlog_dds3_Mat, file="/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Verification_AMA_versus_Young_trans.Count.DESeq2.normalized.raw.rlog.xls",quote=F, sep="\t", row.names=T, col.names=T)
# 只在Linux下能运行，目的是填补表格左上角内容
#system(paste("sed -i '1 s/^/ID\t/'", "ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.xls"))


# 热图绘制R语言与热图分析:https://blog.csdn.net/XIUXIU179/article/details/79975702


#提取差异基因分析结果，构建结果表格
# 定义变量：如果想比较不同的组，只需在这修改即可？能否直接写一个两两比较的循环进行自动处理？
Conditions<-"group"
sampleA <-"AMA"
sampleB <-"Young" 
#results函数：contrast <- c("conditions", sampleA, sampleB)  提取差异基因分析结果，包含log2 fold changes, p values和adjusted p values. 
#输出的log2FoldChange为log2(SampleA/SampleB)。
#如果直接调用results，只会提取出design matrix最后两个变量的 log2 fold changes and p values等结果。 
#results(dds)
#可以指定比较对象,并可以用mcols查看结果存储的元数据，了解列名的含义
contrastV <- c(Conditions, sampleA, sampleB)
res1 <- results(dds3, contrast=contrastV) 
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
res3$baseMean2 <- rowMeans(cbind(baseA, baseB))
head(res3)
# 校正后p-value为NA的复制为1
res3$padj[is.na(res3$padj)] <- 1

# 按pvalue排序, 把差异大的基因放前面
res3 <- res3[order(res3$pvalue),]
head(res3);tail(res3) ;dim(res3)

#差异基因的进一步筛选
# padj<0.05
res_de <- subset(res3, res3$padj<0.05, select=c('ID', 'log2FoldChange', 'padj'))
res_de<-res_de[order(res_de$padj, decreasing = F), ]
res_de_up <- subset(res_de, res_de$log2FoldChange>=0.585)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*0.585)

# 差异基因ID输出
#res_de_up_id = data.frame(ID=res_de_up$ID,type= Conditions)
res_de_up_id = data.frame(ID=res_de_dw$ID,type=paste(sampleA,"_higherThan_", sampleB, sep=""))
res_de_dw_id = data.frame(ID=res_de_dw$ID,type=paste(sampleA,"_lowerThan_", sampleB, sep=""))
de_id_whole = rbind(res_de_up_id, res_de_dw_id)
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_vs_",sampleB,".DEG_ID_padjust005FC1.5",'.xls', sep="")
write.table(as.data.frame(de_id_whole), file=file, sep="\t", quote=F, row.names=F, col.names=T)


res_de <- subset(res3, res3$padj<0.1, select=c('ID', 'log2FoldChange', 'padj'))
res_de<-res_de[order(res_de$padj, decreasing = F), ]
res_de_up <- subset(res_de, res_de$log2FoldChange>=0.585)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*0.585)

# 差异基因ID输出
#res_de_up_id = data.frame(ID=res_de_up$ID,type= Conditions)
res_de_up_id = data.frame(ID=res_de_up$ID,type=paste(sampleA,"_higherThan_", sampleB, sep=""))
res_de_dw_id = data.frame(ID=res_de_dw$ID,type=paste(sampleA,"_lowerThan_", sampleB, sep=""))
de_id_whole = rbind(res_de_up_id, res_de_dw_id)
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_vs_",sampleB,".DEG_ID_padjust01FC1.5",'.xls', sep="")
write.table(as.data.frame(de_id_whole), file=file, sep="\t", quote=F, row.names=F, col.names=T)

#====
#整体分析结果输出到文件
comp <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA, "_vs_", sampleB, sep="")
file_base <- paste(comp,"trans.Count_matrix.DESeq2", "results.xls", sep=".")
write.table(as.data.frame(res3), file=file_base, sep="\t", quote=F, row.names=F,col.names=T)
res3["RPS6KA5",]
#差异基因筛选
# padj<0.05
#res_de_up <- subset(res3, padj<0.05&log2FoldChange>=1, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
#1.5 FC padj<0.05
res_de_up <- subset(res3, padj<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue', 'padj'))
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_higherThan_",sampleB,"_padjust005FC1.5.trans.Count_matrix.DESeq2",'.xls', sep="")
write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F)
res_de_dw <- subset(res3, padj<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue', 'padj'))
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_lowerThan_",sampleB,"_padjust005FC1.5.trans.Count_matrix.xls.DESeq2",'.xls', sep="")
write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)
nrow(res_de_up )#249
nrow(res_de_dw )#136

res_de_up <- subset(res3, padj<0.1&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_higherThan_",sampleB,"_padjust01FC1.5.trans.Count_matrix.DESeq2",'.xls', sep="")
write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F)
res_de_dw <- subset(res3, padj<0.1&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_lowerThan_",sampleB,"_padjust01FC1.5.trans.Count_matrix.xls.DESeq2",'.xls', sep="")
write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)
nrow(res_de_up )#377
nrow(res_de_dw )#220

res_de_up <- subset(res3, pvalue<0.01&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_higherThan_",sampleB,"_pvalue001FC1.5.trans.Count_matrix.DESeq2",'.xls', sep="")
write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F)
res_de_dw <- subset(res3, pvalue<0.01&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_lowerThan_",sampleB,"_pvalue001FC1.5.trans.Count_matrix.xls.DESeq2",'.xls', sep="")
write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)
nrow(res_de_up )#50
nrow(res_de_dw )#46

res_de_up <- subset(res3, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_higherThan_",sampleB,"_pvalue005FC1.5.trans.Count_matrix.DESeq2",'.xls', sep="")
write.table(as.data.frame(res_de_up), file=file, sep="\t", quote=F, row.names=F)
res_de_dw <- subset(res3, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
file <- paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_lowerThan_",sampleB,"_pvalue005FC1.5.trans.Count_matrix.xls.DESeq2",'.xls', sep="")
write.table(as.data.frame(res_de_dw), file=file, sep="\t", quote=F, row.names=F)
nrow(res_de_up)# 207
nrow(res_de_dw)# 162

##以下选作=================
#最后设定阈值，筛选差异基因，导出全部数据。包括标准化后的count数)
##srict criticle
res3 <- res3[order(res3$padj),]
diff_gene <- subset(res3, padj < 0.05 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))
Count_data<-rlog_dds3_Mat#nor_Mat rlogMat,rlog_dds3_Mat
resdata <- merge(as.data.frame(diff_gene), as.data.frame(Count_data),by="row.names",sort=FALSE)
write.csv(resdata,paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_vs_",sampleB,"padj005_FC1.5.DEG_ID_rlog_nor_count",'.xls', sep=""), row.names = F)
dim(resdata)#3 21

res3 <- res3[order(res3$padj),]
diff_gene <- subset(res3, padj < 0.1 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))
Count_data<-rlog_dds3_Mat#nor_Mat rlogMat,rlog_dds3_Mat
resdata <- merge(as.data.frame(diff_gene), as.data.frame(Count_data),by="row.names",sort=FALSE)
write.csv(resdata,paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_vs_",sampleB,"padj01_FC1.5.DEG_ID_rlog_nor_count",'.xls', sep=""), row.names = F)
dim(resdata)#3 21

res3 <- res3[order(res3$padj),]
diff_gene <- subset(res3, pvalue<0.01 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))
Count_data<-rlog_dds3_Mat#nor_Mat rlogMat,rlog_dds3_Mat
resdata <- merge(as.data.frame(diff_gene), as.data.frame(Count_data),by="row.names",sort=FALSE)
write.csv(resdata,paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_vs_",sampleB,"pvalur0.01_FC1.5.DEG_ID_rlog_nor_count",'.xls', sep=""), row.names = F)
dim(resdata)#  96 21

res3 <- res3[order(res3$padj),]
diff_gene <- subset(res3, pvalue<0.05 & (log2FoldChange > 0.585 | log2FoldChange < -0.585))
Count_data<-rlog_dds3_Mat#nor_Mat rlogMat,rlog_dds3_Mat
resdata <- merge(as.data.frame(diff_gene), as.data.frame(Count_data),by="row.names",sort=FALSE)
write.csv(resdata,paste("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/",sampleA,"_vs_",sampleB,"pvalur0.05_FC1.5.DEG_ID_rlog_nor_count",'.xls', sep=""), row.names = F)
dim(resdata)#369  21


####基于结果的可视化

#Counts plot，看看特定基因的count数量。
#最显著基因
topGene <- rownames(res3)[which.min(res3$padj)] #3.840086e-13
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

#绘制火山图1--表达量差异
logCounts <- log2(res3$baseMeanmean+1)
logFC <- res1$log2FoldChange
FDR <- res1$padj
#png(filename=paste(file_base, "Volcano.png", sep="."))
plot(logFC, -1*log10(FDR), col=ifelse(FDR<=0.01, "red", "black"),xlab="logFC", ylab="-1*log1o(FDR)", main="Volcano plot", pch= 20)
#dev.off()

#正式火山图
#png(filename=paste(comp,"Volcano.png", sep=".")) 

res_de_up <- subset(res3, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
res_de_dw <- subset(res3, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
nrow(res_de_up )
nrow(res_de_dw )

filename<-paste(sampleA,"_versus_",sampleB, sep="")
#for pvalue
cols <- rep("#000000", nrow(res1))
cols[res1$log2FoldChange >= 0.585 & res1$pvalue < 0.05] <- "#CC0000"
cols[res1$log2FoldChange <= -0.585 & res1$pvalue < 0.05] <- "#0000CC"
x_lim<-max(res1$log2FoldChange,-res1$log2FoldChange)
y_lim<-max(-log10(na.omit(res1$pvalue)),log10(na.omit(res1$pvalue)))
plot(res1$log2FoldChange, -log10(res1$pvalue),main=filename, pch=16, cex=0.75,col=cols,las=1,xlim=c(-x_lim,x_lim), ylim=c(0,y_lim+2))
mtext(paste("up_gene: ",nrow(res_de_up)," & down_gene: ",nrow(res_de_dw),"  [pvalue <0.05 & FoldChange > 1.5]",sep="") , side=3 , line = 0.3 , outer = FALSE,cex=1,col="black") 
#添加辅助线
abline(h=-1*log10(0.05),lwd=3,lty=3,col="#4C5B61")
abline(v=log2(1.5) ,lwd=3,lty=3,col="#4C5B61")
abline(v=log2(2/3) ,lwd=3,lty=3,col="#4C5B61")

#for padj<0.05
res_de_up <- subset(res3, padj<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
res_de_dw <- subset(res3, padj<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
nrow(res_de_up )
nrow(res_de_dw )
cols <- rep("#000000", nrow(res1))
cols[res1$log2FoldChange >= 0.585 & res1$padj < 0.05] <- "#CC0000"
cols[res1$log2FoldChange <= -0.585 & res1$padj < 0.05] <- "#0000CC"
x_lim<-max(res1$log2FoldChange,-res1$log2FoldChange)
y_lim<-max(-log10(na.omit(res1$padj)),log10(na.omit(res1$padj)))
plot(res1$log2FoldChange, -log10(res1$padj),main=filename, pch=16, cex=0.75,col=cols,las=1,xlim=c(-x_lim,x_lim), ylim=c(0,y_lim+2))
mtext(paste("up_gene: ",nrow(res_de_up)," & down_gene: ",nrow(res_de_dw),"  [padj<0.05 & FoldChange > 1.5]",sep="") , side=3 , line = 0.3 , outer = FALSE,cex=1,col="black") 
#添加辅助线
abline(h=-1*log10(0.05),lwd=3,lty=3,col="#4C5B61")
abline(v=log2(3/2) ,lwd=3,lty=3,col="#4C5B61")
abline(v=log2(2/3) ,lwd=3,lty=3,col="#4C5B61")

#添加统计值
sig <- res1[which(res1$padj < 0.05),]
topGene <- rownames(res1)[which.min(res1$padj)]
with(sig[topGene, ], {
  points(log2FoldChange, -log10(padj), col="grey", cex=2, lwd=2)
  text(log2FoldChange, -log10(padj),paste(topGene,sep="") , pos=2, col="black")
  text(log2FoldChange, -log10(padj),"most significant", pos=1,  cex=0.7,col="black")
})

topGene <- rownames(res1)[which.max(res1$log2FoldChange)]
with(res1[topGene, ], {
  points(log2FoldChange, -log10(padj), col="grey", cex=2, lwd=2)
  text(log2FoldChange, -log10(padj), topGene, pos=2, col="green")
  text(log2FoldChange, -log10(padj),"most UP_change", pos=1,  cex=0.7,col="green")
})

topGene <- rownames(res1)[which.min(res1$log2FoldChange)]
with(res1[topGene, ], {
  points(log2FoldChange, -log10(padj), col="grey", cex=2, lwd=0.5)
  text(log2FoldChange, -log10(padj), topGene, pos=4, col="orange")
  text(log2FoldChange, -log10(padj),"most down_change", pos=1,  cex=0.7,col="orange")
})
#dev.off()


##基因聚类的热图展示。

#前20个样本件差异比较大，然后看他们在不同样本间的表达情况。
library("genefilter")
draw_count<-rlog_dds3_Mat#rlog_dd3_Mat;nor_Mat_var;nor_Mat
topVarGenes <- head(order(rowVars(draw_count), decreasing = TRUE), 20)
#rowVars(x, ...): Row variance and standard deviation of a numeric array
#order返回的这个向量刚好是xx向量中各数的原始的位置
mat  <- draw_count[topVarGenes, ]
#mat  <- mat - rowMeans(mat)
pheatmap(mat,annotation_col = anno1)
pheatmap(mat, cluster_rows=F, cluster_cols =T,scale="row", annotation_col=anno1)

#找样本间全部差异显著基因在不同样本间的表达情况--类似计算
draw_count<-rlog_dds3_Mat#rlog_dd3_Mat;nor_Mat_var;nor_Mat
resSig_RT <- subset(res3, padj< 0.05)
nrow(resSig_RT)
resSig_RT<-resSig_RT[order(resSig_RT$padj, decreasing = F), ]
namm<-rownames(draw_count)
mat  <- draw_count[which(namm %in% rownames(resSig_RT)), ]
#mat  <- mat - rowMeans(mat)
#pheatmap(mat, cluster_row=T, scale="none", annotation_col=anno1,main="scale_Genes_padjust<0.05")
pheatmap(mat, cluster_row=T, scale="row", annotation_col=anno1[2],main="scale_Genes_padjust<0.05")

resSig_RT <- subset(res3, padj< 0.1)
nrow(resSig_RT)
resSig_RT<-resSig_RT[order(resSig_RT$padj, decreasing = F), ]
namm<-rownames(draw_count)
mat  <- draw_count[which(namm %in% rownames(resSig_RT)), ]
#mat  <- mat - rowMeans(mat)
pheatmap(mat, cluster_row=T, scale="row", annotation_col=anno1[2],main="Genes_padjust<0.1")

#绘制各组特异高表达基因热图
res_de_up <- subset(res3, padj<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
res_de_dw <- subset(res3, padj<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
nrow(res_de_up )
nrow(res_de_dw )

# 差异基因ID输出
#res_de_up_id = data.frame(ID=res_de_up$ID,type= Conditions)
res_de_up_id = data.frame(ID=res_de_up$ID,type=paste(sampleA,"_higherThan_", sampleB, sep=""))
res_de_dw_id = data.frame(ID=res_de_dw$ID,type=paste(sampleA,"_lowerThan_", sampleB, sep=""))
de_id_whole = rbind(res_de_up_id, res_de_dw_id)
nrow(de_id_whole)
#by rlog of normal count
draw_count<-rlog_dds3_Mat#rlog_dd3_Mat;nor_Mat_var;nor_Mat
DEG_expr <- draw_count[which(rownames(draw_count) %in% unique(de_id_whole$ID)),]
pheatmap(DEG_expr, cluster_rows=T, cluster_cols =T,scale="row", annotation_col=anno1)
pheatmap(DEG_expr, cluster_rows=T, cluster_cols =T,scale="row", annotation_col=anno1,show_rownames = F)

#for loose
res_de_up <- subset(res3, pvalue <0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
res_de_dw <- subset(res3, pvalue <0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange','pvalue','padj'))
nrow(res_de_up )
nrow(res_de_dw )

# 差异基因ID输出
#res_de_up_id = data.frame(ID=res_de_up$ID,type= Conditions)
res_de_up_id = data.frame(ID=res_de_up$ID,type=paste(sampleA,"_higherThan_", sampleB, sep=""))
res_de_dw_id = data.frame(ID=res_de_dw$ID,type=paste(sampleA,"_lowerThan_", sampleB, sep=""))
de_id_whole = rbind(res_de_up_id, res_de_dw_id)
nrow(de_id_whole)
#by rlog of normal count
draw_count<-rlog_dds3_Mat#rlog_dd3_Mat;nor_Mat_var;nor_Mat
DEG_expr <- draw_count[which(rownames(draw_count) %in% unique(de_id_whole$ID)),]
pheatmap(DEG_expr, cluster_rows=T, cluster_cols =T,scale="row", annotation_col=anno1,show_rownames = F)


#上调与下调最明显的基因分别选取10个，然后看他们在不同样本间的表达情况--normalized_counts
draw_count<-rlog_dds3_Mat#rlog_dd3_Mat;nor_Mat_var;nor_Mat
res_de_up_top20_id <-  rownames(res_de_up[head(order(res_de_up$padj,decreasing = F),10),])
res_de_dw_top20_id <-  rownames(res_de_dw[head(order(res_de_dw$padj,decreasing = F),10),])
red_de_top20 <- c(res_de_up_top20_id, res_de_dw_top20_id)
red_de_top20 
red_de_top20_expr <- draw_count[rownames(draw_count) %in% red_de_top20,]
pheatmap(red_de_top20_expr, cluster_row=T, scale="row", annotation_col=anno1)


#热图修饰
#查看数据范围：
mat<-DEG_expr #pvalue <0.05&log2FoldChange>=0.585
range(mat)#1.638375 18.903036
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
pheatmap(mat, 
         scale = "row", 
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), 
         #         color = c(colorRampPalette(colors = c("blue","black"))(length(bk)/3),colorRampPalette(colors = c("black","red"))(length(bk)/3)), 
         legend_breaks=seq(-2,2,1), 
         breaks=bk,
         show_rownames = F,
         annotation_col = anno1,main="Genes_pvalue <0.05&log2FoldChange>=0.585"
)

#进一步美化：https://shengxin.ren/article/107

#富集分析
library(clusterProfiler)
library(org.Hs.eg.db)#org.Hs.eg.db 的数据类型以及使用简介：可用于数据类型转换
keytypes(org.Hs.eg.db)#查看org.Hs.eg.db数据对象里面包含着各大主流数据库的数据
###padj<0.1&FoldChange>=1.5
library(ReactomePA)
res_de_up<- as.data.frame(read.table("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Arrest_higherThan_Odinopoeia_padjust01FC1.5.trans.Count_matrix.DESeq2.xls",header = T,sep = "\t"))
res_de_dw<-as.data.frame(read.table("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Arrest_lowerThan_Odinopoeia_padjust01FC1.5.trans.Count_matrix.xls.DESeq2.xls",header = T,sep = "\t"))

#res_de_up <- subset(res3, padj<0.1&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
#res_de_dw <- subset(res3, padj<0.1&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
nrow(res_de_up )
nrow(res_de_dw )

#for up genes
gene<-res_de_up# res_de_dw,res_de
#导入数据，作为一个整合数据，里我们要用到的只是entrez ID列和logFC列：
#由于clusterProfiler富集分析推荐的输入文件是Entrez ID，因此这里提取的是Entrez ID

#translate into other types ID
gene.df <- bitr(as.character(gene$ID), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
#for hyper
target<-unique(gene.df$ENTREZID)
str(target)
gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
head(summary(BP))
CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))
write.table(as.data.frame(CC@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust01FC1.5_CC_GO.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust01FC1.5_BP_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust01FC1.5_MF_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(CC_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust01FC1.5_BP_GO_simp.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust01FC1.5_MF_GO_simp.txt", row.names=T, col.names=T) 
head(as.data.frame(BP_simp))
write.table((as.data.frame(kk@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust01FC1.5_KEGG.txt",row.names=T, col.names=T) 
write.table((as.data.frame(rpea@result)), file="//media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust01FC1.5_RPEA.txt",row.names=T, col.names=T) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust01FC1.5_GO_RPEA_result.rds")

##for up genes
gene<-res_de_dw
#translate into other types ID
gene.df <- bitr(as.character(gene$ID), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
#for hyper
target<-unique(gene.df$ENTREZID)
str(target)


gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
head(as.data.frame(BP@result));head(as.data.frame(MF@result));head(as.data.frame(CC@result))

CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result));head(as.data.frame(MF_simp@result));head(as.data.frame(CC_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.05,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))
write.table(as.data.frame(CC@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust01FC1.5_CC_GO.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust01FC1.5_BP_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment//Verification_AMA_versus_Young_Down_gene_padjust01FC1.5_MF_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(CC_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust01FC1.5_CC_GO_simp.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust01FC1.5_BP_GO_simp.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust01FC1.5_MF_GO_simp.txt",row.names=T, col.names=T) 
write.table((as.data.frame(kk@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust01FC1.5_KEGG.txt", row.names=T, col.names=T) 
write.table((as.data.frame(rpea@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust01FC1.5_RPEA.txt",row.names=T, col.names=T) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust01FC1.5_GO_RPEA_result.rds")

###padj<0.05&FoldChange>=1.5
res_de_up<- as.data.frame(read.table("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Arrest_higherThan_Odinopoeia_padjust005FC1.5.trans.Count_matrix.DESeq2.xls",header = T,sep = "\t"))
res_de_dw<-as.data.frame(read.table("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Arrest_lowerThan_Odinopoeia_padjust005FC1.5.trans.Count_matrix.xls.DESeq2.xls",header = T,sep = "\t"))
#res_de_up <- subset(res3, padj<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
#res_de_dw <- subset(res3, padj<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
nrow(res_de_up )
nrow(res_de_dw )

#for up genes
gene<-res_de_up# res_de_dw,res_de
#导入数据，作为一个整合数据，里我们要用到的只是entrez ID列和logFC列：
#由于clusterProfiler富集分析推荐的输入文件是Entrez ID，因此这里提取的是Entrez ID

#translate into other types ID
gene.df <- bitr(as.character(gene$ID), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
#for hyper
target<-unique(gene.df$ENTREZID)
str(target)
gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
head(summary(BP))
CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))
write.table(as.data.frame(CC@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust0.05FC1.5_CC_GO.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust0.05FC1.5_BP_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust0.05FC1.5_MF_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(CC_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust0.05FC1.5_BP_GO_simp.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust0.05FC1.5_MF_GO_simp.txt", row.names=T, col.names=T) 
head(as.data.frame(BP_simp))
write.table((as.data.frame(kk@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust0.05FC1.5_KEGG.txt",row.names=T, col.names=T) 
write.table((as.data.frame(rpea@result)), file="//media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust0.05FC1.5_RPEA.txt",row.names=T, col.names=T) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust0.05FC1.5_GO_RPEA_result.rds")

##for up genes
gene<-res_de_dw
#translate into other types ID
gene.df <- bitr(as.character(gene$ID), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
#for hyper
target<-unique(gene.df$ENTREZID)
str(target)


gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
head(as.data.frame(BP@result));head(as.data.frame(MF@result));head(as.data.frame(CC@result))

CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result));head(as.data.frame(MF_simp@result));head(as.data.frame(CC_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.05,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))
write.table(as.data.frame(CC@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_CC_GO.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_BP_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment//Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_MF_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(CC_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_CC_GO_simp.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_BP_GO_simp.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_MF_GO_simp.txt",row.names=T, col.names=T) 
write.table((as.data.frame(kk@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_KEGG.txt", row.names=T, col.names=T) 
write.table((as.data.frame(rpea@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_RPEA.txt",row.names=T, col.names=T) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_GO_RPEA_result.rds")

###pvalue<0.01&FoldChange>=1.5
res_de_up<- as.data.frame(read.table("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Arrest_higherThan_Odinopoeia_pvalue001FC1.5.trans.Count_matrix.DESeq2.xls",header = T,sep = "\t"))
res_de_dw<-as.data.frame(read.table("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Arrest_lowerThan_Odinopoeia_pvalue001FC1.5.trans.Count_matrix.xls.DESeq2.xls",header = T,sep = "\t"))
#res_de_up <- subset(res3, pvalue<0.01&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
#res_de_dw <- subset(res3, pvalue<0.01&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
nrow(res_de_up)
nrow(res_de_dw)

#for up genes
gene<-res_de_up   ## res_de_dw,res_de_up
#导入数据，作为一个整合数据，里我们要用到的只是entrez ID列和logFC列：
#由于clusterProfiler富集分析推荐的输入文件是Entrez ID，因此这里提取的是Entrez ID

#translate into other types ID
gene.df <- bitr(as.character(gene$ID), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
#for hyper
target<-unique(gene.df$ENTREZID)
str(target)
gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
head(summary(BP))
CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))
write.table(as.data.frame(CC@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.01FC1.5_CC_GO.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.01FC1.5_BP_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.01FC1.5_MF_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(CC_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.01FC1.5_BP_GO_simp.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.01FC1.5_MF_GO_simp.txt", row.names=T, col.names=T) 
head(as.data.frame(BP_simp))
write.table((as.data.frame(kk@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.01FC1.5_KEGG.txt",row.names=T, col.names=T) 
write.table((as.data.frame(rpea@result)), file="//media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.01FC1.5_RPEA.txt",row.names=T, col.names=T) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.01FC1.5_GO_RPEA_result.rds")

##for up genes
gene<-res_de_dw
#translate into other types ID
gene.df <- bitr(as.character(gene$ID), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
#for hyper
target<-unique(gene.df$ENTREZID)
str(target)


gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
head(as.data.frame(BP@result));head(as.data.frame(MF@result));head(as.data.frame(CC@result))

CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result));head(as.data.frame(MF_simp@result));head(as.data.frame(CC_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.05,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))
write.table(as.data.frame(CC@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.01FC1.5_CC_GO.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.01FC1.5_BP_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment//Verification_AMA_versus_Young_Down_gene_pvalue0.01FC1.5_MF_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(CC_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.01FC1.5_CC_GO_simp.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.01FC1.5_BP_GO_simp.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.01FC1.5_MF_GO_simp.txt",row.names=T, col.names=T) 
write.table((as.data.frame(kk@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.01FC1.5_KEGG.txt", row.names=T, col.names=T) 
write.table((as.data.frame(rpea@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.01FC1.5_RPEA.txt",row.names=T, col.names=T) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.01FC1.5_GO_RPEA_result.rds")

###pvalue<0.05&FoldChange>=1.5
res_de_up<- as.data.frame(read.table("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Arrest_higherThan_Odinopoeia_pvalue005FC1.5.trans.Count_matrix.DESeq2.xls",header = T,sep = "\t"))
res_de_dw<-as.data.frame(read.table("/mnt/data/chenwei/huahua/3.mapped_data/valification/trans/4.count/Arrest_lowerThan_Odinopoeia_pvalue005FC1.5.trans.Count_matrix.xls.DESeq2.xls",header = T,sep = "\t"))
#res_de_up <- subset(res3, pvalue<0.05&log2FoldChange>=0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
#res_de_dw <- subset(res3, pvalue<0.05&log2FoldChange<=(-1)*0.585, select=c('ID', colnames(baseMeanA),colnames(baseMeanB), 'log2FoldChange', 'padj'))
nrow(res_de_up );nrow(res_de_dw )

#for up genes
gene<-res_de_up# res_de_dw,res_de
#导入数据，作为一个整合数据，里我们要用到的只是entrez ID列和logFC列：
#由于clusterProfiler富集分析推荐的输入文件是Entrez ID，因此这里提取的是Entrez ID

#translate into other types ID
gene.df <- bitr(as.character(gene$ID), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
#for hyper
target<-unique(gene.df$ENTREZID)
str(target)
gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
head(summary(BP))
CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))
write.table(as.data.frame(CC@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.05FC1.5_CC_GO.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.05FC1.5_BP_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.05FC1.5_MF_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(CC_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.05FC1.5_BP_GO_simp.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.05FC1.5_MF_GO_simp.txt", row.names=T, col.names=T) 
head(as.data.frame(BP_simp))
write.table((as.data.frame(kk@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.05FC1.5_KEGG.txt",row.names=T, col.names=T) 
write.table((as.data.frame(rpea@result)), file="//media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.05FC1.5_RPEA.txt",row.names=T, col.names=T) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_pvalue0.05FC1.5_GO_RPEA_result.rds")

##for up genes
gene<-res_de_dw
#translate into other types ID
gene.df <- bitr(as.character(gene$ID), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
#for hyper
target<-unique(gene.df$ENTREZID)
str(target)


gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
head(as.data.frame(BP@result));head(as.data.frame(MF@result));head(as.data.frame(CC@result))

CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result));head(as.data.frame(MF_simp@result));head(as.data.frame(CC_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.05,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))
write.table(as.data.frame(CC@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.05FC1.5_CC_GO.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.05FC1.5_BP_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment//Verification_AMA_versus_Young_Down_gene_pvalue0.05FC1.5_MF_GO.txt",row.names=T, col.names=T) 
write.table(as.data.frame(CC_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.05FC1.5_CC_GO_simp.txt", row.names=T, col.names=T) 
write.table(as.data.frame(BP_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.05FC1.5_BP_GO_simp.txt",row.names=T, col.names=T) 
write.table(as.data.frame(MF_simp@result), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.05FC1.5_MF_GO_simp.txt",row.names=T, col.names=T) 
write.table((as.data.frame(kk@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.05FC1.5_KEGG.txt", row.names=T, col.names=T) 
write.table((as.data.frame(rpea@result)), file="/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.05FC1.5_RPEA.txt",row.names=T, col.names=T) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_pvalue0.05FC1.5_GO_RPEA_result.rds")



##富集绘图
#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#结果展示 https://www.jianshu.com/p/e133ab3169fa 
gene_erichment_results_Up<-readRDS(file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Up_gene_padjust0.05FC1.5_GO_RPEA_result.rds")
gene_erichment_results_Down<-readRDS(file = "/media/data4/qinmen_BR/qinmen_RNA/6.enrichment/Verification_AMA_versus_Young_Down_gene_padjust0.05FC1.5_GO_RPEA_result.rds")


#for genes in 3kb
#For hyper
gene_erichment_results<-gene_erichment_results_Up #gene_erichment_results_hyper_unselect
CC_plot<-gene_erichment_results[[4]]
BP_plot<-gene_erichment_results[[5]]
MF_plot<-gene_erichment_results[[6]]

#结果展示

merge_three<-as.data.frame(rbind(head(BP_plot@result,n=10),head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
merge_three$qvalue_log10<-c(-log(merge_three$qvalue,10))
merge_three$pvalue_log10<-c(-log(merge_three$pvalue,10))

merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))
pdf("Verification_AMA_versus_Young_UP_genes_enrichment_BP_CC_MF.pdf")
ggbarplot(merge_three, x="Description", y="qvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
                    rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=60)+
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
BP_plot<-gene_erichment_results[[5]]
MF_plot<-gene_erichment_results[[6]]
CC_plot<-gene_erichment_results[[4]]

#结果展示
merge_three<-as.data.frame(rbind(head(BP_plot@result,n=10),head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
dim(merge_three)
#merge_three$qlogtran<-c(-log(merge_three$qvalue,10))
#merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))
#merge_two<-as.data.frame(rbind(head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
merge_three$qvalue_log10<-c(-log(merge_three$qvalue,10))
merge_three$pvalue_log10<-c(-log(merge_three$pvalue,10))
merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))

pdf("Verification_AMA_versus_Young_Down_gene_enrichment_CC_MF_BP.pdf")
ggbarplot(merge_three, x="Description", y="qvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
                    rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=60)+
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
pdf("Arrest_verse_Odinopoeia-Up_Down_gene_KEGG.pdf")

p <- ggplot(PP_plot,aes(qvalue_log10 ,reorder(Description,qvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(p.adjust)",y="Description",title="KEGG enrichment")+
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
str(KEGG_plot)
head(KEGG_plot)

PP_plot<-KEGG_plot
head(PP_plot)
ggplot(PP_plot,aes(pvalue_log10 ,reorder(Description,pvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(pvalue)",y="Description",title="KEGG enrichment")+
  scale_size_continuous(range=c(5,14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=55))+
  geom_text(aes(label=sprintf("%.3f",pvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")
dev.off()

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
barplot(GO_plot_data, showCategory=15,title="EnrichmentGO_ALL")#条状图，按p从小到大排的#（这个图...好丑)

#绘制富集网格图 Enrichment Map 可以将富集条目和重叠的基因集整合为一个网络图，相互重叠的基因集则趋向于成簇，从而易于分辨功能模型。
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

## Reactome pathway enrichment analysis
require(ReactomePA)
rpea <- enrichPathway(gene=gene$ID,pvalueCutoff=0.05, readable=T,organism = "human")
head(as.data.frame(rpea))

file <- paste("H:/ART_Epigenome_Papers/ART项目数据/wholeBlood_RNAseq_20180916/DEseq2_20190809_by_CW/",sampleA,"_vs_",sampleB,sep="")
write.table(as.data.frame(CC_simp@result), file=paste(file,"hyper_RPEA.txt",sep="."),quote=F, row.names=F, col.names=T) 


#visulization for pathway
pathway_data<-kk
barplot(pathway_data,showCategory=10)
dotplot(pathway_data,showCategory=10)
emapplot(pathway_data)  #enrichment map
cnetplot(pathway_data, categorySize="pvalue", foldChange=gene) 
heatplot(pathway_data, gene$log2FoldChange)

#browseKEGG
#函数 browseKEGG 可以帮你打开浏览器，嗯。
browseKEGG(pathway_data, 'hsa04110')

###多个基因集时（eg，WGCNA得到的不同module的基因集），用compareCluster。
#参考：https://www.jianshu.com/p/c01b4cc1b98a




##GSEA分析 Gene Set Enrichment Analysis(GSEA) 参考：https://www.jianshu.com/p/feaefcbdf986
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
