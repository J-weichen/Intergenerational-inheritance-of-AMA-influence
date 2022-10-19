rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')

#library(methylKit)
#载入各种包
library(genomation)
library(reshape2)
library(RColorBrewer)
library(grid)
library(scales)
library(ggsci)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr) 
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#绘制目标基因的目标区域甲基化水平柱状图
#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#def 15%
compare_name<-"mother_AMA_vs_Young"
split_region<- "200bp"
myDiff.hyper<- as.data.frame(read.csv(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_hyper.csv"),header = T,sep = ",",check.names = F))
myDiff.hypo<- as.data.frame(read.csv(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_hypo.csv"),header = T,sep = ",",check.names = F))
myDiff.all<- as.data.frame(read.csv(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_all.csv"),header = T,sep = ",",check.names = F))
head(myDiff.hyper)
#先是读取基因注释信息，用这个注释信息对差异甲基化区域进行注释
gene.obj=readTranscriptFeatures("/mnt/data/chenwei/qinmen_BR/00.ref_data/Regions_bed/RefSeq_Genes.bed")
diffAnn_hyper=annotateWithGeneParts(as(myDiff.hyper,"GRanges"),gene.obj)
diffAnn_hypo=annotateWithGeneParts(as(myDiff.hypo,"GRanges"),gene.obj)
diffAnn_all=annotateWithGeneParts(as(myDiff.all,"GRanges"),gene.obj)
head(getAssociationWithTSS(diffAnn_all))
#将差异甲基化最近的基因名输出到文件里
#def15%
write.csv(getAssociationWithTSS(diffAnn_hyper),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_gene_annotation_",split_region,"_DMRs_myDiff15q005_hyper.csv"))
write.csv(getAssociationWithTSS(diffAnn_hypo),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_gene_annotation_",split_region,"_DMRs_myDiff15q005_hypo.csv"))
write.csv(getAssociationWithTSS(diffAnn_all),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_gene_annotation_",split_region,"_DMRs_myDiff15q005_all.csv"))
#注意：同样的函数methylKit会报错,因为其已经弃用
genomation::getTargetAnnotationStats(diffAnn_hyper,percentage=TRUE,precedence=TRUE)
genomation::getTargetAnnotationStats(diffAnn_hypo,percentage=TRUE,precedence=TRUE)
genomation::getTargetAnnotationStats(diffAnn_all,percentage=TRUE,precedence=TRUE)
#"kids_AMA_vs_Young"
#hyper
#promoter       exon     intron intergenic 
#   12.29       4.66      47.03      36.02 
#hypo
#promoter       exon     intron intergenic 
# 8.57       6.79      43.93      40.71 
#all
#promoter       exon     intron intergenic 
#10.27       5.81      45.35      38.57 

pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_TSS_",split_region,"_DMRs_myDiff15q005_all_gene_element_Annotation.pdf"))
genomation::plotTargetAnnotation(diffAnn_all,precedence=TRUE, main=paste0(compare_name,"_TSS_myDiff15q005_all_DMR annotation"))
dev.off()
pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_TSS_",split_region,"_DMRs_myDiff15q005_hyper_gene_element_Annotation.pdf"))
genomation::plotTargetAnnotation(diffAnn_hyper,precedence=TRUE, main=paste0(compare_name,"_TSS_myDiff15q005_hyper_DMR annotation"))
dev.off()
pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_TSS_Young-AMA_",split_region,"_DMRs_myDiff15q005_hypo_gene_element_Annotation.pdf"))
genomation::plotTargetAnnotation(diffAnn_hypo,precedence=TRUE, main=paste0(compare_name,"_TSS_myDiff15q005_hypo_DMR annotation"))
dev.off()

#读取CpG island的注释信息，用这些注释信息来注释我们差异甲基化的区域
cpg.obj=readFeatureFlank("/mnt/data/chenwei/qinmen_BR/00.ref_data/Regions_bed/7-CpGislands.bed",feature.flank.name=c("CpGi","shores"))

diffCpGann_hyper=annotateWithFeatureFlank(as(myDiff.hyper,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
diffCpGann_hypo=annotateWithFeatureFlank(as(myDiff.hypo,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
diffCpGann_all=annotateWithFeatureFlank(as(myDiff.all,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")

pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_hyper_CpGisland.pdf"))
genomation::plotTargetAnnotation(diffCpGann_hyper,col=c("red","orange","gray"), main=paste0(compare_name,"_myDiff15q005_hyper_DMR annotation"))
dev.off()

pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_hypo_CpGisland.pdf"))
genomation::plotTargetAnnotation(diffCpGann_hypo,col=c("red","orange","gray"), main=paste0(compare_name,"_myDiff15q005_hypo_DMR annotation"))
dev.off()
pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_all_CpGisland.pdf"))
genomation::plotTargetAnnotation(diffCpGann_all,col=c("red","orange","gray"), main=paste0(compare_name,"_myDiff15q005_all_DMR annotation"))
dev.off()

# 根据之前的注释信息，可以得到起始子区域/CpG island区域的位置，
# 然后下面的方法可以总结这些区域的甲基化信息                      
# promoters=regionCounts(myobj,gene.obj$promoters)
# head(promoters[[1]])
genomation::getFeatsWithTargetsStats(diffCpGann_all,percentage=TRUE)
genomation::getTargetAnnotationStats(diffCpGann_all,percentage=TRUE,precedence=TRUE)

write.csv(gene.obj$promoters,paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_promoters.csv"))
write.csv(gene.obj$exons,paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_exons.csv"))
write.csv(gene.obj$introns,paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_introns.csv"))
write.csv(gene.obj$TSSes,paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_TSSs.csv"))

#进一步注释
dim(myDiff.all);dim(myDiff.hyper);dim(myDiff.hypo);range(abs(myDiff.all$meth.diff))
#[1] 516  32
#[1] 236  32
#[1] 280  32
#[1] 15.00851 37.94176

head(myDiff.all);head(myDiff.hyper);head(myDiff.hypo)
DMRs_hyper<-readPeakFile(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_hyper.bed"))
DMRs_hypo<-readPeakFile(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_hypo.bed"))
DMRs_all <- readPeakFile(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_all.bed"))

#3.与TSS距离以及基因注释 ###"Young_versus_AMA"在注释这一步有问题
#ref: https://guangchuangyu.github.io/cn/2017/03/peak注释/
##peaks的注释是用的annotatePeak函数，可以单独对每个peaks文件进行注释，
##也可以将多个peaks制作成一个list，进行比较分析和可视化。
# 制作多个样本比较的list
peaks <- list(DMRs_all=DMRs_all,DMRs_hyper=DMRs_hyper,DMRs_hypo=DMRs_hypo)
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
#查看具体信息
as.GRanges(peakAnnoList$DMRs_all)%>%head(3)
head(as.data.frame(peakAnnoList$DMRs_all))
#两种注释有时候还不够，想看peak上下游某个范围内（比如说-5k到5k的距离）都有什么基因，annotatePeak也可以做到。
#annotatePeak传入annoDb参数,可进行基因ID转换（Entrez，ENSEMBL，SYMBOL，GENENAME）
#peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
#head(peakAnnoList)
#第4步：可视化
#提供了多种可视化方法，如plotAnnoBar(),vennpie(),plotAnnoPie(),plotDistToTSS()等，
#下面展示了两个样本在基因组特征区域的分布以及转录因子在TSS区域的结合
pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_distribution_gene.pdf"))
plotAnnoBar(peakAnnoList,title=paste0(compare_name,": Distribution of_",split_region,"_DMRs_myDiff15q005 \n on genes related element"))
dev.off()
pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMRs_myDiff15q005_distribution_TSS.pdf"))
plotDistToTSS(peakAnnoList,title=paste0(compare_name,":Distribution of_",split_region,"_DMRs_myDiff15q005 \n relative to TSS"))
dev.off()

#注释输出
annotation_DMR<-as.data.frame(peakAnnoList$DMRs_all)
annotation_DMR<-annotation_DMR[,c("seqnames","start","end","width","strand","annotation","distanceToTSS","geneStrand","transcriptId","geneId","ENSEMBL", "SYMBOL", "GENENAME")]
#write.table(annotation_DMR, file="/media/data2/lucunlin/In_Vivo-In_vitro_annotation_for_DMRs_all.xls",sep=" ",quote=F, row.names=F, col.names=T) 
str(annotation_DMR)
head(annotation_DMR);dim(annotation_DMR)
write.table(annotation_DMR, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_all_myDiff15q005.txt"),quote=F,row.names=F,col.names=T,sep = "\t")
write.csv(as.data.frame(annotation_DMR),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_all_myDiff15q005.csv"))

annotation_DMR_hyper<-as.data.frame(peakAnnoList$DMRs_hyper)
annotation_DMR_hyper<-annotation_DMR_hyper[,c("seqnames","start","end","width","strand","annotation","distanceToTSS","geneStrand","transcriptId","geneId","ENSEMBL", "SYMBOL", "GENENAME")]
head(annotation_DMR_hyper);dim(annotation_DMR_hyper)
write.table(annotation_DMR_hyper, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hyper_myDiff15q005.txt"),quote=F,row.names=F,col.names=T,sep = "\t") 
write.csv(as.data.frame(annotation_DMR_hyper),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hyper_myDiff15q005.csv"))

annotation_DMR_hypo<-as.data.frame(peakAnnoList$DMRs_hypo)
annotation_DMR_hypo<-annotation_DMR_hypo[,c("seqnames","start","end","width","strand","annotation","distanceToTSS","geneStrand","transcriptId","geneId","ENSEMBL", "SYMBOL", "GENENAME")]
head(annotation_DMR_hypo);dim(annotation_DMR_hypo)
write.table(annotation_DMR_hypo, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hypo_myDiff15q005.txt"),quote=F,row.names=F,col.names=T,sep = "\t") 
write.csv(as.data.frame(annotation_DMR_hypo),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hypo_myDiff15q005.csv"))

#select DMR within 3KB from TSS
#annotation_DMR2<-read.table(file="/media/data2/lucunlin/In_Vivo-In_vitro_annotation_for_DMRs_all_myDiff15q005.txt",header=T,sep="") 
annotation_DMR_all<-as.data.frame(peakAnnoList$DMRs_all)
range(annotation_DMR_all$distanceToTSS)#-324023  875111
annotation_DMR_all_3kb<-annotation_DMR_all[which(abs(annotation_DMR_all$distanceToTSS)<=3000),]
annotation_DMR_all_3kb<-annotation_DMR_all_3kb[,c("seqnames","start","end","width","strand","annotation","distanceToTSS","geneStrand","transcriptId","geneId","ENSEMBL", "SYMBOL", "GENENAME")]
head(annotation_DMR_all_3kb)
write.table(annotation_DMR_all_3kb, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_all_myDiff15q005_within_TSS_3kb.txt",quote=F,row.names=F,col.names=T) )
write.csv(as.data.frame(annotation_DMR_all_3kb),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_all_myDiff15q005_within_TSS_3kb.csv"))

annotation_DMR_hyper<-as.data.frame(peakAnnoList$DMRs_hyper)
range(annotation_DMR_hyper$distanceToTSS)#-553584  584923
annotation_DMR_hyper_3kb<-annotation_DMR_hyper[which(abs(annotation_DMR_hyper$distanceToTSS)<=3000),]
annotation_DMR_hyper_3kb<-annotation_DMR_hyper_3kb[,c("seqnames","start","end","width","strand","annotation","distanceToTSS","geneStrand","transcriptId","geneId","ENSEMBL", "SYMBOL", "GENENAME")]
head(annotation_DMR_hyper_3kb)
write.table(annotation_DMR_hyper_3kb, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hyper_myDiff15q005_within_TSS_3kb.txt",quote=F,row.names=F,col.names=T))
write.csv(as.data.frame(annotation_DMR_hyper_3kb),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hyper_myDiff15q005_within_TSS_3kb.csv"))

annotation_DMR_hypo<-as.data.frame(peakAnnoList$DMRs_hypo)
range(annotation_DMR_hypo$distanceToTSS)#-324023  366793
annotation_DMR_hypo_3kb<-annotation_DMR_hypo[which(abs(annotation_DMR_hypo$distanceToTSS)<=3000),]
annotation_DMR_hypo_3kb<-annotation_DMR_hypo_3kb[,c("seqnames","start","end","width","strand","annotation","distanceToTSS","geneStrand","transcriptId","geneId","ENSEMBL", "SYMBOL", "GENENAME")]
head(annotation_DMR_hypo_3kb)
write.table(annotation_DMR_hypo_3kb, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hypo_myDiff15q005_within_TSS_3kb.txt",quote=F,row.names=F,col.names=T) )
write.csv(as.data.frame(annotation_DMR_hypo_3kb),paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hypo_myDiff15q005_within_TSS_3kb.csv"))

#注释数据与位点甲基化水平合并
myDiff_data<- read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv", header=T,sep=",",check.names = F)
head(myDiff_data)
myDiff_data$bin_region <- unlist(lapply(strsplit(as.character(myDiff_data$Row.names), "[.]"), function(x) paste(x[1],x[2],x[3],sep="_")))
rownames(myDiff_data)<-myDiff_data$bin_region
head(myDiff_data);dim(myDiff_data)

data_region_merge<-myDiff_data[2: ncol(myDiff_data)]
head(data_region_merge);str(data_region_merge);dim(data_region_merge)

data_region_merge$meth_df2<-(data_region_merge$AMA_Mother_mean -data_region_merge$YOUNG_Mother_mean)
data_region_merge$log10_FC<-log10(data_region_merge$AMA_Mother_mean/data_region_merge$YOUNG_Mother_mean)

head(data_region_merge,n=3);dim(data_region_merge)
length(which(abs(data_region_merge$meth.diff)>=20))#109

#与最近基因以及基因组位置注释合并
annotation_DMR<-read.delim(file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_all_myDiff15q005.txt"),header =T, sep="\t") 
head(annotation_DMR,n=3);dim(annotation_DMR)

annotation_DMR$start<-annotation_DMR$start-1
rownames(annotation_DMR)<-paste0(annotation_DMR$seqnames,"_",annotation_DMR$start,"_",annotation_DMR$end)
length(rownames(annotation_DMR)[as.character(rownames(annotation_DMR)) %in% rownames(data_region_merge)])

data_region_merge2<-data_region_merge[rownames(annotation_DMR),]
dim(data_region_merge2)

head(annotation_DMR);head(data_region_merge2)
data_region_merge_DMR_all<-merge(annotation_DMR,data_region_merge2,by=0)
head(data_region_merge_DMR_all,n=3);dim(data_region_merge_DMR_all)#516  43

write.table(data_region_merge_DMR_all,paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),col.names=T,row.names=T,sep="\t")
#write.table(annotation_DMR,paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_data_region_merge_DMR_defmyDiff15q005_all_annotation_modified.txt"),col.names=T,row.names=T,sep="\t")