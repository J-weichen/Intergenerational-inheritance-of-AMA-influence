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
      
AMA_Mother_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/valification/meth/bismark_cov_file/AMA_mom")
YOUNG_Mother_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/valification/meth/bismark_cov_file/Young_mom")
AMA_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/valification/meth/bismark_cov_file/AMA_Kid")
YOUNG_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/valification/meth/bismark_cov_file/Young_Kid")

AMA_Mother_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/valification/meth/bismark_cov_file/AMA_mom/",AMA_Mother_files))
YOUNG_Mother_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/valification/meth/bismark_cov_file/Young_mom/",YOUNG_Mother_files))
AMA_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/valification/meth/bismark_cov_file/AMA_Kid/",AMA_Kids_files))
YOUNG_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/valification/meth/bismark_cov_file/Young_Kid/",YOUNG_Kids_files))

AMA_Mother_names<- lapply(strsplit(AMA_Mother_files,"_"), function(x) x[1])
YOUNG_Mother_names<- lapply(strsplit(YOUNG_Mother_files,"_"), function(x) x[1])
AMA_Kids_names<- lapply(strsplit(AMA_Kids_files,"_"), function(x) x[1])
YOUNG_Kids_names<- lapply(strsplit(YOUNG_Kids_files,"_"), function(x) x[1])
#covarage 6
myobj=methRead(c(AMA_Kids_files2,YOUNG_Kids_files2,AMA_Mother_files2,YOUNG_Mother_files2),
               sample.id=c(AMA_Kids_names,YOUNG_Kids_names,AMA_Mother_names,YOUNG_Mother_names),
               assembly="hg38",pipeline="bismarkCoverage",#default:'amp'
               treatment=c(rep(0,length(AMA_Kids_names)),rep(1,length(YOUNG_Kids_names)),
                           rep(2,length(AMA_Mother_names)),rep(3,length(YOUNG_Mother_names))),
               context="CpG", resolution = "base", mincov = 6)
dim(as.data.frame(myobj[1]))
saveRDS(myobj, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_myobj_6X_CpG.rds")

#进一步过滤 设置最低覆盖度
myobj<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_myobj_6X_CpG.rds")
depth<-6
filtered.myobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)

#基因组划分区域进行计算
tiles_100=tileMethylCounts(filtered.myobj,win.size=100,step.size=100,cov.bases =1)
meth_100bp=methylKit::unite(tiles_100, destrand=FALSE,min.per.group=0L,mc.cores = 40)

data_raw_CpG<-as.data.frame(percMethylation(meth_100bp,rowids=T))
meth_100bp_2<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY",rownames(data_raw_CpG),invert=TRUE),]
head(meth_100bp_2);tail(meth_100bp_2);dim(meth_100bp_2)#2636827      20

colMeans(meth_100bp_2,na.rm=TRUE)
#   A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M      A4M      A5M 
#  71.06666 76.30424 77.63867 74.99075 76.02178 76.23057 76.01573 76.12958 76.95723 76.59459 75.01546 75.06294 74.80134 75.17037 
#  A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
# 73.53968 75.85593 75.67770 75.91637 74.89375 75.27171 
colMeans(data_raw_CpG,na.rm=TRUE)
#   A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M      A4M      A5M 
# 71.07068 76.30329 77.63787 74.99008 76.03187 76.23718 76.01526 76.12924 76.95645 76.59389 75.01499 75.06229 74.80090 75.16988 
#   A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
#73.53913 75.85512 75.67733 75.91584 74.89308 75.27116 
#基因组划分区域进行计算
tiles_200=tileMethylCounts(filtered.myobj,win.size=200,step.size=200,cov.bases =1)
meth_200bp=methylKit::unite(tiles_200, destrand=FALSE,min.per.group=0L,mc.cores = 40)
data_raw_CpG<-as.data.frame(percMethylation(meth_200bp,rowids=T))
meth_200bp<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY",rownames(data_raw_CpG),invert=TRUE),]
head(meth_200bp);tail(meth_200bp);dim(meth_200bp)
colMeans(meth_200bp,na.rm=TRUE)
# A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M      A4M      A5M 
# 72.53323 77.85149 79.20159 76.55299 77.95672 77.84361 77.48598 78.02655 78.84366 78.14057 76.76078 76.74452 76.27228 76.91838 
# A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
# 75.17069 77.41997 77.43969 77.65633 76.36323 76.96328 
colMeans(data_raw_CpG,na.rm=TRUE)
# A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M      A4M      A5M 
# 72.53603 77.85022 79.20052 76.55219 77.96536 77.84885 77.48534 78.02608 78.84272 78.13958 76.76011 76.74368 76.27154 76.91772 
# A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
# 75.16986 77.41890 77.43904 77.65557 76.36223 76.96265 
#基因组划分区域进行计算
tiles_300=tileMethylCounts(filtered.myobj,win.size=300,step.size=300,cov.bases =1)
meth_300bp=methylKit::unite(tiles_300, destrand=FALSE,min.per.group=0L,mc.cores = 40)
data_raw_CpG<-as.data.frame(percMethylation(meth_300bp,rowids=T))
meth_300bp<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY",rownames(data_raw_CpG),invert=TRUE),]
head(meth_300bp);tail(meth_300bp);dim(meth_300bp)

colMeans(meth_300bp,na.rm=TRUE)
# A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M      A4M      A5M 
# 73.32173 78.69180 80.06222 77.40422 78.97622 78.71771 78.28754 79.01631 79.84646 78.97404 77.70228 77.66476 77.07511 77.87553 
# A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
# 76.03910 78.27113 78.38253 78.59199 77.17272 77.88757 
colMeans(data_raw_CpG,na.rm=TRUE)
#  A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M      A4M      A5M 
# 73.32369 78.69042 80.06119 77.40331 78.98373 78.72216 78.28686 79.01578 79.84550 78.97300 77.70144 77.66372 77.07437 77.87492 
# A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
# 76.03819 78.26988 78.38190 78.59127 77.17162 77.88693 
#基因组划分区域进行计算
tiles_500=tileMethylCounts(filtered.myobj,win.size=500,step.size=500,cov.bases =1)
meth_500bp=methylKit::unite(tiles_500, destrand=FALSE,min.per.group=0L,mc.cores = 40)
data_raw_CpG<-as.data.frame(percMethylation(meth_500bp,rowids=T))
meth_500bp<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY",rownames(data_raw_CpG),invert=TRUE),]
head(meth_500bp);tail(meth_500bp);dim(meth_500bp)

colMeans(meth_500bp,na.rm=TRUE)
# A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M      A4M      A5M 
# 74.25854 79.65323 81.03221 78.37740 80.08254 79.71471 79.20225 80.09645 80.94543 79.94314 78.76612 78.73462 78.03263 78.93025 
# A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
# 77.02354 79.25761 79.42741 79.63876 78.09573 78.90264
colMeans(data_raw_CpG,na.rm=TRUE)
# A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M      A4M      A5M 
# 74.25865 79.65136 81.03048 78.37602 80.08918 79.71776 79.20118 80.09542 80.94397 79.94150 78.76480 78.73315 78.03116 78.92938 
# A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
# 77.02218 79.25598 79.42629 79.63728 78.09404 78.90146 

tile_list<-c(list(tiles_100),list(tiles_200),list(tiles_300),list(tiles_500))
meth_list<-c(list(meth_100bp),list(meth_200bp),list(meth_300bp),list(meth_500bp))
saveRDS(tile_list, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_6X_tile_list.rds")
saveRDS(meth_list, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_6X_meth_list.rds")

##final_determination corverage and bins for downstream calculation
#rm(list = ls())
#library(methylKit)
myobj<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_myobj_6X_CpG.rds")
depth<-6
filtered.myobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
#基因组划分区域进行计算
tiles_200=tileMethylCounts(filtered.myobj,win.size=200,step.size=200,cov.bases =3)
meth_200bp=methylKit::unite(tiles_200, destrand=FALSE,min.per.group=0L,mc.cores = 48)
data_raw_CpG<-as.data.frame(percMethylation(meth_200bp,rowids=T))
meth_200bp_2<-data_raw_CpG[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",rownames(data_raw_CpG),invert=TRUE),]
head(meth_200bp_2);tail(meth_200bp_2);dim(meth_200bp_2)#1042601      20
colMeans(meth_200bp_2,na.rm=TRUE)
# A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M 
# 70.23334 75.04619 76.25282 73.88961 74.96154 74.49886 74.50373 74.72829 73.75739 74.25043 73.87372 73.85539 
# A4M      A5M      A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
# 73.75239 74.13798 72.71597 75.68872 75.11009 75.01015 74.94328 74.78978 
colMeans(data_raw_CpG,na.rm=TRUE)
# A2Q      A3Q      A4Q      A5Q      A6Q      Y1Q      Y2Q      Y4Q      Y5Q      Y6Q      A2M      A3M 
# 70.27821 75.10729 76.31767 73.93032 75.02023 74.56319 74.58581 74.81778 73.83107 74.32584 73.95508 73.93657 
# A4M      A5M      A6M      Y1M      Y2M      Y4M      Y5M      Y6M 
# 73.82391 74.20634 72.77231 75.74430 75.19545 75.07674 75.00338 74.84974 
saveRDS(tiles_200, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_6X_tiles_200.rds")
saveRDS(meth_200bp_2, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_6X_meth_200bp_2.rds")
head(as.data.frame(meth_200bp_2))
write.table(as.data.frame(meth_200bp_2), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_meth_level_200bp_CpG_6X_3site_bismark_no_filter.txt",quote=F, row.names=T, col.names=T) 

#plot for view
#tiles_200_2 <- readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/tiles_200.rds")
meth_200bp_2<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_20_samples_6X_meth_200bp_2.rds")
names(meth_200bp_2)
#"A2Q" "A3Q" "A4Q" "A5Q" "A6Q" "Y1Q" "Y2Q" "Y4Q" "Y5Q" "Y6Q" "A2M" "A3M" "A4M" "A5M" "A6M" "Y1M" "Y2M" "Y4M"
#"Y5M" "Y6M"

colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/verification_AMA_analysis_metadata.csv",row.names=1,header =T)
head(colData_whole)

#count 200bp bin in each samples
#filtering and arrange data
evaluation.table<-meth_200bp_2
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
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.number_of_200bp_bins_detected_in_varification_20_samples_by_6X_CpG_3L.pdf",plot_bin200_detected_number1,width = 10, height = 8)


##使用平均值进行绘图
data_region1<-as.data.frame(colMeans(meth_200bp_2,na.rm = T))
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
#  generation .y.       group1 group2       p  p.adj p.format p.signif method
#1 Mom        meth_mean Young  AMA    0.00184 0.0037 0.0018   **       T-test
#2 Kid        meth_mean Young  AMA    0.807   0.81   0.8074   ns       T-test
compare_means(meth_mean~group,group.by = "generation", data=data_region,method = "wilcox.test")
#Sample.types .y.       group1 group2     p p.adj p.format p.signif method
#1Mom        meth_mean Young  AMA    0.00794 0.016 0.0079   **       Wilcoxon
#2 Kid        meth_mean Young  AMA    0.548   0.55  0.5476   ns       Wilcoxon

global_meth_plot1 <- ggboxplot(data_region, x = "group", y = "meth_mean",facet.by = "generation",
                               color = "group", palette = ppCor[2:1],
                               add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(60,90))
#  Add p-value
pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.AMA_group_global_comparison_mean_varification_20_samples_all_site_by_6X_CpG_200bp_3L.pdf")
global_meth_plot1+stat_compare_means(method = "t.test",comparisons= list(c("Young","AMA")))
global_meth_plot1+stat_compare_means(method = "wilcox.test",comparisons= list(c("Young","AMA")))
#p+stat_compare_means(aes(label=..p.signif..), label.x = 1.5, label.y = 1.1)
global_meth_plot1+stat_compare_means(method = "t.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
global_meth_plot1+stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
dev.off()

#grid.newpage(); #清空画板，开始画新图

#plot only for kids
global_meth_plot_k <- ggboxplot(data_region[which(data_region$generation == "Kid"),], x = "group", y = "meth_mean",
                              color = "group", palette = ppCor[2:1],
                              add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(60,90))+
  stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.Kids_AMA_group_global_comparison_mean_varification_20_samples_all_site_by_6X_CpG_200bp_3L.pdf",global_meth_plot_k,width = 6, height = 6)
#plot only for mom
global_meth_plot_m <- ggboxplot(data_region[which(data_region$generation == "Mom"),], x = "group", y = "meth_mean",
                              color = "group", palette = ppCor[2:1],
                              add = "jitter")+xlab("") +ylab("DNA methylation level(CpG 5mC%)")+ ylim(c(60,90))+
  stat_compare_means(method = "wilcox.test",aes(label=paste0(..method..,": "," p = ",..p.format..,"\n",..p.signif..)), label.x = 1.4, label.y = 85)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/0.Mom_AMA_group_global_comparison_mean_varification_20_samples_all_site_by_6X_CpG_200bp_3L.pdf",global_meth_plot_m,width = 6, height = 6)
