rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(pheatmap) 
library(dplyr)
#library(export)
require(gridExtra)
library(reshape2)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(scales)
library(stringr)
library(RColorBrewer)
library(dendextend)
library(ComplexHeatmap)
library(circlize)
##提供自定义函数
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # 计算长度
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # 以 groupvars 为组,计算每组的长度,均值,以及标准差
  # ddply 就是 dplyr 中的 group_by + summarise
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # 重命名  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  # 计算标准偏差
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  # 计算置信区间
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#读入表达数据矩阵normal count or log2
data_trans<- data.frame(read.csv("/mnt/data/chenwei/huahua/4.related_data/Oocyte_development-FPKM-hg19/Molecule_Cell_Folliculogenesis_FPKM_contain_MII_O.csv",header = T,sep = ","))
str(data_trans)
head(data_trans)
dim(data_trans)
#colData构建
names(data_trans)[names(data_trans) == 'X'] <- 'gene_id'

samplenames<-colnames(data_trans[,2:152])#将要分析的亚集列名字赋予此
Cell_type <- unlist(lapply(strsplit(samplenames,"_"), function(x) paste(x[1],x[2],sep="_")))
unique(Cell_type)
Cell_type <- factor(Cell_type,levels = c("Primordial_follicle","Primary_follicle","Secondary_follicle","Antral_follicle","Preovulatory_follicle"))

colData <- data.frame(subtype = Cell_type,type = Cell_type,row.names = samplenames)
colData<-colData[order(colData$subtype,decreasing = F),]

all(rownames(colData) %in% colnames(data_trans))

#导入基因名
cor_DEG_mom_kids_data<-read.table(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/correlationship_trans_level_of_AMA_common_DEG_between_kids_and_mother.txt",header=T,sep="\t")
head(cor_DEG_mom_kids_data)
gene_target<-unique(as.character(cor_DEG_mom_kids_data[which(cor_DEG_mom_kids_data$group == "candidate_region"),]$Gene))

#gene_target<-c("SMPD1","SMPDL3A","SMPDL3B")
#gene_target<-c("SMPDL3A","SCD")
#gene_target<-c("PSENEN", "SCD", "CD24")
data_expr<-data_trans[which(data_trans$gene_id%in% gene_target),]
rownames(data_expr)=data_expr$gene_id
data_expr<-data_expr[,-(1)]
data_expr<-log2(data_expr+1)
data_expr<-data_expr[,rownames(colData)]

head(data_expr);dim(data_expr)
range(data_expr)

length(Cell_type)
pheatmap(data_expr, cluster_rows=F, cluster_cols =F,scale="row", annotation_col=colData)

#手动scale
DEG_expr_scal<-scale(t(data_expr), center=T,scale=T)
range(DEG_expr_scal)#-2.934484 11.505842
pheatmap(t(DEG_expr_scal), cluster_rows=T, cluster_cols =F,scale="none", annotation_col=colData)
pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row",annotation_col=colData,display_numbers = TRUE,color = colorRampPalette(colors = c("blue","white","red"))(100))

heat_plot<-pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row",annotation_col=colData,color = colorRampPalette(colors = c("blue","white","red"))(100))
#heat_plot<-pheatmap(data_expr, cluster_rows=T, cluster_cols =T,scale="row",annotation_col=colData,color = colorRampPalette(colors = c("blue","white","red"))(100))
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/oocyte_expression_pattern_for_correlated_DEGs_Kid_mom_heatmaps.pdf",width = 25,height =10)
print(heat_plot)
dev.off()
#breaks
bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))
# 做热图：
heat_plot<-pheatmap(data_expr, cluster_rows=T, cluster_cols =F,scale="row",annotation_col=colData,
                    color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                    legend_breaks=seq(-4,4,2),show_colnames = FALSE,
                    breaks=bk)
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/oocyte_expression_pattern_for_correlated_DEGs_Kid_mom_heatmaps2.pdf",width = 20,height =10)
print(heat_plot)
dev.off()


#单一基因箱线图绘制
#整合分组信息至表达数据:
data_expr1<- as.data.frame(t(data_expr))
length(colnames(colData))
data_preplot=merge(data_expr1,colData,by=0)
head(data_preplot)
nrow(data_preplot)


#grid.newpage(); #清空画板，开始画新图
data_plot <- melt(data_preplot)#数据整理
head(data_plot)
#R语言学习笔记-reshape2 https://blog.csdn.net/qq_27755195/article/details/45329161

# summarySE 计算标准差和标准误差以及95%的置信区间.
tgc <- summarySE(data_plot, measurevar="value", groupvars= c("subtype","variable"))
head(tgc);tail(tgc)

genetags<-tgc[which(tgc$subtype=="Primordial_follicle"),]
maintile<-"oocyte_expression_pattern_for_correlated_DEGs_Kid_mom"

line_plot <-ggplot(tgc, aes(x=subtype, y=value, colour=variable,group=variable)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),position=position_dodge(0.1), width=0.5) +
  geom_line(position=position_dodge(0.1)) +
  #geom_point(size=2,position=position_dodge(0.1))+#控制两线各向左向右移0.2的聚类，线也移，点也移，永远不分开。
  ggtitle(maintile) +xlab("Stage") + ylab("log2(FPKM+1)") +
  #expand_limits(x=-.2) +                        # Expand y range
  #  scale_y_continuous(breaks=0:30*1) +# Set tick every 4
  #scale_y_continuous(breaks=0:40*1) +# Set tick every 4
  scale_colour_hue(name="Gene", l=40) +  #  scale_color_manual(values=ppCor[3:5])+
  theme_bw() +
  #  theme(legend.position = 'top')+   ##同理可以取 bottom、right、left
  theme(plot.title = element_text(hjust=0.5,vjust=0.5,size=10,colour = "black",face = "bold"))+ #face="bold"加粗
  #geom_text(data=genetags,aes(x=0.5,y=value,color=variable,label=variable),size=2.5, alpha=.8,hjust=0,check_overlap = TRUE)+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(panel.border = element_rect(colour="black",fill=NA),panel.grid=element_blank(),
     axis.title = element_text(size=15),axis.text.x = element_text(size=10), axis.text.y = element_text(size=15),legend.position = "right")

  line_plot2<-line_plot+facet_wrap(~ variable, scales = "free",ncol=6) 
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/oocyte_expression_pattern_for_correlated_DEGs_Kid_mom_lineplot1.pdf", pointsize=10,width =10,height = 10)
line_plot
dev.off()
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/oocyte_expression_pattern_for_correlated_DEGs_Kid_mom_lineplot2.pdf", pointsize=5,width =20,height = 30)
line_plot2
dev.off()


#趋势分析(STC, Series Test of Cluster)#ref:https://www.jianshu.com/p/ba78fbd05f3f
#https://www.jianshu.com/p/2067864259c5
library(TCseq)#1）调整聚类方法 2）调整类别数目
data_expr<-data_trans[which(data_trans$gene_id%in% gene_target),]
rownames(data_expr)=data_expr$gene_id
data_expr<-data_expr[,-(1)]
data_expr2<-t(data_expr[which(rowSums(data_expr)>0),])
head(colData)
data_expr2<-merge(data_expr2,colData,by=0)
head(data_expr2)

#分组求均值
data_expr3<-aggregate(data_expr2[,2:(ncol(data_expr2)-2)], by=list(data_expr2$subtype), FUN=mean,na.rm= TRUE)
#把第一列设为行名
row.names(data_expr3)<-data_expr3[,1]
data_expr3<-data.frame(t(data_expr3[,-1]))
data_expr3<-log2(data_expr3+1)

head(data_expr3);dim(data_expr3)
plot_matrix <-data.frame(t(scale(t(data_expr3), center = TRUE, scale = TRUE)))
plot_matrix
write.table(plot_matrix, '/mnt/data/chenwei/huahua/manuscript/figure/figure4/table/common_cor_DEGs_zscore_data_oocyte_expression_time_cluster_pattern_for_correlated_DEGs_Kid_mom.txt', sep = '\t', col.names = NA, quote = FALSE)

##绘图
data_expr<-read.table('/mnt/data/chenwei/huahua/manuscript/figure/figure4//table/common_cor_DEGs_zscore_data_oocyte_expression_time_cluster_pattern_for_correlated_DEGs_Kid_mom.txt', header = T,row.names = 1)
head(data_expr)
#plot_matrix<-data_expr[,-ncol(data_expr)]
plot_matrix[1:4,]

ann_colors = list(cell_stage=c(Primordial_follicle=ppCor[1],Primary_follicle=ppCor[2],Secondary_follicle=ppCor[3],Antral_follicle=ppCor[4],Preovulatory_follicle=ppCor[5]))
anno_cell<-c("Primordial_follicle","Primary_follicle","Secondary_follicle", "Antral_follicle", "Preovulatory_follicle")
column_ha = HeatmapAnnotation(cell_stage = anno_cell,col = ann_colors)

htkm0 <- Heatmap(plot_matrix,name= "z-score", border = TRUE,
                 top_annotation = column_ha,
                 cluster_column_slices = FALSE,column_title = "Gene expression in Human oocyte stages for correlated_DEGs_between_mom_kids", 
                 col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                 show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                 row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,
                 row_dend_width = unit(4, "cm")
)
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/mean_oocyte_expression_pattern_for_correlated_DEGs_Kid_mom_heatmaps.pdf",width = 6,height =6)
print(htkm0)
dev.off()


dend = as.dendrogram(hclust(dist(plot_matrix),method="ward.D2"))
d_num<-10;dend = color_branches(dend, k =d_num)
htkm <- Heatmap(plot_matrix,name= "z-score", border = TRUE,
                top_annotation = column_ha,
                cluster_column_slices = FALSE,column_title = "Gene expression in Human oocyte stages for correlated_DEGs_between_mom_kids", 
                col= colorRamp2(seq(from=-5,to=5,length=100),colorRampPalette(colors = c("blue","white","red"))(100)),
                show_row_names= TRUE,show_column_names= FALSE,row_names_gp= gpar(fontsize = 6),
                row_title_rot= 0,cluster_row_slices= FALSE,cluster_columns= FALSE,#cluster_rows= FALSE,
                cluster_rows = dend, 
                row_split = d_num,row_gap = unit(rep(2, d_num-1), "mm"),
                row_dend_width = unit(4, "cm")
)
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/split_mean_oocyte_expression_pattern_for_correlated_DEGs_Kid_mom_heatmaps.pdf",width = 6,height =6)
print(htkm)
dev.off()

###step three :: plot trend line for gene expression of each cluster
#exact genes names in each cluster
clusterlist = row_order(htkm)
names(clusterlist)<-1:d_num
htkm_module <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(plot_matrix[clusterlist[[i]],]),Cluster = paste0("cluster", i),stringsAsFactors = FALSE)
  return(out)}) %>%   do.call(rbind, .)

head(htkm_module);colnames(htkm_module)<-c("row_genes","Cluster")
table(htkm_module$Cluster)
#cluster1 cluster10  cluster2  cluster3  cluster4  cluster5  cluster6  cluster7  cluster8  cluster9 
#   3        10         4         4         7         6         3         4         4         3 
head(plot_matrix)
data_plot<-data.frame(plot_matrix[htkm_module$row_genes,])
data_plot$type<-as.character(htkm_module$Cluster)
cell_names<-colnames(data_plot)
data_plot[1:4,]
data_plot$gene<-rownames(data_plot)
data_preplot <- melt(data_plot,variable.name="cell",value.name = "expression",id.vars = c("type","gene"))
head(data_preplot)

data_preplot3<-data_preplot
range(data_preplot3$expression)#-1.738838  1.787046
data_preplot3<-data_preplot3[order(data_preplot3$cell),]
head(data_preplot3)
data_plot_A<-data_preplot3
#data_plot_B<-aggregate(expression ~ type+cell, data = data_preplot, mean)
data_plot_B<-aggregate(expression ~ type+cell, data = data_preplot, median)

data_plot_B<-data_plot_B[order(data_plot_B$cell),]
head(data_plot_A);head(data_plot_B)
#colnames(data_plot_B)<-c("group","cell","expression")
trend_plot<-ggplot(data_plot_A,aes(x= cell,y=expression))+ 
  geom_point(aes(colour = cell),size=1,shape=19,alpha = 0.5)+
  geom_line(aes(group=gene), color="gray" ,position = position_dodge(0.02)) +
  geom_line(data=data_plot_B,aes(group=type),alpha=0.7,colour="steelblue",size=1.2,linetype=1)+
  scale_colour_manual(values=ppCor)+#ylim(-2,7)+
  theme_bw()+ theme(plot.title = element_text(hjust = 0.5),axis.text.y=element_blank(),
                    axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),
                    panel.grid=element_blank(),axis.ticks = element_blank())+
  xlab("pseudotime") + ylab("gene expression") + facet_wrap(~type,scales="free_y",ncol = 2)
trend_plot2<-trend_plot + stat_summary(fun = "median", color = "red", size = 1, geom = "point",alpha=0.7)
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure4/mean_expression_trends_for_different_cluster_Human_oocyte_stages_for_correlated_DEGs_between_mom_kids.pdf",trend_plot2,width =6, height =8,limitsize = FALSE)

#构建对象
data_expr3a<-as.matrix(data_expr3)
set.seed(19921010)
cluster_num <- 9
tcseq_cluster_cm <- timeclust(data_expr3a, algo = 'cm', k = cluster_num, standardize = TRUE)
#timeclust()是一个整合函数，可执行数据标准化、聚类等多步操作，将上述输入数据中具有相似的时间表达特征的蛋白聚在一类。#standardize 用于 z-score 标准化变量
#基于模糊c均值聚类（timeclust()参数algo='cm'）的原理对蛋白质表达值的时间序列进行了聚类。timeclust()还提供了其它的聚类算法，如层次聚类（参数algo='hc'）、k均值划分（参数algo='km'）、围绕中心点划分（参数algo='pam'）等。
#颜色、线宽、坐标轴、字体等细节可以在函数中调整，具体参数详见函数帮助
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/timecplot_cm_oocyte_develepment_expression_pattern_for_correlated_DEGs_Kid_mom.pdf", width =20,height = 20)
timecplot_cm <- timeclustplot(tcseq_cluster_cm, value = 'z-score', cols = 3, 
                              axis.line.size = 0.6, axis.title.size = 8, axis.text.size =5, 
                              title.size = 8, legend.title.size = 8, legend.text.size = 8)
dev.off()

tcseq_cluster_hc <- timeclust(data_expr3a, algo = 'hc', k = cluster_num, standardize = TRUE)
#timeclust()是一个整合函数，可执行数据标准化、聚类等多步操作，将上述输入数据中具有相似的时间表达特征的蛋白聚在一类。#standardize 用于 z-score 标准化变量
#基于模糊c均值聚类（timeclust()参数algo='cm'）的原理对蛋白质表达值的时间序列进行了聚类。timeclust()还提供了其它的聚类算法，如层次聚类（参数algo='hc'）、k均值划分（参数algo='km'）、围绕中心点划分（参数algo='pam'）等。
#颜色、线宽、坐标轴、字体等细节可以在函数中调整，具体参数详见函数帮助
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/timecplot_hc_oocyte_develepment_expression_pattern_for_correlated_DEGs_Kid_mom.pdf", width =20,height = 20)
timecplot_hc <- timeclustplot(tcseq_cluster_hc, value = 'z-score', cols = 3,
                              axis.line.size = 0.6, axis.title.size = 8, axis.text.size = 8, 
                              title.size = 8, legend.title.size = 8, legend.text.size = 8)
dev.off()
#上述获得了 10 组聚类群 #如果绘制单个的聚类群
timecplot_cm[2]

#tcseq_cluster<-tcseq_cluster_cm
tcseq_cluster<-tcseq_cluster_hc

#查看每个蛋白所属的聚类群，展示前几个为例
head(tcseq_cluster@cluster)
#统计每个聚类群中各自包含的蛋白数量
table(tcseq_cluster@cluster)

#聚类过程通过计算 membership 值判断蛋白质所属的聚类群，以最大的 membership 值为准
head(tcseq_cluster@membership)

#标准化后的表达值(绘制曲线图值）
head(tcseq_cluster@data)

#最后，提取所有蛋白所属的聚类群，并和它们的原始表达值整合在一起
express_cluster <- tcseq_cluster@cluster
data_expr3b <- cbind(data_expr3a[names(express_cluster), ], express_cluster)
head(data_expr3b)
write.table(data_expr3b, '/mnt/data/chenwei/huahua/manuscript/figure/figure4//table/hc_km_raw_data_oocyte_expression_time_cluster_pattern_for_correlated_DEGs_Kid_mom.txt', sep = '\t', col.names = NA, quote = FALSE)
zscore_data<-data.frame(tcseq_cluster@data)
data_expr3a <- cbind(zscore_data[names(express_cluster), ], express_cluster)
write.table(data_expr3a, '/mnt/data/chenwei/huahua/manuscript/figure/figure4//table/hc_km_zscore_data_oocyte_expression_time_cluster_pattern_for_correlated_DEGs_Kid_mom.txt', sep = '\t', col.names = NA, quote = FALSE)
