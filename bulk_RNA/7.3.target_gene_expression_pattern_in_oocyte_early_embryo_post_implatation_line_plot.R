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

gene_target<-c( "HTRA3", "CD24")
#######for early embryo
#step one:: reading gene expresion matrix and cor genes
#读入表达数据矩阵normal count or log2
data_trans<- as.data.frame(read.table("/mnt/data/chenwei/huahua/4.related_data/NSMB-FPKM-hg19/NSMB_sample_gene_FPKM.xls",header = T,sep = "\t"))
str(data_trans)
head(data_trans)
#colData构建
samplenames<-colnames(data_trans[,2:91])#将要分析的亚集列名字赋予此
Cell_type <- unlist(lapply(strsplit(samplenames,"_"), function(x) x[1]))
Cell_type <- factor(Cell_type,levels = c("Oocyte","Zygote","X2cell","X4cell","X8cell","Morula","Bst") )

colData <- data.frame(subtype = Cell_type,row.names = samplenames)
colData$cell<-rownames(colData)

all(rownames(colData) %in% colnames(data_trans))

#导入基因名
data_expr<-data_trans[which(data_trans$gene_id%in% as.character(gene_target)),]
rownames(data_expr)=data_expr$gene_id
data_expr<-data_expr[,-(1)]
data_expr<-data_expr[which(rowSums(data_expr)>0),]
data_expr<-log2(data_expr+1)

head(data_expr);dim(data_expr)
range(data_expr);length(Cell_type)


#step four 单一基因箱线图绘制
#整合分组信息至表达数据:
data_expr1<- as.data.frame(t(data_expr))
length(colnames(colData))
head(colData)
data_preplot=merge(data_expr1,colData[,c("subtype", "cell")],by=0)
head(data_preplot)
nrow(data_preplot)
#grid.newpage(); #清空画板，开始画新图
data_plot <- melt(data_preplot)#数据整理
head(data_plot)
#R语言学习笔记-reshape2 https://blog.csdn.net/qq_27755195/article/details/45329161

# summarySE 计算标准差和标准误差以及95%的置信区间.
tgc <- summarySE(data_plot, measurevar="value", groupvars= c("subtype","variable"))
head(tgc);tail(tgc)
which(tgc$variable == "cell_order")
genetags<-tgc[which(tgc$subtype=="Oocyte"),]
maintile<-"Early_embryo_expression_pattern_for_correlated_DEGs_Kid_mom"


line_plot <-ggplot(tgc, aes(x=subtype, y=value, colour=variable,group=variable)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),position=position_dodge(0.1), width=0.5) +
  geom_line(position=position_dodge(0.1)) +
  geom_point(size=2,position=position_dodge(0.1))+#控制两线各向左向右移0.2的聚类，线也移，点也移，永远不分开。
  ggtitle(maintile) +xlab("Stage") + ylab("log2(FPKM+1)") +
  scale_colour_hue(name="Gene", l=40) +  scale_color_manual(values=ppCor[3:4])+
  theme_bw() +
  #theme(legend.position = 'top')+   ##同理可以取 bottom、right、left
  theme(plot.title = element_text(hjust=0.5,vjust=0.5,size=10,colour = "black",face = "bold"))+ #face="bold"加粗
  #geom_text(data=genetags,aes(x=0.5,y=value,color=variable,label=variable),size=2.5, alpha=.8,hjust=0,check_overlap = TRUE)+
  theme(plot.title = element_text(hjust=0.5,size=15,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(panel.border = element_rect(colour="black",fill=NA),panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), axis.text.y = element_text(size=15),legend.position = "right")

line_plot2<-line_plot+facet_wrap(~ variable, scales = "free",ncol=length(gene_target)) 
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_early_embryo_expression_pattern_lineplot1.pdf", pointsize=5,width =7,height = 6)
line_plot
dev.off()
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_early_embryo_expression_pattern_lineplot2.pdf", pointsize=5,width =6,height = 4)
line_plot2
dev.off()


##for oocyte

#读入表达数据矩阵normal count or log2
data_trans<- data.frame(read.csv("/mnt/data/chenwei/huahua/4.related_data/Oocyte_development-FPKM-hg19/Molecule_Cell_Folliculogenesis_FPKM_contain_MII_O.csv",header = T,sep = ","))
#colData构建
names(data_trans)[names(data_trans) == 'X'] <- 'gene_id'
samplenames<-colnames(data_trans[,2:152])#将要分析的亚集列名字赋予此
Cell_type <- unlist(lapply(strsplit(samplenames,"_"), function(x) paste(x[1],x[2],sep="_")))
unique(Cell_type)
Cell_type <- factor(Cell_type,levels = c("Primordial_follicle","Primary_follicle","Secondary_follicle","Antral_follicle","Preovulatory_follicle"))
colData <- data.frame(subtype = Cell_type,type = Cell_type,row.names = samplenames)
colData<-colData[order(colData$subtype,decreasing = F),]
all(rownames(colData) %in% colnames(data_trans))
###extraxt expression for target genes
data_expr<-data_trans[which(data_trans$gene_id%in% gene_target),]
rownames(data_expr)=data_expr$gene_id
data_expr<-data_expr[,-(1)]
data_expr<-log2(data_expr+1)
data_expr<-data_expr[,rownames(colData)]
#单一基因箱线图绘制
data_expr1<- as.data.frame(t(data_expr))
length(colnames(colData))
data_preplot=merge(data_expr1,colData,by=0)
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
  geom_point(size=2,position=position_dodge(0.1))+#控制两线各向左向右移0.2的聚类，线也移，点也移，永远不分开。
  ggtitle(maintile) +xlab("Stage") + ylab("log2(FPKM+1)") +
  scale_colour_hue(name="Gene", l=40) +  scale_color_manual(values=ppCor[3:4])+
  theme_bw() +
  #theme(legend.position = 'top')+   ##同理可以取 bottom、right、left
  theme(plot.title = element_text(hjust=0.5,vjust=0.5,size=10,colour = "black",face = "bold"))+ #face="bold"加粗
  #geom_text(data=genetags,aes(x=0.5,y=value,color=variable,label=variable),size=2.5, alpha=.8,hjust=0,check_overlap = TRUE)+
  theme(plot.title = element_text(hjust=0.5,size=15,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(panel.border = element_rect(colour="black",fill=NA),panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), axis.text.y = element_text(size=15),legend.position = "right")

line_plot2<-line_plot+facet_wrap(~ variable, scales = "free",ncol=2) 
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_oocyte_expression_pattern_lineplot1.pdf", pointsize=5,width =7,height = 6)
line_plot
dev.off()
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_oocyte_expression_pattern_lineplot2.pdf", pointsize=5,width =7,height = 4)
line_plot2
dev.off()

#######posted implatation

data_trans<- data.frame(read.table("/mnt/data/chenwei/huahua/4.related_data/Post_implantation_Embryo_TPM/GSE109555_All_Embryo_TPM.txt",header = T,sep = "\t"))
str(data_trans);dim(data_trans)# 22359  5911
head(data_trans)
data_trans[1:6,1:6]
#colData构建
#names(data_trans)[names(data_trans) == 'X'] <- 'gene_id'

library(readxl)
colData <-data.frame(read_excel("/mnt/data/chenwei/huahua/4.related_data/Post_implantation_Embryo_TPM/Sample_Information.xlsx",sheet = 1))
colData<-colData[which(colData$CNV == "Normal" & colData$Lineage != "MIX"),]

colData$stage_lineage<-paste0(colData$Day,"_",colData$Lineage)
colData_used<-colData[,c("Sample","Lineage","stage_lineage","Day","Embryo","IVC","Sex")]
table(colData_used$Lineage)
table(colData_used$Day)

colData_used$Day<-factor(colData_used$Day,levels = c("D6","D8","D10","D12","D14"),ordered = T)
colData_used$Lineage<-factor(colData_used$Lineage,levels = c("TE","EPI","PE"),ordered = T)
colData_used$Lineage<-factor(colData_used$Lineage,levels = c("TE","EPI","PE"),ordered = T)
colData_used$stage_lineage<-factor(colData_used$stage_lineage,levels = c("D6_TE","D8_TE","D10_TE","D12_TE","D14_TE","D6_EPI","D8_EPI","D10_EPI","D12_EPI","D14_EPI","D6_PE","D8_PE","D10_PE","D12_PE","D14_PE"),ordered = T)

rownames(colData_used)<-colData_used$Sample;colData_used<-colData_used[,-1]
head(colData_used)
all(rownames(colData_used) %in% colnames(data_trans))
colData_used<-colData_used[order(colData_used$Lineage,colData_used$stage_lineage,decreasing = F),]
str(colData_used)
table(colData_used$stage_lineage)

data_expr<-data_trans[rownames(data_trans) %in% gene_target,rownames(colData_used)];data_expr[,1:6]
data_expr<-data_expr[which(rowSums(data_expr,na.rm = T)>0),]
range(data_expr)# 0.00 6441.81
data_expr<-log2(data_expr+1)
head(data_expr);dim(data_expr)
range(data_expr)#  0.00000 12.65347
data_expr[1:6,1:6]

#箱线图绘制
#整合分组信息至表达数据:
data_expr1<- as.data.frame(t(data_expr))
length(colnames(colData_used))
data_preplot=merge(data_expr1,colData_used,by=0)
head(data_preplot)
nrow(data_preplot)

#grid.newpage(); #清空画板，开始画新图
data_plot <- melt(data_preplot)#数据整理
head(data_plot)
#R语言学习笔记-reshape2 https://blog.csdn.net/qq_27755195/article/details/45329161

# summarySE 计算标准差和标准误差以及95%的置信区间.
tgc <- summarySE(data_plot, measurevar="value", groupvars= c("Lineage","variable"))
head(tgc);tail(tgc)

genetags<-tgc[which(tgc$Lineage=="PE"),]
maintile<-"post_implantation_expression_pattern_for_correlated_DEGs_Kid_mom"

line_plot <-ggplot(tgc, aes(x=Lineage, y=value, colour=variable,group=variable)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),position=position_dodge(0.1), width=0.5) +
  geom_line(position=position_dodge(0.1)) +
  geom_point(size=2,position=position_dodge(0.1))+#控制两线各向左向右移0.2的聚类，线也移，点也移，永远不分开。
  ggtitle(maintile) +xlab("Stage") + ylab("log2(TPM+1)") +
  scale_colour_hue(name="Gene", l=40) +  scale_color_manual(values=ppCor[3:4])+
  theme_bw() +
  #theme(legend.position = 'top')+   ##同理可以取 bottom、right、left
  theme(plot.title = element_text(hjust=0.5,vjust=0.5,size=10,colour = "black",face = "bold"))+ #face="bold"加粗
  #geom_text(data=genetags,aes(x=0.5,y=value,color=variable,label=variable),size=2.5, alpha=.8,hjust=0,check_overlap = TRUE)+
  theme(plot.title = element_text(hjust=0.5,size=15,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(panel.border = element_rect(colour="black",fill=NA),panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), axis.text.y = element_text(size=15),legend.position = "right")

line_plot2<-line_plot+facet_wrap(~ variable, scales = "free",ncol=length(gene_target)) 
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_post_implatation_Lineage_expression_pattern_lineplot1.pdf", pointsize=5,width =7,height = 6)
line_plot
dev.off()
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_post_implatation_Lineage_expression_pattern_lineplot2.pdf", pointsize=5,width =7,height = 4)
line_plot2
dev.off()


#for stage 
head(data_plot)
# summarySE 计算标准差和标准误差以及95%的置信区间.
tgc <- summarySE(data_plot, measurevar="value", groupvars= c("Day","variable"))
head(tgc);tail(tgc)

genetags<-tgc[which(tgc$Day =="D6"),]
maintile<-"Day_post_implantation_expression_pattern_for_correlated_DEGs_Kid_mom"
line_plot <-ggplot(tgc, aes(x=Day, y=value, colour=variable,group=variable))+ 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),position=position_dodge(0.1), width=0.5) +
  geom_line(position=position_dodge(0.1)) +
  geom_point(size=2,position=position_dodge(0.1))+#控制两线各向左向右移0.2的聚类，线也移，点也移，永远不分开。
  ggtitle(maintile) +xlab("Stage") + ylab("log2(TPM+1)") +
  scale_colour_hue(name="Gene", l=40) +  scale_color_manual(values=ppCor[3:4])+
  theme_bw() +
  #theme(legend.position = 'top')+   ##同理可以取 bottom、right、left
  theme(plot.title = element_text(hjust=0.5,vjust=0.5,size=10,colour = "black",face = "bold"))+ #face="bold"加粗
  #geom_text(data=genetags,aes(x=0.5,y=value,color=variable,label=variable),size=2.5, alpha=.8,hjust=0,check_overlap = TRUE)+
  theme(plot.title = element_text(hjust=0.5,size=15,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(panel.border = element_rect(colour="black",fill=NA),panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), axis.text.y = element_text(size=15),legend.position = "right")

line_plot2<-line_plot+facet_wrap(~ variable, scales = "free",ncol=length(gene_target)) 
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_post_implatation_Day_expression_pattern_lineplot1.pdf", pointsize=5,width =7,height = 6)
line_plot
dev.off()
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_post_implatation_Day_expression_pattern_lineplot2.pdf", pointsize=5,width =7,height = 4)
line_plot2
dev.off()

#for stage 
head(data_plot)
table(data_plot$stage_lineage)
# summarySE 计算标准差和标准误差以及95%的置信区间.
tgc <- summarySE(data_plot, measurevar="value", groupvars= c("stage_lineage","variable"))
head(tgc);tail(tgc)
tgc<-tgc[which(tgc$N>=3),]
tgc$stage_lineage<-factor(tgc$stage_lineage,levels = c("D6_EPI","D8_EPI","D10_EPI","D12_EPI",    #"D14_EPI",
                                                       "D6_PE","D8_PE","D10_PE","D12_PE",       #"D14_PE",
                                                       "D6_TE","D8_TE","D10_TE","D12_TE","D14_TE"),ordered = T)

genetags<-tgc[which(tgc$stage_lineage =="D6_EPI"),]
maintile<-"stage_lineage_post_implantation_expression_pattern_for_correlated_DEGs_Kid_mom"

line_plot <-ggplot(tgc, aes(x=stage_lineage, y=value, colour=variable,group=variable)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se),position=position_dodge(0.1), width=0.5) +
  geom_line(position=position_dodge(0.1)) +
  geom_point(size=2,position=position_dodge(0.1))+#控制两线各向左向右移0.2的聚类，线也移，点也移，永远不分开。
  ggtitle(maintile) +xlab("Stage") + ylab("log2(TPM+1)") +
  scale_colour_hue(name="Gene", l=40) +  scale_color_manual(values=ppCor[3:4])+
  theme_bw() +
  #theme(legend.position = 'top')+   ##同理可以取 bottom、right、left
  theme(plot.title = element_text(hjust=0.5,vjust=0.5,size=10,colour = "black",face = "bold"))+ #face="bold"加粗
  #geom_text(data=genetags,aes(x=0.5,y=value,color=variable,label=variable),size=2.5, alpha=.8,hjust=0,check_overlap = TRUE)+
  theme(plot.title = element_text(hjust=0.5,size=15,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(panel.border = element_rect(colour="black",fill=NA),panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), axis.text.y = element_text(size=15),legend.position = "right")

line_plot2<-line_plot+facet_wrap(~ variable, scales = "free",ncol=length(gene_target)) 
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_post_implatation_stage_lineage_expression_pattern_lineplot1.pdf", pointsize=5,width =7,height = 6)
line_plot
dev.off()
pdf("/mnt/data/chenwei/huahua/manuscript/figure/figure4/target_gene_post_implatation_stage_lineage_expression_pattern_lineplot2.pdf", pointsize=5,width =7,height = 4)
line_plot2
dev.off()
