rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(VennDiagram)
library(scales)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(UpSetR)

##提供自定义函数
#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor_all)
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


#FOR kids
#read expression data 
tag_name<-"Kid"; sampleA <-"AMA";sampleB <-"Young" 
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file9_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
kids_merge_data <- read.table(file=file,header=T)
head(kids_merge_data)
rownames(kids_merge_data)<-kids_merge_data$Row.names
kids_merge_data<-kids_merge_data[,-1]
kids_merge_data[1:6,1:6]
kids_draw_count<-kids_merge_data[,1:10]
#read DEGs
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file13_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
kids_AMA_DEGs<-read.table(file)
head(kids_AMA_DEGs)
kids_AMA_DEGs_all<-kids_AMA_DEGs$ID
kids_AMA_DEGs_up_1.5<-kids_AMA_DEGs[which(kids_AMA_DEGs$log2FoldChange>0),]$ID
kids_AMA_DEGs_down_1.5<-kids_AMA_DEGs[which(kids_AMA_DEGs$log2FoldChange<0),]$ID

#DNA_Meth_genes<-c("DNMT3A","DNMT3B","DNMT3C","DNMT3L","DNMT1","TET1","TET2","TET3","CBP","GLP","GNATS","MOZ","YBF2","SAS3","SAS2","TIP60","UTX1")
#DNA_Meth_genes<-c("DNMT3A","DNMT3B","DNMT3C","DNMT3L","DNMT1","TET1","TET2","TET3","UHRF1","TDG")
DNA_Meth_genes<-c("DNMT3A","DNMT3B","DNMT3C","DNMT3L","DNMT1","TET1","TET2","TET3","UHRF1","TDG","DPPA3",
                  "GADD45A","CXXC1","RNF4","GADD45B","KMT2A","CREBBP","SETDB1","PARP1","H2AFX","UHRF1","XRCC1",
                  "MBD4","SMUG1","AICDA","APOBEC2")

Reduce(intersect,list(kids_AMA_DEGs_all,DNA_Meth_genes))
# "DNMT1"  "CREBBP" "UHRF1"  "CXXC1"  "H2AFX" 

kids_AMA_DEGs_all[order(kids_AMA_DEGs_all,decreasing = T)]
kids_AMA_DEGs_all[grepl("HDAC",kids_AMA_DEGs_all)]
#HDAC7  HDAC4  HDAC5  HDAC10

kids_AMA_DEGs_all[grepl("HAT",kids_AMA_DEGs_all)]#HAT1
kids_AMA_DEGs_all[grepl("JARID",kids_AMA_DEGs_all)]#JARID2
kids_AMA_DEGs_all[grepl("KDM",kids_AMA_DEGs_all)]#"KDM2A KDM6B KDM4B KDM3B KDM8 
kids_AMA_DEGs_all[grepl("JMJD",kids_AMA_DEGs_all)]
kids_AMA_DEGs_all[grepl("LSD",kids_AMA_DEGs_all)]
kids_AMA_DEGs_all[grepl("MYST",kids_AMA_DEGs_all)]

#read imprinting genes
imprinting_gene0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/human_imprint_gene_single.csv",header=T)
head(imprinting_gene0)
imprinting_genes<-unique(as.character(imprinting_gene0$Gene))
length(imprinting_genes)#228
Reduce(intersect,list(kids_AMA_DEGs_all,imprinting_genes))
# [1] "USMG5"   "FRG1"    "LPAR6"   "ALKBH3"  "LRRK1"   "SH3BP2"  "DENND3"  "DNMT1"   "TRAPPC9" "PTK2B"  
# [11] "ACCS"    "CUL7"    "AGO1"    "CACNA1E" "ANO8"    "ZFAT"    "ACOT11"  "GNG7"    "PDPR" 
venn_Age_imprinting_1<-venn.diagram(list(imprinted_genes=imprinting_genes,kids_AMA_DEGs_all=kids_AMA_DEGs_all),
                                  alpha=c(0.9,0.9),lwd=1,lty=1,col="black" , 
                                  fill=ppCor[3:4],cex = 1.5,cat.col=ppCor[3:4],cat.fontface=4, cat.cex = 1.5,    
                                  main = "imprinted_genes and Kid AMA DEGs",
                                  main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                  filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_kids_RNA_DEGs_imprinting_gene.pdf",width=6,height = 5)
grid.newpage()
grid.draw(venn_Age_imprinting_1)
dev.off()
#load GenAge database genes
longevity_gene0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/longevity.csv")
longevity_gene<-unique(as.character(unlist(strsplit(as.character(longevity_gene0$Gene.s.),","))))
length(longevity_gene)
Age_gene_lgh0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/human_Age_gene_new_500.csv")
Age_gene_lgh<-unique(as.character(Age_gene_lgh0$Symbol))

genAge_database_extand<-unique(c(Age_gene_lgh,longevity_gene));length(genAge_database_extand) #1223

length(Reduce(intersect,list(kids_AMA_DEGs_up_1.5,Age_gene_lgh)))#21
length(Reduce(intersect,list(kids_AMA_DEGs_down_1.5,Age_gene_lgh)))#81
venn_Age_longevity <-venn.diagram(list(Age_gene=Age_gene_lgh,longevity_gene=longevity_gene,kids_AMA_DEGs_all=kids_AMA_DEGs_all),
                                  alpha=c(0.9,0.9,0.9),lwd=1,lty=1,col="black" , 
                                  fill=ppCor[c(10,11,9)],cex = 1.5,cat.col=ppCor[c(10,11,9)],cat.fontface=4, cat.cex = 1.5,    
                                  main = "Database:: longevity gene and two list of Age gene",
                                  main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                  filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_kids_RNA_DEGs_longevity_Ageing_gene.pdf",width=6,height = 6)
grid.newpage();
grid.draw(venn_Age_longevity)
dev.off()

venn_Age_DEGs_genAge1<-venn.diagram(list(genAge_database=genAge_database_extand,kids_DEGs_up=kids_AMA_DEGs_up_1.5,kids_DEGs_down=kids_AMA_DEGs_down_1.5),
                                   alpha=c(1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(5,1:2)], cex = 1.5, 
                                   cat.col=ppCor[c(5,1:2)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_kids_RNA_DEGs_longevity_Ageing_gene2.pdf",width=6,height = 6)
grid.newpage()
grid.draw(venn_Age_DEGs_genAge1)
dev.off()


## SASP marker genes in DEGs
SASP_gene0 <- read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/SASP_human.csv")
SASP_gene<-unique(as.character(SASP_gene0$V1))
Reduce(intersect,list(kids_AMA_DEGs_all,SASP_gene))
#  [1] "GMFG"     "IGF2R"    "CSF1"     "TGFB1"    "CSF2RB"   "VEGFA"    "TIMP2"    "TNFRSF1B" "PTGES"   
# [10] "MMP9"     "MMP14"    "CXCL2" 
Reduce(intersect,list(kids_AMA_DEGs_up_1.5,SASP_gene))
#[1] "GMFG"
Reduce(intersect,list(kids_AMA_DEGs_down_1.5,SASP_gene))
#[1] "IGF2R"    "CSF1"     "TGFB1"    "CSF2RB"   "VEGFA"    "TIMP2"    "TNFRSF1B" "PTGES"    "MMP9"    
#[10] "MMP14"    "CXCL2" 
venn_Age_DEGs_SASP_gene1<-venn.diagram(list(SASP_gene=SASP_gene,kids_DEGs_all=kids_AMA_DEGs_all),
                                   alpha=c(1,1),lwd=1,lty=1,  col="black" , fill=ppCor[c(10,6)], 
                                   cex = 1.5, cat.col=ppCor[c(10,6)],cat.fontface=4,cat.cex = 1.5,
                                   main.cex = 2, main.fontface = 2, main.fontfamily = 3,  filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_kids_RNA_DEGs_SASP_gene1.pdf",width=6,height = 5)
grid.newpage()
grid.draw(venn_Age_DEGs_SASP_gene1)
dev.off()

##human TFs and TFs_cofactor 
Homo_sapiens_TF <- read.table("/mnt/data/chenwei/qinmen_BR/00.ref_data/Homo_sapiens_TF.txt", sep = "\t", header = T)
Homo_sapiens_coTF <- read.table("/mnt/data/chenwei/qinmen_BR/00.ref_data/Homo_sapiens_TF_cofactors.txt", sep = "\t", header = T)
Hm_TFs<-unique(as.character(Homo_sapiens_TF$Symbol))
Hm_coTFs<-unique(as.character(Homo_sapiens_coTF$Symbol))

Reduce(intersect,list(kids_AMA_DEGs_all,Hm_TFs))
Reduce(intersect,list(kids_AMA_DEGs_all,Hm_coTFs))
Reduce(intersect,list(kids_AMA_DEGs_up_1.5,Hm_coTFs))
Reduce(intersect,list(kids_AMA_DEGs_down_1.5,Hm_coTFs))

grid.newpage(); #清空画板，开始画新图
venn_Age_DEGs_TFs<-venn.diagram(list(kids_DEGs_up=kids_AMA_DEGs_up_1.5,kids_DEGs_down=kids_AMA_DEGs_down_1.5,Homo_sapiens_TF=Hm_TFs,Homo_sapiens_coTF=Hm_coTFs),
                                   alpha=c(1,1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1:3,5)],cex = 1.5, cat.col=ppCor[c(1:3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_kids_RNA_DEGs_TFs_gene1.pdf",width=8,height =6)
grid.newpage()
grid.draw(venn_Age_DEGs_TFs)
dev.off()

#for RNA expression
head(kids_draw_count)
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 5),rep("Young", 5))))
rownames(annotation_col) = colnames(kids_draw_count)
## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))
##展示行或者列的label

Gene_expression<-kids_draw_count[which(rownames(kids_draw_count) %in% DNA_Meth_genes),]
heat_plot1<-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                     #  clustering_distance_rows ="euclidean",#correlation
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors, 
                     #gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                     treeheight_col = 20, #treeheight_row = 30, 
                     # labels_row = labels_row,
                     #border_color ="red", 
                     border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     #   filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps1.pdf"),
                     main ="kids_Gene expression of selected meth_trans_inter Genes_merge:original value",angle_col ="45")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors, 
                      # gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                      treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      #  filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps2.pdf"),
                      main ="kids_Gene expression level of meth_trans_inter Genes: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/5.valificat_result/kids_RNA_gene_epigenetic_heatmap1.pdf",width=6,height = 6)
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/5.valificat_result/kids_RNA_gene_epigenetic_heatmap2.pdf",width=6,height = 6)
print(heat_plot2)
dev.off()

#imprinting_genes
Gene_expression<-kids_draw_count[which(rownames(kids_draw_count) %in% imprinting_genes),]
heat_plot1<-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                     #  clustering_distance_rows ="euclidean",#correlation
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors, 
                     #gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                     treeheight_col = 20, #treeheight_row = 30, 
                     # labels_row = labels_row,
                     #border_color ="red", 
                     border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     #   filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps1.pdf"),
                     main ="kids_Gene expression of selected imprinting_genes:original value",angle_col ="45")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =T,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors, 
                      # gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                      treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      #  filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps2.pdf"),
                      main ="kids_Gene expression level of imprinting_genes: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/5.valificat_result/kids_imprinting_genes_heatmap1.pdf",width=12,height =15)
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/5.valificat_result/kids_imprinting_genes_heatmap2.pdf",width=12,height =15)
print(heat_plot2)
dev.off()

#SASP_gene
Gene_expression<-kids_draw_count[which(rownames(kids_draw_count) %in% SASP_gene),]
heat_plot1<-pheatmap(Gene_expression, cluster_row =T,cluster_col =FALSE,na_col = "black",
                     #  clustering_distance_rows ="euclidean",#correlation
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors, 
                     #gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                     treeheight_col = 20, #treeheight_row = 30, 
                     # labels_row = labels_row,
                     #border_color ="red", 
                     border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     #   filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps1.pdf"),
                     main ="kids_Gene expression of  SASP_gene:original value",angle_col ="45")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =T,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors, 
                      # gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                      treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      #  filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps2.pdf"),
                      main ="kids_Gene expression level of SASP_gene : scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/5.valificat_result/kids_SASP_gene_heatmap1.pdf",width=10,height =12)
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/5.valificat_result/kids_SASP_gene_heatmap2.pdf",width=10,height =12)
print(heat_plot2)
dev.off()

#FOR Mom
#read expression data 
tag_name<-"Mom"; sampleA <-"AMA";sampleB <-"Young" 
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file9_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
Mom_merge_data <- read.table(file=file,header=T)
head(Mom_merge_data)
rownames(Mom_merge_data)<-Mom_merge_data$Row.names
Mom_merge_data<-Mom_merge_data[,-1]
Mom_merge_data[1:6,1:6]
Mom_draw_count<-Mom_merge_data[,1:10]

#for methylation related genes 
#read DEGs  
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file13_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
Mom_AMA_DEGs<-read.table(file)
head(Mom_AMA_DEGs)
Mom_AMA_DEGs_all<-Mom_AMA_DEGs$ID
Mom_AMA_DEGs_up_1.5<-Mom_AMA_DEGs[which(Mom_AMA_DEGs$log2FoldChange>0),]$ID
Mom_AMA_DEGs_down_1.5<-Mom_AMA_DEGs[which(Mom_AMA_DEGs$log2FoldChange<0),]$ID


#DNA_Meth_genes<-c("DNMT3A","DNMT3B","DNMT3C","DNMT3L","DNMT1","TET1","TET2","TET3","CBP","GLP","GNATS","MOZ","YBF2","SAS3","SAS2","TIP60","UTX1")
#DNA_Meth_genes<-c("DNMT3A","DNMT3B","DNMT3C","DNMT3L","DNMT1","TET1","TET2","TET3","UHRF1","TDG")
DNA_Meth_genes<-c("DNMT3A","DNMT3B","DNMT3C","DNMT3L","DNMT1","TET1","TET2","TET3","UHRF1","TDG","DPPA3",
                  "GADD45A","CXXC1","RNF4","GADD45B","KMT2A","CREBBP","SETDB1","PARP1","H2AFX","UHRF1","XRCC1",
                  "MBD4","SMUG1","AICDA","APOBEC2")


Reduce(intersect,list(Mom_AMA_DEGs_all,DNA_Meth_genes))#"GADD45B"
Mom_AMA_DEGs_all[order(Mom_AMA_DEGs_all,decreasing = T)]
Mom_AMA_DEGs_all[grepl("HDAC",Mom_AMA_DEGs_all)]
Mom_AMA_DEGs_all[grepl("HAT",Mom_AMA_DEGs_all)]
Mom_AMA_DEGs_all[grepl("JARID",Mom_AMA_DEGs_all)]
Mom_AMA_DEGs_all[grepl("KDM",Mom_AMA_DEGs_all)]
Mom_AMA_DEGs_all[grepl("JMJD",Mom_AMA_DEGs_all)]
Mom_AMA_DEGs_all[grepl("LSD",Mom_AMA_DEGs_all)]
Mom_AMA_DEGs_all[grepl("MYST",Mom_AMA_DEGs_all)]

#read imprinting genes
Reduce(intersect,list(Mom_AMA_DEGs_all,imprinting_genes))
# "IRF7"    "CYP2J2"  "CACNA1E" "SEPT4"   "ZNF714"  "NAV2"    "ZC3H12C" "MDGA1"  
grid.newpage(); #清空画板，开始画新图
venn_Age_imprinting_genes <-venn.diagram(list(imprinted_genes=imprinting_genes,Mom_AMA_DEGs_all=Mom_AMA_DEGs_all),
                                  alpha=c(0.9,0.9),lwd=1,lty=1,col="black" , 
                                  fill=ppCor[3:4],cex = 1.5,cat.col=ppCor[3:4],cat.fontface=4, cat.cex = 1.5,    
                                  main = "imprinted_genes and AMA DEGs",
                                  main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                  filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_mom_RNA_DEGs_venn_Age_imprinting_genes.pdf",width=6,height = 5)
grid.newpage()
grid.draw(venn_Age_imprinting_genes)
dev.off()
#load GenAge database genes
Reduce(intersect,list(Mom_AMA_DEGs_up_1.5,Age_gene_lgh))
#"TNFSF13B" "LMNB1"    "CD55"     "NRG1"     "GRN"      "BCL2A1"   "HOXC4"    "TXN" 
Reduce(intersect,list(Mom_AMA_DEGs_down_1.5,Age_gene_lgh))
# [1] "EGR1"   "NOG"    "IRS1"   "PDGFB"  "CDK1"   "MMP1"   "TOP2A"  "S100B"  "CAMK4"  "A2M"    "EPS8"  
# [12] "FOXM1"  "IGFBP2" 

grid.newpage(); #清空画板，开始画新图
venn_Age_longevity <-venn.diagram(list(Age_gene=Age_gene_lgh,longevity_gene=longevity_gene,Mom_AMA_DEGs_all=Mom_AMA_DEGs_all),
                                  alpha=c(0.9,0.9,0.9),lwd=1,lty=1,col="black" , 
                                  fill=ppCor[c(10,11,9)],cex = 1.5,cat.col=ppCor[c(10,11,9)],cat.fontface=4, cat.cex = 1.5,    
                                  main = "Database:: longevity gene and two list of Age gene",
                                  main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                  filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_mom_RNA_DEGs_venn_Age_longevity_agedatebase_genes.pdf",width=6,height =6)
grid.newpage()
grid.draw(venn_Age_longevity)
dev.off()
venn_Age_DEGs_genAge<-venn.diagram(list(genAge_database=genAge_database_extand,Mom_DEGs_up=Mom_AMA_DEGs_up_1.5,Mom_DEGs_down=Mom_AMA_DEGs_down_1.5),
                                   alpha=c(1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(5,1:2)], cex = 1.5, 
                                   cat.col=ppCor[c(5,1:2)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_mom_RNA_DEGs_venn_Age_longevity_agedatebase_genes2.pdf",width=6,height =5)
grid.newpage()
grid.draw(venn_Age_DEGs_genAge)
dev.off()

## SASP marker genes in DEGs
Reduce(intersect,list(Mom_AMA_DEGs_all,SASP_gene))
# "ETS2"   "CD55"   "GMFG"   "MMP1"   "ITPKA"  "IGFBP2"
Reduce(intersect,list(Mom_AMA_DEGs_up_1.5,SASP_gene))
#[1] "ETS2" "CD55" "GMFG"
Reduce(intersect,list(Mom_AMA_DEGs_down_1.5,SASP_gene))
#[1] "MMP1"   "ITPKA"  "IGFBP2"

venn_Age_DEGs_SASP_gene<-venn.diagram(list(SASP_gene=SASP_gene,Mom_DEGs_all=Mom_AMA_DEGs_all),
                                   alpha=c(1,1),lwd=1,lty=1,  col="black" , fill=ppCor[c(10,6)], 
                                   cex = 1.5, cat.col=ppCor[c(10,6)],cat.fontface=4,cat.cex = 1.5,
                                   main.cex = 2, main.fontface = 2, main.fontfamily = 3,  filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_mom_RNA_DEGs_venn_Age_SASP_gene.pdf",width=6,height =5)
grid.newpage()
grid.draw(venn_Age_DEGs_SASP_gene)
dev.off()
##human TFs and TFs_cofactor 
Reduce(intersect,list(Mom_AMA_DEGs_all,Hm_TFs))
Reduce(intersect,list(Mom_AMA_DEGs_all,Hm_coTFs))
Reduce(intersect,list(Mom_AMA_DEGs_up_1.5,Hm_coTFs))
Reduce(intersect,list(Mom_AMA_DEGs_down_1.5,Hm_coTFs))

venn_Age_DEGs_TFs<-venn.diagram(list(Mom_DEGs_up=Mom_AMA_DEGs_up_1.5,Mom_DEGs_down=Mom_AMA_DEGs_down_1.5,Homo_sapiens_TF=Hm_TFs,Homo_sapiens_coTF=Hm_coTFs),
                                   alpha=c(1,1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1:3,5)],cex = 1.5, cat.col=ppCor[c(1:3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_mom_RNA_DEGs_venn_Age_TFs_gene.pdf",width=8,height =6)
grid.newpage();
grid.draw(venn_Age_DEGs_TFs)
dev.off()

#for RNA expression
head(Mom_draw_count)
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 5),rep("Young", 5))))
rownames(annotation_col) = colnames(Mom_draw_count)
## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))
##展示行或者列的label

Gene_expression<-Mom_draw_count[which(rownames(Mom_draw_count) %in% DNA_Meth_genes),]
heat_plot1<-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                     #  clustering_distance_rows ="euclidean",#correlation
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors, 
                     #gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                     treeheight_col = 20, #treeheight_row = 30, 
                     # labels_row = labels_row,
                     #border_color ="red", 
                     border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     #   filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps1.pdf"),
                     main ="Mom_Gene expression of epigenetic Genes:original value",angle_col ="45")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =T,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors, 
                      # gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                      treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      main ="Mom_Gene expression level of epigenetic Genes: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Mom_RNA_gene_epigenetic_heatmap1.pdf",width=6,height =8)
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Mom_RNA_gene_epigenetic_heatmap2.pdf",width=6,height =8)
print(heat_plot2)
dev.off()

#imprinting_genes
Gene_expression<-Mom_draw_count[which(rownames(Mom_draw_count) %in% imprinting_genes),]
heat_plot1<-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                     #  clustering_distance_rows ="euclidean",#correlation
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors, 
                     #gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                     treeheight_col = 20, #treeheight_row = 30, 
                     # labels_row = labels_row,
                     #border_color ="red", 
                     border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     #   filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps1.pdf"),
                     main ="Mom_Gene expression of selected meth_trans_inter Genes_merge:original value",angle_col ="45")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =T,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors, 
                      # gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                      treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      #  filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps2.pdf"),
                      main ="Mom_Gene expression level of imprinting_genes: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Mom_imprinting_genes_heatmap1.pdf",width=10,height =12)
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Mom_imprinting_genes_heatmap2.pdf",width=10,height =12)
print(heat_plot2)
dev.off()

#SASP_gene
Gene_expression<-Mom_draw_count[which(rownames(Mom_draw_count) %in% SASP_gene),]
heat_plot1<-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                     #  clustering_distance_rows ="euclidean",#correlation
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors, 
                     #gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                     treeheight_col = 20, #treeheight_row = 30, 
                     # labels_row = labels_row,
                     #border_color ="red", 
                     border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     #   filename =paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps1.pdf"),
                     main ="Mom_Gene expression of SASP_gene:original value",angle_col ="45")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =T,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors, 
                      # gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
                      treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      main ="Mom_Gene expression level of SASP_gene: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Mom_SASP_gene_heatmap1.pdf",width=8,height =10)
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Mom_SASP_gene_heatmap2.pdf",width=8,height =10)
print(heat_plot2)
dev.off()

#For overlapped genes
Mom_kid_overlap_DEGs<-Reduce(intersect,list(Mom_AMA_DEGs_all,kids_AMA_DEGs_all))
common_UP_DEGs<-Reduce(intersect,list(kids_AMA_DEGs_up_1.5,Mom_AMA_DEGs_up_1.5))
common_Down_DEGs<-Reduce(intersect,list(kids_AMA_DEGs_down_1.5,Mom_AMA_DEGs_down_1.5))
common_DEGs<-c(common_UP_DEGs,common_Down_DEGs)
#in  Mom
Reduce(intersect,list(common_DEGs,DNA_Meth_genes))#"GADD45B"
common_DEGs[order(common_DEGs,decreasing = T)]
common_DEGs[grepl("HDAC",common_DEGs)]
common_DEGs[grepl("HAT",common_DEGs)]
common_DEGs[grepl("JARID",common_DEGs)]
common_DEGs[grepl("KDM",common_DEGs)]
common_DEGs[grepl("JMJD",common_DEGs)]
common_DEGs[grepl("LSD",common_DEGs)]
common_DEGs[grepl("MYST",common_DEGs)]

#read imprinting genes
Reduce(intersect,list(common_DEGs,imprinting_genes))

#load GenAge database genes
Reduce(intersect,list(common_DEGs,Age_gene_lgh))
#"NRG1"   "BCL2A1" "FOXM1"  "EGR1" 
Reduce(intersect,list(common_DEGs,longevity_gene))
#[1] "MS4A6A"  "AIF1"    "RAB44"   "RPS6KA2" "PLXNA4"

grid.newpage(); #清空画板，开始画新图
venn_Age_longevity1 <-venn.diagram(list(Age_gene=Age_gene_lgh,longevity_gene=longevity_gene,common_DEGs=common_DEGs),
                                  alpha=c(0.9,0.9,0.9),lwd=1,lty=1,col="black" , 
                                  fill=ppCor[c(10,11,9)],cex = 1.5,cat.col=ppCor[c(10,11,9)],cat.fontface=4, cat.cex = 1.5,    
                                  main = "Database:: longevity gene and two list of Age gene",
                                  main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                  filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_common_RNA_DEGs_longevity_Aged_gene.pdf",width=8,height =8)
grid.newpage()
grid.draw(venn_Age_longevity1)
dev.off()


venn_Age_DEGs_genAge0<-venn.diagram(list(genAge_database=genAge_database_extand,common_UP_DEGs=common_UP_DEGs,common_Down_DEGs=common_Down_DEGs),
                                   alpha=c(1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(5,1:2)], cex = 1.5, 
                                   cat.col=ppCor[c(5,1:2)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_common_RNA_DEGs_longevity_Aged_gene1.pdf",width=6,height =6)
grid.newpage()
grid.draw(venn_Age_DEGs_genAge0)
dev.off()


venn_Age_DEGs_genAge1<-venn.diagram(list(genAge_database=genAge_database_extand,kids_DEGs_down=kids_AMA_DEGs_down_1.5,kids_DEGs_up=kids_AMA_DEGs_up_1.5,Mom_DEGs_down=Mom_AMA_DEGs_down_1.5,Mom_DEGs_up=Mom_AMA_DEGs_up_1.5),
                                   alpha=c(1,1,1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(5,7,1:3)], cex = 1.5, 
                                   cat.col=ppCor[c(5,7,1:3)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
grid.newpage()
grid.draw(venn_Age_DEGs_genAge1)
venn_Age_DEGs_genAge2<-venn.diagram(list(genAge_database=genAge_database_extand,kids_DEGs=c(kids_AMA_DEGs_down_1.5,kids_AMA_DEGs_up_1.5),Mom_DEGs=c(Mom_AMA_DEGs_down_1.5,Mom_AMA_DEGs_up_1.5)),
                                   alpha=c(0.7,0.7,0.7),lwd=1,lty=1,col="black" ,fill=ppCor[c(2:4)], cex = 1.5, 
                                   cat.col=ppCor[c(2:4)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
grid.newpage()
grid.draw(venn_Age_DEGs_genAge2)

pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_common_RNA_DEGs_longevity_Aged_gene3.pdf",width=6,height =6)
grid.draw(venn_Age_DEGs_genAge1)
grid.newpage()
grid.draw(venn_Age_DEGs_genAge2)
dev.off()

## SASP marker genes in DEGs
Reduce(intersect,list(common_DEGs,SASP_gene))
#"GMFG"

venn_Age_DEGs_SASP_gene<-venn.diagram(list(SASP_gene=SASP_gene,common_DEGs=common_DEGs),
                                   alpha=c(1,1),lwd=1,lty=1,  col="black" , fill=ppCor[c(10,6)], 
                                   cex = 1.5, cat.col=ppCor[c(10,6)],cat.fontface=4,cat.cex = 1.5,
                                   main.cex = 2, main.fontface = 2, main.fontfamily = 3,  filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_common_RNA_DEGs_SASP_gene.pdf",width=8,height =6)
grid.newpage()
grid.draw(venn_Age_DEGs_SASP_gene)
dev.off()
##human TFs and TFs_cofactor 
Reduce(intersect,list(common_DEGs,Hm_TFs))
#"FOXM1"  "LTF"    "EGR1"   "ERG"    "CEBPE"  "E2F8"   "ZNF471" "TWIST2"
Reduce(intersect,list(common_DEGs,Hm_coTFs))
#"NRG1"   "BUD31"  "NEDD8"  "MED12L" "PLK1"   "ELANE" 

venn_Age_DEGs_TFs<-venn.diagram(list(common_DEGs=common_DEGs,Homo_sapiens_TF=Hm_TFs,Homo_sapiens_coTF=Hm_coTFs),
                                   alpha=c(1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1,3,5)],cex = 1.5, cat.col=ppCor[c(1,3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)

pdf("/mnt/data/chenwei/huahua/5.valificat_result/venn_common_RNA_DEGs_TFs_gene.pdf",width=8,height =6)
grid.newpage();
grid.draw(venn_Age_DEGs_genAge)
dev.off()

#for RNA expression  需要时候进行
#读入表达数据矩阵normal count or log2
target_gene<-c("GADD45B","DNMT1","ETS2")
target_gene<-c( "CEACAM8","DDX11","HTRA3","IFITM1","LPO","MED12L","MS4A6A" )
target_gene<-c("GMFG")
kids_data<-kids_draw_count[target_gene,]
Mom_data<-Mom_draw_count[target_gene,]
kids_mom_data<-merge(kids_data,Mom_data,by=0)
rownames(kids_mom_data)<-kids_mom_data$Row.names;kids_mom_data<-kids_mom_data[,-1]

## coldata
colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/file3_all_no_filter_verification_AMA_analysis_metadata.csv",row.names=1,header =T)
colData_whole<-colData_whole[which(colData_whole$sample %in% colnames(kids_mom_data)),]
head(colData_whole);dim(colData_whole)

#整合分组信息至表达数据:
data_expr<- as.data.frame(t(kids_mom_data))
data_preplot=merge(data_expr,colData_whole,by=0)
head(data_preplot)
data_plot <- melt(data_preplot,variable.name="Genes",value.name = "expression",id.vars = c("Row.names","sample","group","generation", "group_type","type"))
head(data_plot)

#挑选绘制的亚组
grid.newpage(); #清空画板，开始画新图
colnames(data_plot)
# summarySE 计算标准差和标准误差以及95%的置信区间.
tgc1 <- summarySE(data_plot, measurevar="expression", groupvars= c("group_type","Genes"))
head(tgc1)
colnames(tgc1)<-c("subtype","Genes","N","expression","sd","se","ci")
tgc1$subtype<-factor(tgc1$subtype,levels=c( "Young_Kid","AMA_Kid","Young_Mom","AMA_Mom"))
genetags<-tgc1[which(tgc1$subtype=="Young_Kid"),]
maintile<-"mRNA level of AMA_associated mom and kids common Genes"

ggplot(tgc1, aes(x=subtype, y=expression, colour=Genes,group=Genes)) + 
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),position=position_dodge(0.2), width=1) +
  geom_line(position=position_dodge(0.2)) +
  geom_point(size=2,position=position_dodge(0.2))+
  ggtitle(maintile) + xlab("Subtype") +ylab("log2(expression+1)") +
  expand_limits(x=-.5) +                        # Expand y range
  scale_y_continuous(breaks=0:30*1) +# Set tick every 4
  theme_bw() + #  theme(legend.position = 'top')+   ##同理可以取 bottom、right、left
  theme(plot.title = element_text(hjust=0.5,vjust=0.5,size=20,colour = "black",face = "bold"))+#face="bold"加粗
  geom_text(data=genetags,aes(x=0,y=expression,color=Genes,label=genetags$Genes),size=4, alpha=.8,hjust=0)

#优化：https://blog.csdn.net/g_r_c/article/details/19673625
#gene_target4<- bitr(gene_target3, fromType ="REFSEQ", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
ggplot(tgc1, aes(x=subtype, y=expression, colour=Genes,group=Genes)) + 
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=0.05) +
  geom_line() +
  geom_point(size=2)+#控制两线各向左向右移0.2的聚类，线也移，点也移，永远不分开。
  ggtitle(maintile) +
  xlab("Subtype") +
  ylab("expression") +
  expand_limits(x=-.5) +                        # Expand y range
  scale_y_continuous(breaks=0:30*1) +# Set tick every 4
  theme_bw() +
  #  theme(legend.position = 'top')+   ##同理可以取 bottom、right、left
  theme(plot.title = element_text(hjust=0.5,vjust=0.5,size=20,colour = "black",face = "bold"))+#face="bold"加粗
  geom_text(data=genetags,aes(x=0,y=expression,color=Genes,label=genetags$Genes),size=4, alpha=.5,hjust=0)


###plot for seperated genes
unique(data_plot$group_type)
data_plot$group_type<-factor(data_plot$group_type, levels=c( "Young_Kid","AMA_Kid","Young_Mom","AMA_Mom"),ordered=TRUE)
data_plot<-data_plot[order(data_plot$group_type),]
unique(data_plot$group_type)

data_plot_mean <- ddply(data_plot,c("group_type","Genes"),summarise,mean = mean(expression,na.rm =T),length = length(!is.na(expression)))
head(data_plot_mean);range(data_plot_mean$length)## 6 8

###个体分开绘制
my_comparisons <- list(c("Young_Kid","AMA_Kid"),c("Young_Mom","AMA_Mom"))
head(data_plot)
P_mean_site<-ggboxplot(data_plot, x = "group_type", y = "expression", color = "group_type", palette = ppCor, add = "jitter")+#ylim(c(70,85))+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),
        axis.title.y = element_text(size=10,colour = "black",face = "bold"),
        axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))+
  ylim(0,max(data_plot$expression)+1)+facet_wrap(~ Genes, scales = "free",ncol=4)
P_mean_site2<-P_mean_site+stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format",label.y=max(data_plot$expression)+0.5) +labs(title="mRNA expression level for AMA related genes(wilcox.test)", x ="sample_group", y ="log2(count+1)")
P_mean_site3<-P_mean_site+stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.format",label.y=max(data_plot$expression)+0.5) + labs(title="mRNA expression level for AMA related genes(t.test)", x ="sample_group", y ="log2(count+1)")

gene_names<-"GMFG"
ggsave(P_mean_site,file=paste0("/mnt/data/chenwei/huahua/5.valificat_result/box_plot_trans_level_of_siginificant_AMA_common_DEGs_between_kids_and_mother_seprated_",gene_names,".pdf"),width = 8, height =6)

##intergrated with AMA DEGs identified in human oocyte
## AMA oocyte DEGs
hO_AMA_DEGs <- read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/human_oocyte_AMA_related_DEGs_extend.csv")
head(hO_AMA_DEGs)
hO_AMA_DEGs2<-hO_AMA_DEGs[which(hO_AMA_DEGs$p_val<0.05 & abs(hO_AMA_DEGs$avg_diff) > log2(1.5)),]
hO_AMA_DEGs_all<-unique(as.character(hO_AMA_DEGs2$gene))
hO_AMA_DEGs_Up<-unique(as.character(hO_AMA_DEGs2[which(hO_AMA_DEGs2$avg_diff > 0),]$gene))
hO_AMA_DEGs_Down<-unique(as.character(hO_AMA_DEGs2[which(hO_AMA_DEGs2$avg_diff < 0),]$gene))
length(hO_AMA_DEGs_all);length(hO_AMA_DEGs_Up);length(hO_AMA_DEGs_Down)#1157 964  193

length(Reduce(intersect,list(kids_AMA_DEGs_all,hO_AMA_DEGs_all)))# 214
length(Reduce(intersect,list(Mom_AMA_DEGs_all,hO_AMA_DEGs_all)))#30
length(Reduce(intersect,list(common_DEGs,hO_AMA_DEGs_all)))# 4
Reduce(intersect,list(common_DEGs,hO_AMA_DEGs_all))
# "PSENEN" "PLK1"   "SCD"    "CD24"

venn_Age_DEGs<-venn.diagram(list(Mom_DEGs_up=Mom_AMA_DEGs_up_1.5,Mom_DEGs_down=Mom_AMA_DEGs_down_1.5,hO_AMA_DEGs_Up=hO_AMA_DEGs_Up,hO_AMA_DEGs_Down=hO_AMA_DEGs_Down),
                                   alpha=c(1,1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1:3,5)],cex = 1.5, cat.col=ppCor[c(1:3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
grid.newpage()
grid.draw(venn_Age_DEGs)
venn_Age_DEGs<-venn.diagram(list(kids_DEGs_up=kids_AMA_DEGs_up_1.5,kids_DEGs_down=kids_AMA_DEGs_down_1.5,hO_AMA_DEGs_Up=hO_AMA_DEGs_Up,hO_AMA_DEGs_Down=hO_AMA_DEGs_Down),
                            alpha=c(1,1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1:3,5)],cex = 1.5, cat.col=ppCor[c(1:3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                            cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                            filename = NULL)
grid.newpage()
grid.draw(venn_Age_DEGs)
venn_Age_DEGs<-venn.diagram(list(common_UP_DEGs=common_UP_DEGs,common_Down_DEGs=common_Down_DEGs,hO_AMA_DEGs_Up=hO_AMA_DEGs_Up,hO_AMA_DEGs_Down=hO_AMA_DEGs_Down),
                            alpha=c(1,1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1:3,5)],cex = 1.5, cat.col=ppCor[c(1:3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                            cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                            filename = NULL)
grid.newpage()
grid.draw(venn_Age_DEGs)


#upset plot for all six lists of DEGs
listinput<-c(list(kids_AMA_DEGs_up_1.5),list(kids_AMA_DEGs_down_1.5),list(Mom_AMA_DEGs_up_1.5),list(Mom_AMA_DEGs_down_1.5),list(hO_AMA_DEGs_Up),list(hO_AMA_DEGs_Down))
names(listinput)<-c("kids_DEGs_up","kids_DEGs_down","Mom_DEGs_up","Mom_DEGs_down","hO_AMA_DEGs_Up","hO_AMA_DEGs_Down")
upset_dataframe<-as.data.frame(fromList(listinput))
dim(upset_dataframe)
upset_dataframe[1:10,1:5]
upset_dataframe_rowSums<-rowSums(upset_dataframe)
range(upset_dataframe_rowSums)#1 3
upset_dataframe_colSums<-colSums(upset_dataframe)
range(upset_dataframe_colSums)#193 2029
upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]
upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400)

list2_q1<-list(query = intersects, params = list("kids_DEGs_up","Mom_DEGs_up","hO_AMA_DEGs_Up"), color = ppCor_all[1], active = T,query.name = "Common Up DEGs")
list2_q2<-list(query = intersects, params = list("kids_DEGs_down","Mom_DEGs_down","hO_AMA_DEGs_Down"), color = ppCor_all[1], active = T,query.name = "Common Down DEGs")
list2_q3<-list(query = intersects, params = list("kids_DEGs_up","hO_AMA_DEGs_Up"), color = ppCor_all[2], active = T,query.name = "Kids_common_hyper DMRs")
list2_q4<-list(query = intersects, params = list("kids_DEGs_down","hO_AMA_DEGs_Down"), color = ppCor_all[2], active = T,query.name = "mom_common_hyper DMRs")
list2_q5<-list(query = intersects, params = list("Mom_DEGs_up","hO_AMA_DEGs_Up"), color = ppCor_all[3], active = T,query.name = "pre_common hypo DMRs")
#list2_q6<-list(query = intersects, params = list("Mom_DEGs_down","hO_AMA_DEGs_Down"), color = ppCor_all[3], active = T,query.name = "val_common hypo DMRs")
list2_q7<-list(query = intersects, params = list("kids_DEGs_up","Mom_DEGs_up"), color = ppCor_all[5], active = T,query.name = "Kids_common_hypo DMRs")
list2_q8<-list(query = intersects, params = list("kids_DEGs_down","Mom_DEGs_down"), color = ppCor_all[5], active = T,query.name = "mom_common_hypo DMRs")

upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), #nintersects = 400,
      sets=c("kids_DEGs_up","Mom_DEGs_up","hO_AMA_DEGs_Up","kids_DEGs_down","Mom_DEGs_down","hO_AMA_DEGs_Down"),
      keep.order = TRUE,
      query.legend = "top",
      sets.bar.color = "brown",
      #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
      shade.color="pink",
      #      matrix.color="purple",
      order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
      point.size = 3,line.size = 1.3,
      mainbar.y.label = "DMRs number Intersections", sets.x.label = "DMRs number per subset",
      mb.ratio = c(0.60, 0.40),
      #text.scale = c(1,1,1,1,1,1,1,1),
      show.numbers = 'yes',
      queries = list(list2_q1,list2_q2,list2_q3,list2_q4,list2_q5,list2_q7,list2_q8)
)

hO_AMA_DEGs[which(hO_AMA_DEGs$gene %in% c("PSENEN","SCD","CD24")),]

setdiff(Reduce(intersect,list(common_UP_DEGs,hO_AMA_DEGs_Up)),common_Down_DEGs)
setdiff(Reduce(intersect,list(common_Down_DEGs,hO_AMA_DEGs_Down)),common_UP_DEGs)
setdiff(Reduce(intersect,list(Mom_AMA_DEGs_up_1.5,hO_AMA_DEGs_Up)),Mom_AMA_DEGs_down_1.5)
setdiff(Reduce(intersect,list(Mom_AMA_DEGs_down_1.5,hO_AMA_DEGs_Down)),Mom_AMA_DEGs_up_1.5)
setdiff(Reduce(intersect,list(kids_AMA_DEGs_up_1.5,hO_AMA_DEGs_Up)),kids_AMA_DEGs_down_1.5)
setdiff(Reduce(intersect,list(kids_AMA_DEGs_down_1.5,hO_AMA_DEGs_Down)),kids_AMA_DEGs_up_1.5)