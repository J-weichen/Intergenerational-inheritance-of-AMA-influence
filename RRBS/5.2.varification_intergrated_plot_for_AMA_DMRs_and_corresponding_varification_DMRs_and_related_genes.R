rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(scales)
library(ggsci)
library(VennDiagram)

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

#for kids 
##meth file
kids_meth_data<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_kids_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth_PVALUE.csv",header=T,check.names = F)
kids_meth_data$Row.names <- unlist(lapply(strsplit(as.character(kids_meth_data$Row.names), "[.]"), function(x) paste(x[1],x[2],x[3],sep="_")))

## expression file
tag_name<-"Kid"; sampleA <-"AMA";sampleB <-"Young" 
kids_exp_file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file9_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
kids_log2_normalized_Mat<-read.table(kids_exp_file,sep="\t", header=T)
head(kids_log2_normalized_Mat);head(kids_meth_data)

#甲基化对应基因列表的读取
#data preparation
compare_name1<-"kids_AMA_vs_Young";split_region<- "200bp"
kids_data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t",check.names = F)
DMR_hyper<-kids_data_region_merge_DMR_all_big[kids_data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-kids_data_region_merge_DMR_all_big[kids_data_region_merge_DMR_all_big$meth.diff <0,]

kids_hyper_gene<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
kids_hypo_gene<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))
kids_DMR_hyper<-unique(as.character(DMR_hyper$Row.names))
kids_DMR_hypo<-unique(as.character(DMR_hypo$Row.names))

dim(kids_data_region_merge_DMR_all_big)#423  43
length(kids_DMR_hyper);length(kids_DMR_hypo)#186 #237
length(kids_hyper_gene);length(kids_hypo_gene)#180 #211

#绘韦恩图
venn1 <-venn.diagram(list(verif_DMRs_gene_names=unique(as.character(kids_meth_data$Row.names)),
                          kids_DMR_hyper=kids_DMR_hyper, kids_DMR_hypo= kids_DMR_hypo),
                     alpha=c(0.7,0.7,0.7),lwd=1,lty=1,col="white",fill=ppCor[c(2:4)], 
                     cex = 1.5,cat.col=ppCor[c(2:4)], 
                     cat.fontface=4, cat.cex = 1.5, main=paste0("Kids:","varification bins and pre_DMRs"), 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)
grid.newpage(); 
grid.draw(venn1)

venn1 <-venn.diagram(list(RNA_gene_names=unique(as.character(kids_log2_normalized_Mat$Row.names)),
                          kids_hyper_gene=kids_hyper_gene, kids_hypo_gene= kids_hypo_gene),
                     alpha=c(0.7,0.7,0.7),lwd=1,lty=1,col="white",fill=c("#20854EFF","#FFDC91FF", "#4DBBD5FF"), 
                     cex = 1.5,cat.col=c("#20854EFF","#FFDC91FF", "#4DBBD5FF"),#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=paste0("Kids:","gene and pre_DMR related gene"), 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)
grid.newpage(); 
grid.draw(venn1)
grid.newpage(); 

#building new dataframe for varification
kids_meth_data_new<-kids_meth_data[which(kids_meth_data$Row.names %in% c(kids_DMR_hyper,kids_DMR_hypo)),]
kids_trans_data_new<-kids_log2_normalized_Mat[which(kids_log2_normalized_Mat$Row.names %in% c(kids_hyper_gene,kids_hypo_gene)),]
dim(kids_meth_data_new);dim(kids_trans_data_new)# 315  17 #240  20
table(kids_meth_data_new$pvalue<0.05)
#FALSE  TRUE 
# 132    183 
table(kids_trans_data_new$pvalue<0.05)
#FALSE  TRUE 
# 226    14 
kids_meth_hyper_new<-kids_meth_data[which(kids_meth_data$Row.names %in% c(kids_DMR_hyper)),]
table(kids_meth_hyper_new$pvalue<0.05 & kids_meth_hyper_new$meth.diff >0)
#FALSE  TRUE 
# 101    33 
kids_meth_hypo_new<-kids_meth_data[which(kids_meth_data$Row.names %in% c(kids_DMR_hypo)),]
table(kids_meth_hypo_new$pvalue<0.05 & kids_meth_hypo_new$meth.diff <0)
#FALSE  TRUE 
# 128    53 
kids_meth_varified<-rbind(kids_meth_hyper_new[which(kids_meth_hyper_new$pvalue<0.05 & kids_meth_hyper_new$meth.diff >0),],
                          kids_meth_hypo_new[which(kids_meth_hypo_new$pvalue<0.05 & kids_meth_hypo_new$meth.diff <0),])

head(kids_meth_varified)

kids_meth_varified_genes<-unique(kids_data_region_merge_DMR_all_big[which(kids_data_region_merge_DMR_all_big$Row.names %in% kids_meth_varified$Row.names),]$SYMBOL)
kids_trans_varified_genes<-unique(kids_trans_data_new[which(kids_trans_data_new$pvalue<0.05),]$Row.names)

venn1 <-venn.diagram(list(kids_meth_varified_genes=kids_meth_varified_genes,
                          kids_trans_varified_genes=kids_trans_varified_genes),
                     alpha=c(0.7,0.7),lwd=1,lty=1,col="white",fill=ppCor[c(5,10)], 
                     cex = 1.5,cat.col=ppCor[c(5,10)], 
                     cat.fontface=4, cat.cex = 1.5, main=paste0("Kids:","varification bins and pre_DMRs"), 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)
grid.newpage()
grid.draw(venn1)
grid.newpage(); 

genelist<-intersect(kids_meth_varified_genes,kids_trans_varified_genes)
genelist
#"CD247"    "PPP1R12B" "CITED4"   "CNNM2"    "ITGA2B"   "EMILIN2"  "ZNF407"   "POLRMT"   "FAM95C"  
key_kids_expression<-kids_log2_normalized_Mat[which(kids_log2_normalized_Mat$Row.names %in% genelist),]
key_kids_meth<-kids_data_region_merge_DMR_all_big[which(kids_data_region_merge_DMR_all_big$SYMBOL %in% genelist),]
key_kids_meth_new<-kids_meth_data[which(kids_meth_data$Row.names %in% key_kids_meth$Row.names),]

#绘图
plot_trans<-key_kids_expression
plot_trans$genetype<-as.character(plot_trans$Row.names)
rownames(plot_trans)<-plot_trans$Row.names
plot_trans<-plot_trans[,-1]

head(plot_trans)
#for RNA expression
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 5),rep("Young", 5))))
rownames(annotation_col) = colnames(plot_trans[1:10])
## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))
##展示行或者列的label
Gene_expression<-plot_trans[,1:10]

#pdf(paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps.pdf"))

heat_plot1<-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,annotation_colors = ann_colors, 
                     treeheight_col = 20,  border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     main ="Kids Gene expression level of key changed Genes :original value",angle_col ="45")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col, annotation_colors = ann_colors, treeheight_col = 20, 
                      scale ="row", border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      main ="Kids Gene expression level of key changed Genes: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Kid_RNA_key_changed_gene_in_pre_and_valif_DMRs_heatmaps1.pdf")
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Kid_RNA_key_changed_gene_in_pre_and_valif_DMRs_heatmaps2.pdf")
print(heat_plot2)
dev.off()
#dev.new()
anno<-key_kids_meth[,c("Row.names","SYMBOL")]
key_kids_meth_new2<-merge(key_kids_meth_new,anno,by="Row.names")

plot_DMRs<-key_kids_meth
plot_DMRs$DMRs<-plot_DMRs$Row.names
plot_DMRs$SYMBOL<-as.character(plot_DMRs$SYMBOL)

plot_DMRs$SYMBOL[duplicated(plot_DMRs$SYMBOL)]
#[1] "FAM95C"
plot_DMRs$SYMBOL2<-plot_DMRs$SYMBOL
plot_DMRs[plot_DMRs$SYMBOL == "FAM95C",]$SYMBOL2<-c("FAM95C(1)","FAM95C(2)")
rownames(plot_DMRs)<-plot_DMRs$SYMBOL2

sample_list<-c("10-3-1E","16-3-4D","17-3-2F","18-3","19-3-4E","20Q-4H","6-3","7-3","8-3","9-3-1C",
               "11-3-3F","12-3-4A","13-3-1G","14-3","15-3-2C","1Q","2-3-3C","3-3","4-3","5-3")
heat_plot<-plot_DMRs[,c(sample_list,"SYMBOL2","SYMBOL")]

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
rownames(huahua_meta)<-huahua_meta$library_code
huahua_meta2<-huahua_meta[sample_list,]
head(huahua_meta2)
colnames(heat_plot)<-c(huahua_meta2$analysis_name,"SYMBOL2","SYMBOL")
as.character(huahua_meta2$analysis_name[order(huahua_meta2$analysis_name,decreasing = F)])
sample_list2 <-c(paste0("AMA_K_",1:10),paste0("YOUNG_K_",1:10))
plot_DMRs2<-heat_plot[,c(sample_list2,"SYMBOL2","SYMBOL")]
head(plot_DMRs2)

#for DNA methylation
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 10),rep("Young", 10))))
rownames(annotation_col) =sample_list2

## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))

grid.newpage(); 
heat_plot3 <-pheatmap(plot_DMRs2[,sample_list2], cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_colors = ann_colors, 
                      treeheight_col = 20, border=FALSE,
                      color = colorRampPalette(c("navy","pink","orange","firebrick3"))(50),
                      main ="original Kids DNA methylation level of pre-DMRs in key changed Genes: original value",angle_col ="90")

#scale by row
heat_plot4 <-pheatmap(plot_DMRs2[,sample_list2], cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col, annotation_colors = ann_colors, 
                      treeheight_col = 20, scale ="row", border=FALSE,
                      color = colorRampPalette(c("purple","white","orange"))(40),
                      main ="Kids DNA methylation level of pre-DMRs in key changed Genes:scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Kid_meth_pre_DMRs_in_key_changed_gene_heatmaps1.pdf")
print(heat_plot3)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Kid_meth_pre_DMRs_in_key_changed_gene_heatmaps2.pdf")
print(heat_plot4)
dev.off()

##for validation
plot_DMRs<-key_kids_meth_new2
plot_DMRs$DMRs<-plot_DMRs$Row.names
plot_DMRs$SYMBOL<-as.character(plot_DMRs$SYMBOL)

plot_DMRs$SYMBOL[duplicated(plot_DMRs$SYMBOL)]
#[1] "FAM95C"
plot_DMRs$SYMBOL2<-plot_DMRs$SYMBOL
plot_DMRs[plot_DMRs$SYMBOL == "FAM95C",]$SYMBOL2<-c("FAM95C(1)","FAM95C(2)")
rownames(plot_DMRs)<-plot_DMRs$SYMBOL2
head(plot_DMRs)

#heat_plot<-plot_DMRs[,c(2:11)]
sample_list<-c("A2Q","A3Q","A4Q","A5Q","A6Q","Y1Q","Y2Q","Y4Q","Y5Q","Y6Q")
heat_plot<-plot_DMRs[,sample_list]
#for DNA methylation
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 5),rep("Young", 5))))
rownames(annotation_col) =sample_list

## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))

grid.newpage(); 
heat_plot3 <-pheatmap(heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_colors = ann_colors, 
                      treeheight_col = 20, border=FALSE,
                      color = colorRampPalette(c("navy","pink","orange","firebrick3"))(50),
                      main ="varification Kids DNA methylation level of  of key changed Genes: original value",angle_col ="90")

#scale by row
heat_plot4 <-pheatmap(heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col, annotation_colors = ann_colors, 
                      treeheight_col = 20, scale ="row", border=FALSE,
                      color = colorRampPalette(c("purple","white","orange"))(40),
                      main ="varification Kids DNA methylation level of  of key changed Genes:scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Kid_meth_varif_DMRs_in_key_changed_gene_heatmaps1.pdf")
print(heat_plot3)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Kid_meth_varif_DMRs_in_key_changed_gene_heatmaps2.pdf")
print(heat_plot4)
dev.off()

#for Mom 
##meth file
Mom_meth_data<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/varification_mother_10_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_4sample_meth_PVALUE.csv",header=T,check.names = F)
Mom_meth_data$Row.names <- unlist(lapply(strsplit(as.character(Mom_meth_data$Row.names), "[.]"), function(x) paste(x[1],x[2],x[3],sep="_")))

## expression file
tag_name<-"Mom"; sampleA <-"AMA";sampleB <-"Young" 
Mom_exp_file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file9_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
Mom_log2_normalized_Mat<-read.table(Mom_exp_file,sep="\t", header=T)

head(Mom_log2_normalized_Mat);head(Mom_meth_data)
#甲基化对应基因列表的读取
compare_name1<-"mother_AMA_vs_Young";split_region<- "200bp"
Mom_data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t",check.names = F)
DMR_hyper<-Mom_data_region_merge_DMR_all_big[Mom_data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-Mom_data_region_merge_DMR_all_big[Mom_data_region_merge_DMR_all_big$meth.diff <0,]

Mom_hyper_gene<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
Mom_hypo_gene<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))
Mom_DMR_hyper<-unique(as.character(DMR_hyper$Row.names))
Mom_DMR_hypo<-unique(as.character(DMR_hypo$Row.names))

dim(Mom_data_region_merge_DMR_all_big)#519  41
length(Mom_DMR_hyper);length(Mom_DMR_hypo)#250 #269
length(Mom_hyper_gene);length(Mom_hypo_gene)#223 #251

#绘韦恩图
grid.newpage(); 
venn1 <-venn.diagram(list(verif_DMRs_gene_names=unique(as.character(Mom_meth_data$Row.names)),
                          Mom_DMR_hyper=Mom_DMR_hyper, Mom_DMR_hypo= Mom_DMR_hypo),
                     alpha=c(0.7,0.7,0.7),lwd=1,lty=1,col="white",fill=ppCor[c(2:4)], 
                     cex = 1.5,cat.col=ppCor[c(2:4)], 
                     cat.fontface=4, cat.cex = 1.5, main=paste0("Mom:","varification bins and pre_DMRs"), 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)

grid.draw(venn1)
grid.newpage(); 
venn1 <-venn.diagram(list(RNA_gene_names=unique(as.character(Mom_log2_normalized_Mat$Row.names)),
                          Mom_hyper_gene=Mom_hyper_gene, Mom_hypo_gene= Mom_hypo_gene),
                     alpha=c(0.7,0.7,0.7),lwd=1,lty=1,col="white",fill=c("#20854EFF","#FFDC91FF", "#4DBBD5FF"), 
                     cex = 1.5,cat.col=c("#20854EFF","#FFDC91FF", "#4DBBD5FF"),#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=paste0("Mom:","gene and pre_DMR related gene"), 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)

grid.draw(venn1)
grid.newpage(); 

#building new dataframe
Mom_meth_data_new<-Mom_meth_data[which(Mom_meth_data$Row.names %in% c(Mom_DMR_hyper,Mom_DMR_hypo)),]
Mom_trans_data_new<-Mom_log2_normalized_Mat[which(Mom_log2_normalized_Mat$Row.names %in% c(Mom_hyper_gene,Mom_hypo_gene)),]
dim(Mom_meth_data_new);dim(Mom_trans_data_new)# 350  17 #313  20
table(Mom_meth_data_new$pvalue<0.05)
#FALSE  TRUE 
# 121   229 
table(Mom_trans_data_new$pvalue<0.05)
#FALSE  TRUE 
# 304     9
Mom_meth_hyper_new<-Mom_meth_data[which(Mom_meth_data$Row.names %in% c(Mom_DMR_hyper)),]
table(Mom_meth_hyper_new$pvalue<0.05 & Mom_meth_hyper_new$meth.diff >0)
#FALSE  TRUE 
# 109    52 
Mom_meth_hypo_new<-Mom_meth_data[which(Mom_meth_data$Row.names %in% c(Mom_DMR_hypo)),]
table(Mom_meth_hypo_new$pvalue<0.05 & Mom_meth_hypo_new$meth.diff <0)
#FALSE  TRUE 
# 123    66 
Mom_meth_varified<-rbind(Mom_meth_hyper_new[which(Mom_meth_hyper_new$pvalue<0.05 & Mom_meth_hyper_new$meth.diff >0),],
                         Mom_meth_hypo_new[which(Mom_meth_hypo_new$pvalue<0.05 & Mom_meth_hypo_new$meth.diff <0),])

head(Mom_meth_varified)

Mom_meth_varified_genes<-unique(Mom_data_region_merge_DMR_all_big[which(Mom_data_region_merge_DMR_all_big$Row.names %in% Mom_meth_varified$Row.names),]$SYMBOL)
Mom_trans_varified_genes<-unique(Mom_trans_data_new[which(Mom_trans_data_new$pvalue<0.05),]$Row.names)

grid.newpage(); 
venn1 <-venn.diagram(list(Mom_meth_varified_genes=Mom_meth_varified_genes,
                          Mom_trans_varified_genes=Mom_trans_varified_genes),
                     alpha=c(0.7,0.7),lwd=1,lty=1,col="white",fill=ppCor[c(5,10)], 
                     cex = 1.5,cat.col=ppCor[c(5,10)], 
                     cat.fontface=4, cat.cex = 1.5, main=paste0("Mom:","varification bins and pre_DMRs"), 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)

grid.draw(venn1)
grid.newpage(); 


genelist<-intersect(Mom_meth_varified_genes,Mom_trans_varified_genes)
genelist
#"SLC28A3"
key_Mom_expression<-Mom_log2_normalized_Mat[which(Mom_log2_normalized_Mat$Row.names %in% genelist),]
key_Mom_meth<-Mom_data_region_merge_DMR_all_big[which(Mom_data_region_merge_DMR_all_big$SYMBOL %in% genelist),]
key_Mom_meth_new<-Mom_meth_data[which(Mom_meth_data$Row.names %in% key_Mom_meth$Row.names),]

#绘图
plot_trans<-key_Mom_expression
plot_trans$genetype<-as.character(plot_trans$Row.names)
rownames(plot_trans)<-plot_trans$Row.names
plot_trans<-plot_trans[,-1]

head(plot_trans)
#for RNA expression
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 5),rep("Young", 5))))
rownames(annotation_col) = colnames(plot_trans[1:10])
## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))
##展示行或者列的label
Gene_expression<-plot_trans[,1:10]

#pdf(paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps.pdf"))

heat_plot1<-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,annotation_colors = ann_colors, 
                     treeheight_col = 20,  border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     main ="Mom Gene expression level of key changed Genes :original value",angle_col ="45")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col, annotation_colors = ann_colors, treeheight_col = 20, 
                      scale ="row", border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      main ="Mom Gene expression level of key changed Genes: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Mom_RNA_level_changed_in_both_pre_and_varifi_DMRs_heatmaps1.pdf")
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Mom_RNA_level_changed_in_both_pre_and_varifi_DMRs_heatmaps2.pdf")
print(heat_plot2)
dev.off()
#dev.new()
anno<-key_Mom_meth[,c("Row.names","SYMBOL")]
key_Mom_meth_new2<-merge(key_Mom_meth_new,anno,by="Row.names")

plot_DMRs<-key_Mom_meth
plot_DMRs$DMRs<-plot_DMRs$Row.names
plot_DMRs$SYMBOL<-as.character(plot_DMRs$SYMBOL)

plot_DMRs$SYMBOL[duplicated(plot_DMRs$SYMBOL)]
plot_DMRs$SYMBOL2<-plot_DMRs$SYMBOL
#plot_DMRs[plot_DMRs$SYMBOL == "FAM95C",]$SYMBOL2<-c("FAM95C(1)","FAM95C(2)")
rownames(plot_DMRs)<-plot_DMRs$SYMBOL2

sample_list<-c("10-1-in-9","16-1","17-1-2D","18-1","19-1-2G","20M-4F","6-1","7-1","8-1","9-1-1A",
               "11-1","12-1-3G","13-1-1F","14-1","15-1","1M","2-1","3M","4-1","5-1")
heat_plot<-plot_DMRs[,c(sample_list,"SYMBOL2","SYMBOL")]

huahua_meta2<-huahua_meta[sample_list,]
head(huahua_meta2)
colnames(heat_plot)<-c(huahua_meta2$analysis_name,"SYMBOL2","SYMBOL")
as.character(huahua_meta2$analysis_name[order(huahua_meta2$analysis_name,decreasing = F)])

sample_list2 <-c(paste0("AMA_M_",1:10),paste0("YOUNG_M_",1:10))
plot_DMRs2<-heat_plot[,c(sample_list2,"SYMBOL2","SYMBOL")]
head(plot_DMRs2)

#for DNA methylation
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 10),rep("Young", 10))))
rownames(annotation_col) =sample_list2

## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))

grid.newpage(); 
heat_plot3 <-pheatmap(plot_DMRs2[,sample_list2], cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_colors = ann_colors, 
                      treeheight_col = 20, border=FALSE,
                      color = colorRampPalette(c("navy","pink","orange","firebrick3"))(50),
                      main ="original Mom DNA methylation level of  pre_DMRs of key genes: original value",angle_col ="90")

#scale by row
heat_plot4 <-pheatmap(plot_DMRs2[,sample_list2], cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col, annotation_colors = ann_colors, 
                      treeheight_col = 20, scale ="row", border=FALSE,
                      color = colorRampPalette(c("purple","white","orange"))(40),
                      main ="Mom DNA methylation level of  pre_DMRs of key genes:scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Mom_meth_level_changed_in_pre_DMRs_of_key_gene_heatmaps1.pdf")
print(heat_plot3)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Mom_meth_level_changed_in_pre_DMRs_of_key_gene_heatmaps2.pdf")
print(heat_plot4)
dev.off()

##for validation
plot_DMRs<-key_Mom_meth_new2
plot_DMRs$DMRs<-plot_DMRs$Row.names
plot_DMRs$SYMBOL<-as.character(plot_DMRs$SYMBOL)

plot_DMRs$SYMBOL[duplicated(plot_DMRs$SYMBOL)]
plot_DMRs$SYMBOL2<-plot_DMRs$SYMBOL
rownames(plot_DMRs)<-plot_DMRs$SYMBOL2
head(plot_DMRs)

#heat_plot<-plot_DMRs[,c(2:11)]
sample_list<-c("A2M","A3M","A4M","A5M","A6M","Y1M","Y2M","Y4M","Y5M","Y6M")
heat_plot<-plot_DMRs[,sample_list]
#for DNA methylation
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 5),rep("Young", 5))))
rownames(annotation_col) =sample_list

## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))

grid.newpage(); 
heat_plot3 <-pheatmap(heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_colors = ann_colors, 
                      treeheight_col = 20, border=FALSE,
                      color = colorRampPalette(c("navy","pink","orange","firebrick3"))(50),
                      main ="varification Mom DNA methylation level of varify bin of changed Gene: original value",angle_col ="90")

#scale by row
heat_plot4 <-pheatmap(heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col, annotation_colors = ann_colors, 
                      treeheight_col = 20, scale ="row", border=FALSE,
                      color = colorRampPalette(c("purple","white","orange"))(40),
                      main ="varification Mom DNA methylation level of varify bin of changed Genes:scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Mom_meth_level_changed_in_varify_bin_of_key_gene_heatmaps1.pdf")
print(heat_plot3)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/varification_Mom_meth_level_changed_in_varify_bin_of_key_gene_heatmaps2.pdf")
print(heat_plot4)
dev.off()