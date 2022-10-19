rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(scales)
library(ggsci)
library(pheatmap)
library(VennDiagram)

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

#甲基化对应基因列表的读取
#data preparation
#for first
compare_name1<-"kids_AMA_vs_Young"
split_region<- "200bp"
kids_data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t",check.names = F)
DMR_hyper<-kids_data_region_merge_DMR_all_big[kids_data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-kids_data_region_merge_DMR_all_big[kids_data_region_merge_DMR_all_big$meth.diff <0,]
kids_hyper_gene<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
kids_hypo_gene<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))
kids_DMR_hyper<-unique(as.character(DMR_hyper$Row.names))
kids_DMR_hypo<-unique(as.character(DMR_hypo$Row.names))

dim(kids_data_region_merge_DMR_all_big);length(kids_hyper_gene);length(kids_hypo_gene)
#[1] 451  43
#[1] 170
#[1] 236
#转录组对应基因列表的读取

#read DEGs
tag_name<-"Kid"; sampleA <-"AMA";sampleB <-"Young" 
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file13_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
kids_AMA_DEGs<-read.table(file)
head(kids_AMA_DEGs)
kids_AMA_DEGs_all<-kids_AMA_DEGs$ID
kids_AMA_DEGs_up_1.5<-kids_AMA_DEGs[which(kids_AMA_DEGs$log2FoldChange>0),]$ID
kids_AMA_DEGs_down_1.5<-kids_AMA_DEGs[which(kids_AMA_DEGs$log2FoldChange<0),]$ID
length(kids_AMA_DEGs_up_1.5);length(kids_AMA_DEGs_down_1.5)

#绘韦恩图
tag<- paste0("Kids:"," pre_DMRs versus DEGs")

mainname<-tag
venn1 <-venn.diagram(list(kids_DEGs_up_1.5=kids_AMA_DEGs_up_1.5,kids_DEGs_down_1.5=kids_AMA_DEGs_down_1.5,
                          kids_hypo_gene= kids_hypo_gene ,kids_hyper_gene=kids_hyper_gene),
                     alpha=c(0.9,0.9,0.9,0.9),
                     lwd=1,lty=1,col="white",fill=c("#20854EFF","#FFDC91FF", "#4DBBD5FF", "#F39B7FFF"), 
                     cex = 1.5,cat.col=c("#20854EFF","#FFDC91FF", "#4DBBD5FF", "#F39B7FFF"),#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=mainname, 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)

grid.newpage(); #清空画板，开始画新图
grid.draw(venn1)

nameA<-"kids_Up_gene"
nameB<-"kids_Down_gene"
nameC<-"kids_hypo_DMRs"
nameD<-"kids_hyper_DMRs"

#求具体值
listA<-kids_AMA_DEGs_up_1.5
listB<-kids_AMA_DEGs_down_1.5
listC<-kids_hypo_gene
listD<-kids_hyper_gene

#sink("C:/Users/WeiChen/Desktop/Venn.result.txt",append=FALSE,split=TRUE)

#求交集
genelist<-intersect(listD,listB)
mm<-paste(nameD,"_inter_",nameB,":","  ",genelist,sep="")
mm
#"kids_hyper_DMRs_inter_kids_Down_gene:  
#[1] "PEX14"    "SUFU"     "CNNM2"    "CD82"     "CDC42BPB" "USP36"    "MAP3K10"  "PLA2G4C"  "MED16"    "MYO7B"   
#[11] "PDXK"     "DIP2A"    "LZTR1"    "GRAMD4"   "SYNE1"    "ITPR3"    "GRB10"    "CLIP2"    "MROH6"    "NRBP2"   
#[21] "RALGPS1"  "SLC28A3"

genelist<-intersect(listC,listA)
mm<-paste(nameC,"_inter_",nameA,":","  ",genelist,sep="")
mm 
#"kids_hypo_DMRs_inter_kids_Up_gene: "SLC7A7"    "LINC01550" "RPL23"     "MYCT1"  

genelist<-intersect(listD,listA)
mm<-paste(nameD,"_inter_",nameA,":","  ",genelist,sep="")
mm
#kids_hyper_DMRs_inter_kids_Up_gene:"PRDX1"    "MRPS35"   "MIPEP"    "HLA-DPA1"
genelist<-intersect(listC,listB)
mm<-paste(nameC,"_inter_",nameB,":","  ",genelist,sep="")
mm
#kids_hypo_DMRs_inter_kids_Down_gene
#[1] "ACAP3"    "KAZN"     "PPP1R12B" "CITED4"   "HIVEP3"   "TNFRSF9"  "PFKFB3"   "CTSD"     "SCARB1"   "GLT1D1"  
#[11] "ATP11A"   "ABCC6"    "ITGA2B"   "ITGB3"    "DNAH17"   "ZNF407"   "CTDP1"    "ARHGEF18" "DYSF"     "AATBC"   
#[21] "IQSEC1"   "SLC6A6"   "TTC21A"   "MOSPD3"   "SDK1"     "NCS1"     "VAV2"     "FAM95C"  
###确定目的差异基因关联DMRs位置
plot_gene<-intersect(listC,listA)
plot_gene2<-intersect(listD,listB)
plot_gene1<-intersect(listD,listA)
plot_gene3<-intersect(listC,listB)
inter_genes<-c(plot_gene,plot_gene1,plot_gene2,plot_gene3)
unique(inter_genes)#58

kids_data_region_merge_DMR_all_big[which(kids_data_region_merge_DMR_all_big$SYMBOL %in% inter_genes),c("SYMBOL","annotation","meth.diff")]

kid_exp_file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file9_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
kid_log2_normalized_Mat<-read.table(kid_exp_file,sep="\t", header=T)
head(kid_log2_normalized_Mat);dim(kid_log2_normalized_Mat)#17505    20
range(kid_log2_normalized_Mat[,2:11])# 0.00000 23.06456
plot_trans<-kid_log2_normalized_Mat[which(kid_log2_normalized_Mat$Row.names %in% inter_genes),]
plot_trans$genetype<-as.character(plot_trans$Row.names)
plot_trans[plot_trans$genetype %in% plot_gene, ]$genetype<-"UP_hypo"
plot_trans[plot_trans$genetype %in% plot_gene1, ]$genetype<-"UP_hyper"
plot_trans[plot_trans$genetype %in% plot_gene2, ]$genetype<-"Down_hyper"
plot_trans[plot_trans$genetype %in% plot_gene3, ]$genetype<-"Down_hypo"
rownames(plot_trans)<-plot_trans$Row.names
plot_trans<-plot_trans[,-1]
plot_trans<-plot_trans[inter_genes,]
#plot_trans<-normalized_counts[plot_gene,]
plot_trans$genetype<-factor(plot_trans$genetype,level = unique(as.character(plot_trans$genetype)))
Gene_expression<-plot_trans[order(plot_trans$genetype,decreasing = T),]
head(Gene_expression)
#for RNA expression
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 5),rep("Young", 5))))
rownames(annotation_col) = colnames(Gene_expression[1:10])
#for singlegenes
annotation_row = data.frame(Class = factor(plot_trans$genetype))
rownames(annotation_row) = rownames(plot_trans)

## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(UP_hypo =ppCor[6],UP_hyper =ppCor[3], Down_hyper =ppCor[4],Down_hypo=ppCor[5]))
##展示行或者列的label
Gene_expression<-plot_trans[,1:10]

#pdf(paste0("D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/Intergrated_",compare_name,"_",split_region,"_DEGs_DMRs_heatmaps.pdf"))

heat_plot1<-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                     #  clustering_distance_rows ="euclidean",#correlation
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,annotation_row=annotation_row,
                     annotation_colors = ann_colors, 
                     gaps_row =c(8,12,32),gaps_col =c(5),cutree_col = 2,
                     # treeheight_col = 20, #treeheight_row = 30, 
                     # labels_row = labels_row,
                     #border_color ="red", 
                     border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     main ="Kids Gene expression of selected meth_trans_inter Genes_merge:original value",angle_col ="90")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_row=annotation_row,
                      annotation_colors = ann_colors, 
                      gaps_row =c(8,12,32),gaps_col =c(5),cutree_col = 2,
                      # treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      main ="Kids Gene expression level of meth_trans_inter Genes: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/Kid_RNA_DEGs_pre_DMRs_heatmaps1.pdf")
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/Kid_RNA_DEGs_pre_DMRs_heatmaps2.pdf")
print(heat_plot2)
dev.off()
#dev.new()

head(kids_data_region_merge_DMR_all_big);dim(kids_data_region_merge_DMR_all_big)
kids_data_region_merge_DMR_all_big[which(kids_data_region_merge_DMR_all_big$SYMBOL %in% inter_genes),]

plot_DMRs<-kids_data_region_merge_DMR_all_big[which(kids_data_region_merge_DMR_all_big$SYMBOL %in% inter_genes),]
plot_DMRs$DMRs<-plot_DMRs$Row.names
plot_DMRs$SYMBOL<-as.character(plot_DMRs$SYMBOL)

plot_DMRs$SYMBOL[duplicated(plot_DMRs$SYMBOL)]
#[1] "USP36"  "SDK1"   "IQSEC1" "SLC28A3"
plot_DMRs$SYMBOL2<-plot_DMRs$SYMBOL
plot_DMRs[plot_DMRs$SYMBOL == "SLC28A3",]$SYMBOL2<-c("SLC28A3(1)","SLC28A3(2)")
plot_DMRs[plot_DMRs$SYMBOL == "USP36",]$SYMBOL2<-c("USP36(1)","USP36(2)")
plot_DMRs[plot_DMRs$SYMBOL == "SDK1",]$SYMBOL2<-c("SDK1(1)","SDK1(2)")
plot_DMRs[plot_DMRs$SYMBOL == "IQSEC1",]$SYMBOL2<-c("IQSEC1(1)","IQSEC1(2)")
rownames(plot_DMRs)<-plot_DMRs$SYMBOL2

plot_DMRs[plot_DMRs$SYMBOL %in% plot_gene, ]$SYMBOL<-"UP_hypo"
plot_DMRs[plot_DMRs$SYMBOL %in% plot_gene1, ]$SYMBOL<-"UP_hyper"
plot_DMRs[plot_DMRs$SYMBOL %in% plot_gene2, ]$SYMBOL<-"Down_hyper"
plot_DMRs[plot_DMRs$SYMBOL %in% plot_gene3, ]$SYMBOL<-"Down_hypo"
#ICR_meth_level<-plot_DMRs[1:7]

sample_list<-c("10-3-1E","16-3-4D","17-3-2F","18-3","19-3-4E","20Q-4H","6-3","7-3","8-3","E16C",
               "11-3-3F","12-3-4A","13-3-1G","14-3","15-3-2C","1Q","2-3-3C","3-3","4-3","5-3")
heat_plot<-plot_DMRs[,c(sample_list,"SYMBOL2","SYMBOL")]
head(heat_plot)
huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
huahua_meta$analysis_name<-as.character(huahua_meta$analysis_name)
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

plot_DMRs2$SYMBOL<-factor(plot_DMRs2$SYMBOL,levels = c("UP_hypo","UP_hyper","Down_hyper","Down_hypo"))
plot_DMRs2<-plot_DMRs2[order(plot_DMRs2$SYMBOL),]

annotation_row = data.frame( Class = factor(plot_DMRs2$SYMBOL))
rownames(annotation_row) = rownames(plot_DMRs2)

## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(UP_hypo =ppCor[6],UP_hyper =ppCor[3], Down_hyper =ppCor[4],Down_hypo=ppCor[5]))
gaps_rownumber<-as.numeric(table(plot_DMRs2$SYMBOL))
grid.newpage(); #清空画板，开始画新图
heat_plot3 <-pheatmap(plot_DMRs2[,sample_list2], cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_row=annotation_row,
                      annotation_colors = ann_colors, 
                      gaps_row =c(4,8,32),gaps_col =c(10),cutree_col = 2,
                      #treeheight_col = 20, 
                      #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","pink","orange","firebrick3"))(50),
                      # filename = "D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/DMR_methlevel_genes_over_DEGs_single_heatmap.pdf",
                      main ="Kids DNA methylation level of meth_trans_inter genes original value",angle_col ="90")

#scale by row
heat_plot4 <-pheatmap(plot_DMRs2[,sample_list2], cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_row=annotation_row,
                      annotation_colors = ann_colors, 
                      gaps_row =c(4,8,32),gaps_col =c(10),cutree_col = 2,
                      #cutree_col = 2,
                      #treeheight_col = 20, #treeheight_row = 30, 
                      #labels_row = labels_row,
                      #border_color ="lightgrey", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("purple","white","orange"))(40),
                      #color = colorRampPalette(c("blue","lightgrey","red"))(50),
                      main ="Kids DNA methylation level of selected  DMRs :scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/Kid_meth_DEGs_pre_DMRs_heatmaps1.pdf")
print(heat_plot3)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/Kid_meth_DEGs_pre_DMRs_heatmaps2.pdf")
print(heat_plot4)
dev.off()

#甲基化对应基因列表的读取
#data preparation
#for first
compare_name1<-"mother_AMA_vs_Young"
split_region<- "200bp"
Mom_data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t",check.names = F)

DMR_hyper<-Mom_data_region_merge_DMR_all_big[Mom_data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-Mom_data_region_merge_DMR_all_big[Mom_data_region_merge_DMR_all_big$meth.diff <0,]
Mom_hyper_gene<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
Mom_hypo_gene<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))
Mom_DMR_hyper<-unique(as.character(DMR_hyper$Row.names))
Mom_DMR_hypo<-unique(as.character(DMR_hypo$Row.names))

dim(Mom_data_region_merge_DMR_all_big);length(Mom_hyper_gene);length(Mom_hypo_gene)
#[1] 516  43
#[1] 212
#[1] 253

#mother 转录组对应基因列表的读取
#read DEGs
tag_name<-"Mom"; sampleA <-"AMA";sampleB <-"Young" 
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file13_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
Mom_AMA_DEGs <- read.table(file=file,header=T)
head(Mom_AMA_DEGs)

Mom_AMA_DEGs_all<-Mom_AMA_DEGs$ID
Mom_AMA_DEGs_up_1.5<-Mom_AMA_DEGs[which(Mom_AMA_DEGs$log2FoldChange>0),]$ID
Mom_AMA_DEGs_down_1.5<-Mom_AMA_DEGs[which(Mom_AMA_DEGs$log2FoldChange<0),]$ID
length(Mom_AMA_DEGs_up_1.5);length(Mom_AMA_DEGs_down_1.5)

#绘韦恩图
tag<- paste0("Mom:"," pre_DMRs versus DEGs")

mainname<-tag
venn1 <-venn.diagram(list(Mom_DEGs_up_1.5=Mom_AMA_DEGs_up_1.5,Mom_DEGs_down_1.5=Mom_AMA_DEGs_down_1.5,
                          Mom_hypo_gene= Mom_hypo_gene ,Mom_hyper_gene=Mom_hyper_gene),
                     alpha=c(0.9,0.9,0.9,0.9),
                     lwd=1,lty=1,col="white",fill=c("#20854EFF","#FFDC91FF", "#4DBBD5FF", "#F39B7FFF"), 
                     cex = 1.5,cat.col=c("#20854EFF","#FFDC91FF", "#4DBBD5FF", "#F39B7FFF"),#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=mainname, 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)
grid.newpage();
grid.draw(venn1)

nameA<-"Mom_Up_gene"
nameB<-"Mom_Down_gene"
nameC<-"Mom_hypo_DMRs"
nameD<-"Mom_hyper_DMRs"

#求具体值
listA<-Mom_AMA_DEGs_up_1.5
listB<-Mom_AMA_DEGs_down_1.5
listC<-Mom_hypo_gene
listD<-Mom_hyper_gene

#sink("C:/Users/WeiChen/Desktop/Venn.result.txt",append=FALSE,split=TRUE)

#求交集
genelist<-intersect(listD,listB)
mm<-paste(nameD,"_inter_",nameB,":","  ",genelist,sep="")
mm
#"Mom_hyper_DMRs_inter_Mom_Down_gene :"NAV2"    "ZNF229"  "SLC28A3"
genelist<-intersect(listC,listA)
mm<-paste(nameC,"_inter_",nameA,":","  ",genelist,sep="")
mm
#Mom_hypo_DMRs_inter_Mom_Up_gene:"CD177"  "IL20RB" "AGPAT4"
genelist<-intersect(listD,listA)
mm<-paste(nameD,"_inter_",nameA,":","  ",genelist,sep="")
mm
#Mom_hyper_DMRs_inter_Mom_Up_gene: "PRR22"     "LINC00299" "ETS2"  
genelist<-intersect(listC,listB)
mm<-paste(nameC,"_inter_",nameB,":","  ",genelist,sep="")
mm
#Mom_hypo_DMRs_inter_Mom_Down_gene: "ANTXRL" "UGT2A3" "CAMK4"  "SDK1"

###确定目的差异基因关联DMRs位置
plot_gene<-intersect(listC,listA)
plot_gene2<-intersect(listD,listB)
plot_gene1<-intersect(listD,listA)
plot_gene3<-intersect(listC,listB)
inter_genes<-c(plot_gene,plot_gene1,plot_gene2,plot_gene3)
unique(inter_genes)
# [1] "CD177"     "IL20RB"    "AGPAT4"    "PRR22"     "LINC00299" "ETS2"      "NAV2"      "ZNF229"    "SLC28A3"  
#[10] "ANTXRL"    "UGT2A3"    "CAMK4"     "SDK1"  
Mom_data_region_merge_DMR_all_big[which(Mom_data_region_merge_DMR_all_big$SYMBOL %in% inter_genes),c("SYMBOL","annotation","meth.diff")]
#      SYMBOL                                                  annotation meth.diff
# 63     ANTXRL                                                      3' UTR -21.42057
# 83       NAV2                Exon (ENST00000360655.8/89797, exon 5 of 38)  15.09471
#250     CD177                                            Promoter (<=1kb) -16.41116
#252    ZNF229              Intron (ENST00000588655.1/7772, intron 4 of 5)  15.08807
#256     PRR22                                            Promoter (<=1kb)  21.45220
#290 LINC00299 Intron (ENST00000669954.1/ENST00000669954.1, intron 1 of 7)  16.26291
#313      ETS2                                           Distal Intergenic  19.96044
#356    IL20RB                                           Distal Intergenic -17.44573
#378    UGT2A3 Intron (ENST00000608365.1/ENST00000608365.1, intron 1 of 1) -16.01507
#384     CAMK4              Intron (ENST00000512453.5/814, intron 6 of 11) -15.52398
#422    AGPAT4             Intron (ENST00000320285.9/56895, intron 7 of 8) -16.99776
#449      SDK1                                           Distal Intergenic -37.94176
#509   SLC28A3 Intron (ENST00000650453.1/ENST00000650453.1, intron 1 of 6)  16.68969
#510   SLC28A3 Intron (ENST00000650453.1/ENST00000650453.1, intron 1 of 6)  18.62711
#511   SLC28A3 Intron (ENST00000650453.1/ENST00000650453.1, intron 1 of 6)  17.52775

Mom_exp_file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file9_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
Mom_log2_normalized_Mat<-read.table(Mom_exp_file,sep="\t", header=T)
head(Mom_log2_normalized_Mat);dim(Mom_log2_normalized_Mat)# 18415    20
range(Mom_log2_normalized_Mat[,2:11])#0.00000 18.55641
plot_trans<-Mom_log2_normalized_Mat[which(Mom_log2_normalized_Mat$Row.names %in% inter_genes),]
plot_trans$genetype<-as.character(plot_trans$Row.names)
plot_trans[plot_trans$genetype %in% plot_gene, ]$genetype<-"UP_hypo"
plot_trans[plot_trans$genetype %in% plot_gene1, ]$genetype<-"UP_hyper"
plot_trans[plot_trans$genetype %in% plot_gene2, ]$genetype<-"Down_hyper"
plot_trans[plot_trans$genetype %in% plot_gene3, ]$genetype<-"Down_hypo"
rownames(plot_trans)<-plot_trans$Row.names
plot_trans<-plot_trans[,-1]
plot_trans<-plot_trans[inter_genes,]
#plot_trans<-normalized_counts[plot_gene,]
plot_trans$genetype<-factor(plot_trans$genetype,level = unique(as.character(plot_trans$genetype)))
plot_trans<-plot_trans[order(plot_trans$genetype,decreasing = T),]
Gene_expression<-plot_trans
head(Gene_expression)

#for RNA expression
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 5),rep("Young", 5))))
rownames(annotation_col) = colnames(Gene_expression[1:10])
#for singlegenes
annotation_row = data.frame(Class = factor(plot_trans$genetype))
rownames(annotation_row) = rownames(plot_trans)

## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(UP_hypo =ppCor[6],UP_hyper =ppCor[3], Down_hyper =ppCor[4],Down_hypo=ppCor[5]))
##展示行或者列的label
Gene_expression<-plot_trans[,1:10]

heat_plot1<-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                     #  clustering_distance_rows ="euclidean",#correlation
                     show_rownames = T,show_colnames = T,
                     annotation_col = annotation_col,annotation_row=annotation_row,
                     annotation_colors = ann_colors, 
                     gaps_row =c(4,7,10),gaps_col =c(5),cutree_col = 2,
                     #treeheight_col = 20, #treeheight_row = 30, 
                     # labels_row = labels_row,
                     #border_color ="red", 
                     border=FALSE,
                     color = colorRampPalette(c("lightblue","white","red"))(50),
                     main ="Moms Gene expression of selected meth_trans_inter Genes_merge:original value",angle_col ="45")
#scale by row
heat_plot2 <-pheatmap(Gene_expression, cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_row=annotation_row,
                      annotation_colors = ann_colors, 
                      gaps_row =c(4,7,10),gaps_col =c(5),cutree_col = 2,
                      #treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      main ="Moms Gene expression level of meth_trans_inter Genes: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/Mom_RNA_DEGs_pre_DMRs_heatmaps1.pdf")
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/Mom_RNA_DEGs_pre_DMRs_heatmaps2.pdf")
print(heat_plot2)
dev.off()
#dev.new()

head(Mom_data_region_merge_DMR_all_big);dim(Mom_data_region_merge_DMR_all_big)
Mom_data_region_merge_DMR_all_big[which(Mom_data_region_merge_DMR_all_big$SYMBOL %in% inter_genes),]

plot_DMRs<-Mom_data_region_merge_DMR_all_big[which(Mom_data_region_merge_DMR_all_big$SYMBOL %in% inter_genes),]
plot_DMRs$DMRs<-plot_DMRs$Row.names
plot_DMRs$SYMBOL<-as.character(plot_DMRs$SYMBOL)

plot_DMRs$SYMBOL[duplicated(plot_DMRs$SYMBOL)]
#"SLC28A3" "SLC28A3"
plot_DMRs$SYMBOL2<-plot_DMRs$SYMBOL
plot_DMRs[plot_DMRs$SYMBOL == "SLC28A3",]$SYMBOL2<-c("SLC28A3(1)","SLC28A3(2)","SLC28A3(3)")

rownames(plot_DMRs)<-plot_DMRs$SYMBOL2
plot_DMRs[plot_DMRs$SYMBOL %in% plot_gene, ]$SYMBOL<-"UP_hypo"
plot_DMRs[plot_DMRs$SYMBOL %in% plot_gene1, ]$SYMBOL<-"UP_hyper"
plot_DMRs[plot_DMRs$SYMBOL %in% plot_gene2, ]$SYMBOL<-"Down_hyper"
plot_DMRs[plot_DMRs$SYMBOL %in% plot_gene3, ]$SYMBOL<-"Down_hypo"
#ICR_meth_level<-plot_DMRs[1:7]

sample_list<-c("10-1-in-9","16-1","17-1-2D","18-1","19-1-2G","20M-4F","6-1","7-1","8-1","E16M",
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

plot_DMRs2$SYMBOL<-factor(plot_DMRs2$SYMBOL,levels = rev(c("UP_hypo","UP_hyper","Down_hyper","Down_hypo")))
plot_DMRs2<-plot_DMRs2[order(plot_DMRs2$SYMBOL),]

annotation_row = data.frame( Class = factor(plot_DMRs2$SYMBOL))
rownames(annotation_row) = rownames(plot_DMRs2)

## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(UP_hypo =ppCor[6],UP_hyper =ppCor[3], Down_hyper =ppCor[4],Down_hypo=ppCor[5]))

grid.newpage(); #清空画板，开始画新图
heat_plot3 <-pheatmap(plot_DMRs2[,sample_list2], cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_row=annotation_row,
                      annotation_colors = ann_colors, 
                      gaps_row =c(4,9,12),gaps_col =c(10),cutree_col = 2,
                      treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","pink","orange","firebrick3"))(50),
                      # filename = "D:/3.秦萌高龄甲基化/1.RRBS_result/1.DMR_call/DMR_methlevel_genes_over_DEGs_single_heatmap.pdf",
                      main ="Moms DNA methylation level of meth_trans_inter genes original value",angle_col ="90")

#scale by row
heat_plot4 <-pheatmap(plot_DMRs2[,sample_list2], cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_row=annotation_row,
                      annotation_colors = ann_colors, 
                      gaps_row =c(4,9,12),gaps_col =c(10),cutree_col = 2,
                      treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("purple","white","orange"))(40),
                      #color = colorRampPalette(c("blue","lightgrey","red"))(50),
                      main ="Moms DNA methylation level of selected  DMRs :scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/Mom_meth_DEGs_pre_DMRs_heatmaps1.pdf")
print(heat_plot3)
dev.off()
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/Mom_meth_DEGs_pre_DMRs_heatmaps2.pdf")
print(heat_plot4)
dev.off()