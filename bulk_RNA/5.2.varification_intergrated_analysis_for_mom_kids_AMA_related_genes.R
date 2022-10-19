rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(VennDiagram)
library(scales)
library(ggsci)
library(pheatmap)
library(ggpubr)

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

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

#read DEGs  
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file13_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
Mom_AMA_DEGs<-read.table(file)
head(Mom_AMA_DEGs)
Mom_AMA_DEGs_all<-Mom_AMA_DEGs$ID
Mom_AMA_DEGs_up_1.5<-Mom_AMA_DEGs[which(Mom_AMA_DEGs$log2FoldChange>0),]$ID
Mom_AMA_DEGs_down_1.5<-Mom_AMA_DEGs[which(Mom_AMA_DEGs$log2FoldChange<0),]$ID

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

##disscussion relationship
Mom_kid_overlap_DEGs<-Reduce(intersect,list(Mom_AMA_DEGs_all,kids_AMA_DEGs_all))
common_UP_DEGs<-Reduce(intersect,list(kids_AMA_DEGs_up_1.5,Mom_AMA_DEGs_up_1.5))
common_Down_DEGs<-Reduce(intersect,list(kids_AMA_DEGs_down_1.5,Mom_AMA_DEGs_down_1.5))
length(Mom_kid_overlap_DEGs)#151
length(common_UP_DEGs);length(common_Down_DEGs)#39 #70
Mom_kid_common_DEGs<-c(common_UP_DEGs,common_Down_DEGs)

Mom_kid_common_DEGsdata<-data.frame(MK_DEGs=Mom_kid_common_DEGs,trend=c(rep("common_UP_DEGs",length(common_UP_DEGs)),rep("common_Down_DEGs",length(common_Down_DEGs))))
write.table(Mom_kid_common_DEGsdata, file="/mnt/data/chenwei/huahua/5.valificat_result/AMA_common_DEG_between_kids_and_mother_trend_add.txt",row.names=T, col.names=T,sep="\t") 

venn_Kids_Mom_DEGs1 <-venn.diagram(list(Mom_AMA_DEGs_all=Mom_AMA_DEGs_all,kids_AMA_DEGs_all=kids_AMA_DEGs_all),
                                  alpha=c(0.9,0.9),lwd=1,lty=1,col="black" , 
                                  fill=ppCor[3:4],cex = 1.5,cat.col=ppCor[3:4],cat.fontface=4, cat.cex = 1.5,    
                                  main = "Kids and Mom AMA related DEGs",
                                  main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                  filename = NULL)

venn_Kids_Mom_DEGs2<-venn.diagram(list(Mom_DEGs_up=Mom_AMA_DEGs_up_1.5,Mom_DEGs_down=Mom_AMA_DEGs_down_1.5,
                                        kids_AMA_DEGs_up=kids_AMA_DEGs_up_1.5,kids_AMA_DEGs_down=kids_AMA_DEGs_down_1.5),
                                   alpha=c(1,1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1:3,5)],cex = 1.5, cat.col=ppCor[c(1:3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, 
                                   main = "Kids and Mom AMA related DEGs",main.fontfamily = 3,filename = NULL)
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Venn_Mom_RNA_Mom_kid_commom_DEGs.pdf",width = 8,height = 6)
grid.newpage();
grid.draw(venn_Kids_Mom_DEGs1)
grid.newpage();
grid.draw(venn_Kids_Mom_DEGs2)
dev.off()

#for mom RNA expression
head(Mom_draw_count)
annotation_col<-data.frame(Treatment = factor(c(rep("AMA", 5),rep("Young", 5))))
rownames(annotation_col) = colnames(Mom_draw_count)
## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))
##展示行或者列的label

Gene_expression<-Mom_draw_count[which(rownames(Mom_draw_count) %in% Mom_kid_common_DEGs),]
Gene_expression<-Gene_expression[Mom_kid_common_DEGs,]

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
                     main ="Mom_kid_overlap_DEGs mom_expression:original value",angle_col ="45")
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
                      main ="Mom_kid_overlap_DEGs mom_expression: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Mom_RNA_Mom_kid_commom_DEGs_heatmap1.pdf",width = 12,height = 15)
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Mom_RNA_Mom_kid_commom_DEGs_heatmap2.pdf",width = 12,height = 15)
print(heat_plot2)
dev.off()

#for Kids RNA expression
head(kids_draw_count)
annotation_col<-data.frame(Treatment = factor(c(rep("AMA",5),rep("Young", 5))))
rownames(annotation_col) = colnames(kids_draw_count)
## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]))
##展示行或者列的label

Gene_expression<-kids_draw_count[which(rownames(kids_draw_count) %in% Mom_kid_common_DEGs),]
Gene_expression<-Gene_expression[Mom_kid_common_DEGs,]
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
                     main ="Mom_kid_overlap_DEGs kids_expression:original value",angle_col ="45")
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
                      main ="Mom_kid_overlap_DEGs kids_expression: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Kid_RNA_Mom_kid_common_DEGs_heatmap1.pdf",width = 12,height = 15)
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Kid_RNA_Mom_kid_common_DEGs_heatmap2.pdf",width = 12,height = 15)
print(heat_plot2)
dev.off()

length(Mom_kid_common_DEGs)
kids_Gene_expression<-kids_draw_count[which(rownames(kids_draw_count) %in% Mom_kid_common_DEGs),]
Mom_Gene_expression<-Mom_draw_count[which(rownames(Mom_draw_count) %in% Mom_kid_common_DEGs),]

MK_common_DEG_data<-merge(kids_Gene_expression,Mom_Gene_expression,by=0)
rownames(MK_common_DEG_data)<-MK_common_DEG_data$Row.names;MK_common_DEG_data<-MK_common_DEG_data[,-1]
head(kids_Gene_expression);head(Mom_Gene_expression);head(MK_common_DEG_data)
MK_common_DEG_data2 <- as.data.frame(t(MK_common_DEG_data))
head(MK_common_DEG_data2)
MK_common_DEG_data2$sample<-rownames(MK_common_DEG_data2)
MK_common_DEG_data_long<- melt(MK_common_DEG_data2,id.vars=c("sample"),variable.name="gene_name",value.name = "trans_level")
MK_common_DEG_data_long$library_code <-MK_common_DEG_data_long$sample
head(MK_common_DEG_data_long)

colData_whole<-read.csv(file="/mnt/data/chenwei/huahua/file3_all_no_filter_verification_AMA_analysis_metadata.csv",row.names=1,header =T)
colData_whole<-colData_whole[which(!(colData_whole$sample %in% c("Y5Q","Y6Q","A2Q","Y5M","Y6M","A2M"))),]
head(colData_whole)
colData_whole$library_code<-colData_whole$sample

MK_common_DEG_data3<-merge(colData_whole,MK_common_DEG_data_long,by="library_code")
MK_common_DEG_data3$Family<-substr(MK_common_DEG_data3$library_code,1,2)
head(MK_common_DEG_data3)

MK_common_DEG_data4<-MK_common_DEG_data3[,c("Family","gene_name","generation","trans_level")]
MK_common_DEG_data5 <- dcast(MK_common_DEG_data4,Family+gene_name ~ generation , value.var = "trans_level")
head(MK_common_DEG_data5)

#calculation the correlationship for kid and mother 
data_analysis<-MK_common_DEG_data5[,c("gene_name","Family","Kid","Mom")]
colnames(data_analysis)<-c("Gene","Family","kids","parent")
##配对样本相关性分析
Gene_name<-c();cor_r<-c();cor_pvalue <- c()
for (i in unique(data_analysis$Gene)){
  #i="ACTB"
  print(i)
  data_cal<-na.omit(data_analysis[which(data_analysis$Gene == i),c("Gene","kids","parent")])
  #dim(data_cal)
  c_r   <- cor(as.numeric(data_cal$kids),as.numeric(data_cal$parent),method="spearman")
  p_vlue<- cor.test(as.numeric(data_cal$kids),as.numeric(data_cal$parent),method="spearman")[[3]]
  Gene_name<-c(Gene_name,i)
  cor_r<-c(cor_r,c_r)
  cor_pvalue<-c(cor_pvalue,p_vlue)
}
length(cor_r);length(Gene_name);length(cor_pvalue)
#padj=p.adjust(cor_pvalue, method = "bonferroni", n = length(cor_pvalue))
cor_data_df<-data.frame(Gene_name=Gene_name,correlation=cor_r,pvalue=cor_pvalue)

names(cor_data_df) <- c("Gene","correlation","pvalue")
head(cor_data_df)
hist(cor_data_df$correlation)
cor_data_df$group<-"no_sig"
cor_data_df[which(abs(cor_data_df$correlation) >= 0.600 & cor_data_df$pvalue < 0.05),]$group <- "candidate_region"
table(cor_data_df$group)
#candidate_region           no_sig 
#          48               61 
data_analysis2<-merge(data_analysis,cor_data_df,by="Gene")
head(data_analysis2)
cor_data_df[which(cor_data_df$group =="candidate_region"),]
write.table(data_analysis2, file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/correlationship_trans_level_of_AMA_common_DEG_between_kids_and_mother.txt",row.names=T, col.names=T,sep="\t") 
write.table(cor_data_df, file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/correlated_pvalue_of_AMA_common_DEG_between_kids_and_mother.txt",row.names=T, col.names=T,sep="\t") 

data_analysis <-read.table(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/correlationship_trans_level_of_AMA_common_DEG_between_kids_and_mother.txt",header=T,sep="\t")
head(data_analysis)
table(data_analysis$group)
#cor_data_df2<-cor_data_df %>% filter(pvalue < 0.05)
#data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region" & data_analysis$pvalue<0.01),]
data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region"),]
#formal plotting 
plot_cordata<-data_analysis_sig[,c("Gene","kids","parent","Family")]
plot_cordata$Gene <-as.character(plot_cordata$Gene)
#plot_cordata<-data_analysis_sig[which(data_analysis_sig$Gene %in% c("SCD","CD24")),]
ceiling(length(unique(plot_cordata$Gene))/8)
cor_plot0<-ggplot(data=plot_cordata, aes(x=kids, y=parent))+geom_point()+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=plot_cordata, method = "spearman")+  facet_wrap(~ Gene, scales = "free",ncol = 8)
#ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/correlationship_trans_level_of_siginificant_AMA_common_DEGs_between_kids_and_mother_seprated0.pdf",cor_plot0,width = 3*8, height =3*ceiling(length(unique(plot_cordata$DMR_region))/8))
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/correlationship_trans_level_of_siginificant_AMA_common_DEGs_between_kids_and_mother_seprated0.pdf",cor_plot0,width = 24, height =20)

#for selected genes
pos_Gene_selected <- unique(as.character(data_analysis_sig[which(data_analysis_sig$correlation>0),]$Gene))
neg_Gene_selected <- unique(as.character(data_analysis_sig[which(data_analysis_sig$correlation<0),]$Gene))
plot_cordata$type<-"postive_relation"
plot_cordata[which(plot_cordata$Gene %in% neg_Gene_selected),]$type<-"negtive_relation"

cor_plot2<-ggscatter(plot_cordata[which(plot_cordata$Gene %in% c(pos_Gene_selected,neg_Gene_selected)),], 
                     x = "kids", y = "parent",
                     add = "reg.line",                         # Add regression line
                     conf.int = TRUE,                          # Add confidence interval
                     color = "Gene")+stat_cor(aes(color = Gene), method = "spearman", label.x = 3)           # Add correlation coefficient
cor_plot2
cor_plot4<-cor_plot2+ facet_wrap(~ Gene, scales = "free",ncol=8)+
  theme(legend.title = element_text(color = "black", size = 7),legend.position = "right")
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/correlationship_trans_level_of_siginificant_AMA_DMR_between_kids_and_mother_seprated.pdf",cor_plot2,width = 400, height =300, units = "mm")
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/correlationship_trans_level_of_siginificant_AMA_DMR_between_kids_and_mother_seprated3.pdf",cor_plot4,width = 40, height =40)

###plot for seperated genes
data_plot<-MK_common_DEG_data3[which(MK_common_DEG_data3$gene_name %in% unique(plot_cordata$Gene)),]
unique(data_plot$group_type)
data_plot$group_type<-factor(data_plot$group_type, levels=c( "Young_Kid","AMA_Kid","Young_Mom","AMA_Mom"),ordered=TRUE)
data_plot<-data_plot[order(data_plot$group_type),]
unique(data_plot$group_type)

###个体分开绘制
my_comparisons <- list(c("Young_Kid","AMA_Kid"),c("Young_Mom","AMA_Mom"))
head(data_plot)
P_mean_site<-ggboxplot(data_plot, x = "group_type", y = "trans_level", color = "group_type", palette = ppCor, add = "jitter")+#ylim(c(70,85))+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),
        axis.title.y = element_text(size=10,colour = "black",face = "bold"),
        axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))+
  ylim(0,max(data_plot$trans_level)+1)+facet_wrap(~ gene_name, scales = "free",ncol=8)
P_mean_site2<-P_mean_site+stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format",label.y=max(data_plot$trans_level)+0.5) +labs(title="mRNA expression level for AMA related genes(wilcox.test)", x ="sample_group", y ="log2(count+1)")
P_mean_site3<-P_mean_site+stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.format",label.y=max(data_plot$trans_level)+0.5) + labs(title="mRNA expression level for AMA related genes(t.test)", x ="sample_group", y ="log2(count+1)")
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/simple_all_trans_level_of_siginificant_AMA_common_DEGs_between_kids_and_mother_seprated.pdf",P_mean_site,width = 32, height =25)
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/wilcox_test_for_all_trans_level_of_siginificant_AMA_common_DEGs_between_kids_and_mother_seprated.pdf",P_mean_site2,width = 32, height =25)
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/T_test_for_all_trans_level_of_siginificant_AMA_common_DEGs_between_kids_and_mother_seprated.pdf",P_mean_site3,width = 32, height =25)

#heatmap plot for kids or Mom
MK_common_DEG_heat<-MK_common_DEG_data[which(rownames(MK_common_DEG_data) %in% unique(plot_cordata$Gene)),]
head(MK_common_DEG_heat)
MK_common_DEG_heat$genetype<-"Common_Up"
MK_common_DEG_heat[which(rownames(MK_common_DEG_heat) %in% common_Down_DEGs), ]$genetype<-"Common_Down"
MK_common_DEG_heat$genetype<-factor(MK_common_DEG_heat$genetype,level = c("Common_Down","Common_Up"),ordered = T)
MK_common_DEG_heat<-MK_common_DEG_heat[order(MK_common_DEG_heat$genetype,decreasing = T),]

#for Mom
Gene_expression<-MK_common_DEG_heat[,c(colnames(Mom_draw_count),"genetype")]
annotation_col<-data.frame(Treatment = factor(c(rep("AMA",5),rep("Young", 5))))
rownames(annotation_col) = colnames(Mom_draw_count)
annotation_row = data.frame(Class = factor(Gene_expression$genetype))
rownames(annotation_row) = rownames(Gene_expression)
length(which(Gene_expression$genetype=="Common_Down"))
## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(Common_Up =ppCor[3],Common_Down =ppCor[6]))
##展示行或者列的label
heat_plot1 <-pheatmap(Gene_expression[,colnames(Mom_draw_count)],cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_row=annotation_row,
                      annotation_colors = ann_colors, 
                      gaps_row =c(length(which(Gene_expression$genetype=="Common_Up"))),gaps_col =c(5),cutree_col = 2,
                      #treeheight_col = 20, #treeheight_row = 30, 
                      # labels_row = labels_row,
                      #border_color ="red", 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      main ="Moms Gene expression level of significant correlated DEGs: scale by row",angle_col ="90")
#for Kids
Gene_expression<-MK_common_DEG_heat[,c(colnames(kids_draw_count),"genetype")]
annotation_col<-data.frame(Treatment = factor(c(rep("AMA",5),rep("Young", 5))))
rownames(annotation_col) = colnames(kids_draw_count)
annotation_row = data.frame(Class = factor(Gene_expression$genetype))
rownames(annotation_row) = rownames(Gene_expression)
length(which(Gene_expression$genetype=="Common_Down"))
## 自定义分组的颜色
ann_colors = list(Treatment=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(Common_Up =ppCor[3],Common_Down =ppCor[6]))
##展示行或者列的label
heat_plot2 <-pheatmap(Gene_expression[,colnames(kids_draw_count)],cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
                      clustering_distance_rows ="euclidean",#correlation
                      show_rownames = T,show_colnames = T,
                      annotation_col = annotation_col,annotation_row=annotation_row,
                      annotation_colors = ann_colors, 
                      gaps_row =c(length(which(Gene_expression$genetype=="Common_Up"))),gaps_col =c(5),cutree_col = 2, 
                      scale ="row", 
                      border=FALSE,
                      color = colorRampPalette(c("navy","white","firebrick3"))(50),
                      main ="Kids Gene expression level of significant correlated DEGs: scale by row",angle_col ="90")
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Mom_RNA_Mom_kid_corplot_common_sig_DEGs_heatmap1.pdf",width = 10,height = 10)
print(heat_plot1)
dev.off()
pdf("/mnt/data/chenwei/huahua/5.valificat_result/Kid_RNA_Mom_kid_corplot_common_sig_DEGs_heatmap2.pdf",width = 10,height = 10)
print(heat_plot2)
dev.off()