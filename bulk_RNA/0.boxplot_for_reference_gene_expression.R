rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
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
##not log2
file <- c("/mnt/data/chenwei/huahua/5.valificat_result/file5_all_filter_Mom_verification_AMA_deseq2_normalized_count.txt")
Mom_merge_data <- read.table(file=file,header=T)
head(Mom_merge_data)
dim(Mom_merge_data)
Mom_draw_count<-Mom_merge_data

#FOR kids
#read expression data 
tag_name<-"Kid"; sampleA <-"AMA";sampleB <-"Young" 
file <- c("/mnt/data/chenwei/huahua/5.valificat_result/file5_all_filter_Kid_verification_AMA_deseq2_normalized_count.txt")
kids_merge_data <- read.table(file=file,header=T)
head(kids_merge_data)
kids_draw_count<-kids_merge_data
##reference genes
Mom_kid_common_DEGs<-c("ACTB","GAPDH","PPIB")

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

#colData_whole<-colData_whole[which(!(colData_whole$sample %in% c("Y5Q","Y6Q","A2Q","Y5M","Y6M","A2M"))),]
head(colData_whole)
colData_whole$library_code<-colData_whole$sample

MK_common_DEG_data3<-merge(colData_whole,MK_common_DEG_data_long,by="library_code")
MK_common_DEG_data3$Family<-substr(MK_common_DEG_data3$library_code,1,2)
head(MK_common_DEG_data3)

data_plot<-MK_common_DEG_data3
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
  ylim(0,max(data_plot$trans_level)+20000)+facet_wrap(~ gene_name, scales = "free",ncol=3)
P_mean_site2<-P_mean_site+stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format",label.y=max(data_plot$trans_level)+5000) +labs(title="mRNA expression level for AMA related genes(wilcox.test)", x ="sample_group", y ="log2(count+1)")
P_mean_site3<-P_mean_site+stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.format",label.y=max(data_plot$trans_level)+5000) + labs(title="mRNA expression level for AMA related genes(t.test)", x ="sample_group", y ="log2(count+1)")
P_mean_site2
P_mean_site3
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/nolog2_count_simple_all_trans_level_of_reference_gene_between_kids_and_mother_seprated.pdf",P_mean_site,width = 8, height =5)
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/nolog2_count_wilcox_test_for_all_trans_level_of_reference_gene_between_kids_and_mother_seprated.pdf",P_mean_site2,width = 8, height =5)
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/nolog2_count_T_test_for_all_trans_level_of_reference_gene_between_kids_and_mother_seprated.pdf",P_mean_site3,width = 8, height =5)

##answer for reviewer
P_mean_site<-ggboxplot(data_plot[which(data_plot$gene_name =="ACTB"),], x = "group_type", y = "trans_level", color = "group_type", palette = ppCor, add = "jitter")+#ylim(c(70,85))+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=45,hjust=0.5, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(plot.title = element_text(size=10,colour = "black",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),
        axis.title.y = element_text(size=10,colour = "black",face = "bold"),
        axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))+
         ylim(5000,350000)
P_mean_site3<-P_mean_site+stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif",label.y=300000) + labs(title="mRNA expression level for AMA related genes(t.test)", x ="sample_group", y =" normalized count")
P_mean_site3
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/Response_nolog2_count_T_test_for_all_trans_level_of_ACTB_between_kids_and_mother_seprated.pdf",P_mean_site3,width = 5, height =5)

###for log2

#FOR Mom
#read expression data 
tag_name<-"Mom"; sampleA <-"AMA";sampleB <-"Young" 
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file9_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".log2_normalized_count_add_gene_information.txt")
##not log2
file <- c("/mnt/data/chenwei/huahua/5.valificat_result/file5_all_filter_Mom_verification_AMA_deseq2_normalized_count.txt")
Mom_merge_data <- read.table(file=file,header=T)
head(Mom_merge_data)
rownames(Mom_merge_data)<-Mom_merge_data$Row.names
Mom_merge_data<-Mom_merge_data[,-1]
Mom_merge_data[1:6,1:6]
Mom_draw_count<-Mom_merge_data[,1:10]


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
##reference genes
Mom_kid_common_DEGs<-c("ACTB","GAPDH","PPIB")

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

#colData_whole<-colData_whole[which(!(colData_whole$sample %in% c("Y5Q","Y6Q","A2Q","Y5M","Y6M","A2M"))),]
head(colData_whole)
colData_whole$library_code<-colData_whole$sample

MK_common_DEG_data3<-merge(colData_whole,MK_common_DEG_data_long,by="library_code")
MK_common_DEG_data3$Family<-substr(MK_common_DEG_data3$library_code,1,2)
head(MK_common_DEG_data3)

data_plot<-MK_common_DEG_data3
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
P_mean_site2
P_mean_site3
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/simple_all_trans_level_of_reference_gene_between_kids_and_mother_seprated.pdf",P_mean_site,width = 8, height =5)
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/wilcox_test_for_all_trans_level_of_reference_gene_between_kids_and_mother_seprated.pdf",P_mean_site2,width = 8, height =5)
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/T_test_for_all_trans_level_of_reference_gene_between_kids_and_mother_seprated.pdf",P_mean_site3,width = 8, height =5)

