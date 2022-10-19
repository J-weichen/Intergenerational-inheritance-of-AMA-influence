#甲基化对应基因列表的读取
#data preparation
#read kid DMRs
compare_name1<-"kids_AMA_vs_Young"
split_region<- "200bp"
kids_data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t",check.names = F)
kids_DMR_hyper<-kids_data_region_merge_DMR_all_big[kids_data_region_merge_DMR_all_big$meth.diff >0,]
kids_DMR_hypo<-kids_data_region_merge_DMR_all_big[kids_data_region_merge_DMR_all_big$meth.diff <0,]
#read all bins  filter 8 sample per group among all sample
all_sampel_data_region_merge_DMR_all_big<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.csv",row.names = 1,header=T,sep=",",check.names = F)
all_sampel_data_region_merge_DMR_all_big[1:3,]
all_sampel_data_region_merge_DMR_all_big$bin_region<-unlist(lapply(strsplit(rownames(all_sampel_data_region_merge_DMR_all_big),"[.]"), function(x) paste(x[1],x[2],x[3],sep = "_")))

kids_data_region_merge_DMR_all_big[1:3,]
kids_DMR_hyper$Row.names <- as.character(kids_DMR_hyper$Row.names)
kids_DMR_hypo$Row.names <- as.character(kids_DMR_hypo$Row.names)

#construct analysis data frame
data_three_kids_DMR_hyper<-all_sampel_data_region_merge_DMR_all_big[which(all_sampel_data_region_merge_DMR_all_big$bin_region  %in%  kids_DMR_hyper$Row.names),]
data_three_kids_DMR_hypo <-all_sampel_data_region_merge_DMR_all_big[which(all_sampel_data_region_merge_DMR_all_big$bin_region  %in%  kids_DMR_hypo$Row.names),]

data_three_kids_DMR <- rbind(data_three_kids_DMR_hyper,data_three_kids_DMR_hypo)
rownames(data_three_kids_DMR)<-data_three_kids_DMR$bin_region
data_three_kids_DMR2<-data_three_kids_DMR[,c(1:60)]

data_three_kids_DMR3 <- as.data.frame(t(data_three_kids_DMR2))
head(data_three_kids_DMR3)
data_three_kids_DMR3$sample<-rownames(data_three_kids_DMR3)
methy_plot_long<- melt(data_three_kids_DMR3,id.vars=c("sample"),variable.name="Sample",value.name = "meth_level")
head(methy_plot_long)
methy_plot_long$library_code <-methy_plot_long$sample
#reading meta data 
huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
rownames(huahua_meta)<-huahua_meta$library_code
head(huahua_meta)

data_three_kids_DMR4<-merge(huahua_meta,methy_plot_long,by="library_code")
head(data_three_kids_DMR4)
data_three_kids_DMR5<-data_three_kids_DMR4[,c("Family","Sample","Sample.types","meth_level")]
data_three_kids_DMR6 <- dcast(data_three_kids_DMR5,Family+Sample ~ Sample.types , value.var = "meth_level")
head(data_three_kids_DMR6)

#calculation the correlationship for kid and father 
data_analysis<-data_three_kids_DMR6[,c("Sample","Family","kids","father")]
colnames(data_analysis)<-c("DMR_region","Family","kids","parent")

DMR_name<-c();cor_r<-c();cor_pvalue <- c()
for (i in unique(data_analysis$DMR_region)){
  #i="chr1_1954801_1955000"
  print(i)
  data_cal<-na.omit(data_analysis[which(data_analysis$DMR_region == i),c("DMR_region","kids","parent")])
  #dim(data_cal)
  c_r   <- cor(as.numeric(data_cal$kids),as.numeric(data_cal$parent),method="spearman")
  p_vlue<- cor.test(as.numeric(data_cal$kids),as.numeric(data_cal$parent),method="spearman")[[3]]
  DMR_name<-c(DMR_name,i)
  cor_r<-c(cor_r,c_r)
  cor_pvalue<-c(cor_pvalue,p_vlue)
}
length(cor_r);length(DMR_name);length(cor_pvalue)
cor_data_df<-data.frame(DMRs_code=DMR_name,correlation=cor_r,pvalue=cor_pvalue)
names(cor_data_df) <- c("DMR_region","correlation","pvalue")
head(cor_data_df)
hist(cor_data_df$correlation)
cor_data_df$group<-"no_sig"
cor_data_df[which(abs(cor_data_df$correlation) >= 0.600 & cor_data_df$pvalue < 0.05),]$group <- "candidate_region"
table(cor_data_df$group)
#candidate_region           no_sig 
#           58              393 
data_analysis2<-merge(data_analysis,cor_data_df,by="DMR_region")
head(data_analysis2)
write.table(data_analysis2, file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_father.txt",row.names=T, col.names=T,sep="\t") 
write.table(cor_data_df, file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlated_pvalue_for_kid_AMA_DMR_between_kids_and_father.txt",row.names=T, col.names=T,sep="\t") 

#calculation the correlationship for kid and mother 
data_analysis<-data_three_kids_DMR6[,c("Sample","Family","kids","mother")]
colnames(data_analysis)<-c("DMR_region","Family","kids","parent")

DMR_name<-c();cor_r<-c();cor_pvalue <- c()
for (i in unique(data_analysis$DMR_region)){
  #i="chr1_1954801_1955000"
  print(i)
  data_cal<-na.omit(data_analysis[which(data_analysis$DMR_region == i),c("DMR_region","kids","parent")])
  #dim(data_cal)
  c_r   <- cor(as.numeric(data_cal$kids),as.numeric(data_cal$parent),method="spearman")
  p_vlue<- cor.test(as.numeric(data_cal$kids),as.numeric(data_cal$parent),method="spearman")[[3]]
  DMR_name<-c(DMR_name,i)
  cor_r<-c(cor_r,c_r)
  cor_pvalue<-c(cor_pvalue,p_vlue)
}
length(cor_r);length(DMR_name);length(cor_pvalue)#451
cor_data_df<-data.frame(DMRs_code=DMR_name,correlation=cor_r,pvalue=cor_pvalue)
names(cor_data_df) <- c("DMR_region","correlation","pvalue")
head(cor_data_df)
hist(cor_data_df$correlation)
cor_data_df$group<-"no_sig"
cor_data_df[which(abs(cor_data_df$correlation) >= 0.600 & cor_data_df$pvalue < 0.05),]$group <- "candidate_region"
table(cor_data_df$group)
#candidate_region           no_sig 
#            66              385 
data_analysis2<-merge(data_analysis,cor_data_df,by="DMR_region")
head(data_analysis2)
cor_data_df[which(cor_data_df$group =="candidate_region"),]
write.table(data_analysis2, file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_mother.txt",row.names=T, col.names=T,sep="\t") 
write.table(cor_data_df, file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlated_pvalue_for_kid_AMA_DMR_between_kids_and_mother.txt",row.names=T, col.names=T,sep="\t") 

#plot for correlated data 
#for father                                                                         
data_analysis <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_father.txt",header=T,sep="\t")
head(data_analysis);
table(data_analysis$group)
#cor_data_df2<-cor_data_df %>% filter(pvalue < 0.05)
#data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region" & data_analysis$pvalue<0.01),]
data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region"),]
#formal plotting 
library(ggpubr)
#meth_tran_cor4<-data_analysis[which(data_analysis$DMRs_code %in% cor_data_df2$DMRs_code),]
#meth_tran_cor4$gene_DMRs<-paste(meth_tran_cor4$gene,meth_tran_cor4$DMRs_region,sep="-")
plot_cordata<-data_analysis_sig[,c("DMR_region","kids","parent","Family")]
plot_cordata$DMR_region <-as.character(plot_cordata$DMR_region)

ceiling(length(unique(plot_cordata$DMR_region))/8)
cor_plot0<-ggplot(data=plot_cordata, aes(x=kids, y=parent))+geom_point()+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=plot_cordata, method = "spearman")+  facet_wrap(~ DMR_region, scales = "free",ncol = 8)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_siginificant_AMA_DMR_between_kids_and_father_seprated0.pdf",cor_plot0,width = 60*8, height =60*ceiling(length(unique(plot_cordata$DMR_region))/8), units = "mm")

#for selected DMRs
pos_DMRs_selected <- unique(as.character(data_analysis_sig[which(data_analysis_sig$correlation>0),]$DMR_region))
neg_DMRs_selected <- unique(as.character(data_analysis_sig[which(data_analysis_sig$correlation<0),]$DMR_region))
plot_cordata$type<-"postive_relation"
plot_cordata[which(plot_cordata$DMR_region %in% neg_DMRs_selected),]$type<-"negtive_relation"

cor_plot2<-ggscatter(plot_cordata[which(plot_cordata$DMR_region %in% c(pos_DMRs_selected,neg_DMRs_selected)),], 
                     x = "kids", y = "parent",
                     add =  "reg.line",                         # Add regression line
                     conf.int = TRUE,                          # Add confidence interval
                     color = "DMR_region" # palette =ppCor_all2,           # Color by groups "cyl"
                     #   facet.by = "gene_DMRs",
                     #)+stat_cor(aes(color = gene_DMRs), label.x = 3)           # Add correlation coefficient
)+stat_cor(aes(color = DMR_region), method = "spearman", label.x = 3)           # Add correlation coefficient
cor_plot2
#cor_plot3<-cor_plot2+ facet_wrap(~ type, scales = "free",ncol=2)
cor_plot4<-cor_plot2+ facet_wrap(~ DMR_region, scales = "free",ncol=8)+
  theme(legend.title = element_text(color = "black", size = 7),legend.position = "right")
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_siginificant_AMA_DMR_between_kids_and_father_seprated.pdf",cor_plot2,width = 400, height =300, units = "mm")
#ggsave(file="D:/1.huahua/3.correlated_DMRs/correlationship_methlation_level_of_siginificant_AMA_DMR_between_kids_and_father_seprated2.pdf",cor_plot3,width = 800, height =300, units = "mm")
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_siginificant_AMA_DMR_between_kids_and_father_seprated3.pdf",cor_plot4,width = 60*8, height =60*ceiling(length(unique(plot_cordata$DMR_region))/8), units = "mm")

#for mother 
data_analysis <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_mother.txt",header=T,sep="\t")
head(data_analysis);
table(data_analysis$group)
#cor_data_df2<-cor_data_df %>% filter(pvalue < 0.05)
data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region"),]
#data_analysis_sig<-data_analysis[which(data_analysis$correlation >= 0.7 & data_analysis$pvalue<0.01),]

#formal plotting 
library(ggpubr)
#meth_tran_cor4<-data_analysis[which(data_analysis$DMRs_code %in% cor_data_df2$DMRs_code),]
#meth_tran_cor4$gene_DMRs<-paste(meth_tran_cor4$gene,meth_tran_cor4$DMRs_region,sep="-")
plot_cordata<-data_analysis_sig[,c("DMR_region","kids","parent","Family")]
plot_cordata$DMR_region <-as.character(plot_cordata$DMR_region)

ceiling(length(unique(plot_cordata$DMR_region))/8)
cor_plot0<-ggplot(data=plot_cordata, aes(x=kids, y=parent))+geom_point()+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=plot_cordata, method = "spearman")+  facet_wrap(~ DMR_region, scales = "free",ncol = 9)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_siginificant_AMA_DMR_between_kids_and_mother_seprated0.pdf",cor_plot0,width = 60*9, height =60*ceiling(length(unique(plot_cordata$DMR_region))/9), units = "mm")

#for selected DMRs
pos_DMRs_selected <- unique(as.character(data_analysis_sig[which(data_analysis_sig$correlation>0),]$DMR_region))
neg_DMRs_selected <- unique(as.character(data_analysis_sig[which(data_analysis_sig$correlation<0),]$DMR_region))
plot_cordata$type<-"postive_relation"
plot_cordata[which(plot_cordata$DMR_region %in% neg_DMRs_selected),]$type<-"negtive_relation"

cor_plot2<-ggscatter(plot_cordata[which(plot_cordata$DMR_region %in% c(pos_DMRs_selected,neg_DMRs_selected)),], 
                     x = "kids", y = "parent",
                     add = "reg.line",                         # Add regression line
                     conf.int = TRUE,                          # Add confidence interval
                     color = "DMR_region" # palette =ppCor_all2,           # Color by groups "cyl"
                     #   facet.by = "gene_DMRs",
                     #)+stat_cor(aes(color = gene_DMRs), label.x = 3)           # Add correlation coefficient
)+stat_cor(aes(color = DMR_region), method = "spearman", label.x = 3)           # Add correlation coefficient
cor_plot2
cor_plot3<-cor_plot2+ facet_wrap(~ type, scales = "free",ncol=2)
cor_plot4<-cor_plot2+ facet_wrap(~ DMR_region, scales = "free",ncol=8)+
  theme(legend.title = element_text(color = "black", size = 7),legend.position = "right")
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_siginificant_AMA_DMR_between_kids_and_mother_seprated.pdf",cor_plot2,width = 400, height =300, units = "mm")
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_siginificant_AMA_DMR_between_kids_and_mother_seprated2.pdf",cor_plot3,width = 800, height =300, units = "mm")
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_siginificant_AMA_DMR_between_kids_and_mother_seprated3.pdf",cor_plot4,width = 60*8, height =60*ceiling(length(unique(plot_cordata$DMR_region))/8), units = "mm")

#show genes for cor_DMRs
data_analysis_MK <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_mother.txt",header=T,sep="\t")
data_analysis_FK <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_father.txt",header=T,sep="\t")
data_analysis_all<-rbind(data_analysis_MK,data_analysis_FK)
data_analysis_sig<-data_analysis_all[which(data_analysis_all$group == "candidate_region"),]
data_region_merge_DMR_all_big<-read.table("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/kids_AMA_vs_Young_200bp_DMR_myDiff15q005_merge_data.txt",header=T,sep="\t")
DMRs_unique<-unique(as.character(data_analysis_sig$DMR_region))
cor_DMR_gene<-data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$Row.names %in% DMRs_unique),c("bin_region","distanceToTSS","annotation","SYMBOL")]
dim(cor_DMR_gene)#118   4
write.table(cor_DMR_gene, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/all_cor_DMR_gene_kids_AMA_vs_Young_200bp_DMR_myDiff15q005.txt"),quote=F,row.names=F,col.names=T,sep = "\t")
