rm(list = ls())
#载入各种包
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
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)#org.Hs.eg.db 的数据类型以及使用简介：可用于数据类型转换
keytypes(org.Hs.eg.db)#查看org.Hs.eg.db数据对象里面包含着各大主流数据库的数据
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#绘制目标基因的目标区域甲基化水平柱状图
#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)


#reading kids DMRs
compare_name1<-"kids_AMA_vs_Young"
split_region<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("D:/huahua/1.DMR_called/",compare_name1,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t",check.names = F)
kids_alldata<-data_region_merge_DMR_all_big
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
DMR_hyper_kid<-unique(as.character(DMR_hyper$Row.names))
DMR_hypo_kid<-unique(as.character(DMR_hypo$Row.names))
length(DMR_hyper_kid);length(DMR_hypo_kid)#186 #237

#read correlated DMRs
mother_cordata <-read.table(file="D:/huahua/3.correlated_DMRs/correlated_pvalue_for_kid_AMA_DMR_between_kids_and_mother.txt",header=T,sep="\t")
MK_cor_DMR_sig<-as.character(mother_cordata[which(mother_cordata$group == "candidate_region"),]$DMR_region)
father_cordata <-read.table(file="D:/huahua/3.correlated_DMRs/correlated_pvalue_for_kid_AMA_DMR_between_kids_and_father.txt",header=T,sep="\t")
FK_cor_DMR_sig<-as.character(father_cordata[which(father_cordata$group == "candidate_region"),]$DMR_region)

#seperated hyper and hypo cor DMRs
##filter pure related DMRs but not cor DMRs
MK_cor_DMR_sig_hyper_raw<-intersect(DMR_hyper_kid,MK_cor_DMR_sig)
MK_cor_DMR_sig_hypo_raw<-intersect(DMR_hypo_kid,MK_cor_DMR_sig)
FK_cor_DMR_sig_hyper_raw<-intersect(DMR_hyper_kid,FK_cor_DMR_sig)
FK_cor_DMR_sig_hypo_raw<-intersect(DMR_hypo_kid,FK_cor_DMR_sig)

MF_cor_DMRs <- intersect(MK_cor_DMR_sig,FK_cor_DMR_sig)
##find out both parental cor DMRs seprated by hyper and hypo
MF_cor_DMRs_hyper<-Reduce(intersect,c(list(DMR_hyper_kid),list(MF_cor_DMRs)))
MF_cor_DMRs_hypo<-Reduce(intersect,c(list(DMR_hypo_kid),list(MF_cor_DMRs)))
length(MF_cor_DMRs_hyper);length(MF_cor_DMRs_hypo)

MK_cor_DMR_sig_hyper<-setdiff(MK_cor_DMR_sig_hyper_raw,MF_cor_DMRs_hyper)
MK_cor_DMR_sig_hypo<-setdiff(MK_cor_DMR_sig_hypo_raw,MF_cor_DMRs_hypo)
FK_cor_DMR_sig_hyper<-setdiff(FK_cor_DMR_sig_hyper_raw,MF_cor_DMRs_hyper)
FK_cor_DMR_sig_hypo<-setdiff(FK_cor_DMR_sig_hypo_raw,MF_cor_DMRs_hypo)
length(MK_cor_DMR_sig_hyper);length(MK_cor_DMR_sig_hypo);length(FK_cor_DMR_sig_hyper);length(FK_cor_DMR_sig_hypo)
#22 28 22 26

#remove cor_DMRs_from_kids
#for MF_cor
hyper_no_cor_DMRs<-setdiff(DMR_hyper_kid, c(MK_cor_DMR_sig,FF_cor_DMR_sig))
hypo_no_cor_DMRs<-setdiff(DMR_hypo_kid, c(MK_cor_DMR_sig,FF_cor_DMR_sig))
length(hyper_no_cor_DMRs);length(hypo_no_cor_DMRs)

#reading parental related DMRs
#for mather
compare_name2<-"mother_AMA_vs_Young"
split_region2<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("D:/huahua/1.DMR_called/",compare_name2,"_",split_region2,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t",check.names = F)
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,];DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
DMR_hyper_mom<-unique(as.character(DMR_hyper$Row.names));DMR_hypo_mom<-unique(as.character(DMR_hypo$Row.names));length(DMR_hyper_mom);length(DMR_hypo_mom)
# 250 #269
#for father
compare_name2<-"father_AMA_vs_Young"
split_region2<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("D:/huahua/1.DMR_called/",compare_name2,"_",split_region2,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t",check.names = F)
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
DMR_hyper_dad<-unique(as.character(DMR_hyper$Row.names));DMR_hypo_dad<-unique(as.character(DMR_hypo$Row.names));length(DMR_hyper_dad);length(DMR_hypo_dad)
#217 #268
#find out pure Kid AMA
hyper_DMRs_kids_pure<-setdiff(hyper_no_cor_DMRs,c(DMR_hyper_mom,DMR_hyper_dad))
hypo_DMRs_kids_pure<-setdiff(hypo_no_cor_DMRs, c(DMR_hypo_mom,DMR_hypo_dad))
length(hyper_DMRs_kids_pure);length(hypo_DMRs_kids_pure)

##find mother related kids DMRs 
MK_related_hyper<-Reduce(intersect,c(list(DMR_hyper_kid),list(DMR_hyper_mom)))
MK_related_hypo<-Reduce(intersect,c(list(DMR_hypo_kid),list(DMR_hypo_mom)))
##find father related kids DMRs
FK_related_hyper<-Reduce(intersect,c(list(DMR_hyper_kid),list(DMR_hyper_dad)))
FK_related_hypo<-Reduce(intersect,c(list(DMR_hypo_kid),list(DMR_hypo_dad)))

##filter pure related DMRs but not cor DMRs
MK_related_hyper_pure_raw<-setdiff(MK_related_hyper,c(MK_cor_DMR_sig,FK_cor_DMR_sig))
MK_related_hypo_pure_raw<-setdiff(MK_related_hypo, c(MK_cor_DMR_sig,FK_cor_DMR_sig))
FK_related_hyper_pure_raw<-setdiff(FK_related_hyper, c(MK_cor_DMR_sig,FK_cor_DMR_sig))
FK_related_hypo_pure_raw <-setdiff(FK_related_hypo, c(MK_cor_DMR_sig,FK_cor_DMR_sig))
length(MK_related_hyper_pure_raw);length(MK_related_hypo_pure_raw);length(FK_related_hyper_pure_raw);length(FK_related_hypo_pure_raw);
#17;36;23;27

##find out both parental related kids DMRs 
FMK_common_hyper<-Reduce(intersect,c(list(DMR_hyper_kid),list(DMR_hyper_mom),list(DMR_hyper_dad)))
FMK_common_hypo<-Reduce(intersect,c(list(DMR_hypo_kid),list(DMR_hypo_mom),list(DMR_hypo_dad)))

FMK_common_hyper_no_cor<-setdiff(FMK_common_hyper,c(MK_cor_DMR_sig,FK_cor_DMR_sig))
FMK_common_hypo_no_cor<-setdiff(FMK_common_hypo,c(MK_cor_DMR_sig,FK_cor_DMR_sig))
length(FMK_common_hyper_no_cor);length(FMK_common_hypo_no_cor)
#1 4 
###find out mother or father related kids DMRs no common in both parents
MK_related_hyper_pure<-setdiff(MK_related_hyper_pure_raw,FMK_common_hyper_no_cor)
MK_related_hypo_pure<-setdiff(MK_related_hypo_pure_raw,FMK_common_hypo_no_cor)
FK_related_hyper_pure<-setdiff(FK_related_hyper_pure_raw,FMK_common_hyper_no_cor)
FK_related_hypo_pure <-setdiff(FK_related_hypo_pure_raw,FMK_common_hypo_no_cor)
length(MK_related_hyper_pure);length(MK_related_hypo_pure);length(FK_related_hyper_pure);length(FK_related_hypo_pure);
#16 32 22 23

##filter both related DMRs and cor DMRs
MK_hyper_high<-intersect(MK_related_hyper, MK_cor_DMR_sig)
MK_hypo_high <-intersect(MK_related_hypo, MK_cor_DMR_sig)
FK_hyper_high<-intersect(FK_related_hyper, FK_cor_DMR_sig)
FK_hypo_high <-intersect(FK_related_hypo, FK_cor_DMR_sig)
length(MK_hyper_high);length(MK_hypo_high);length(FK_hyper_high);length(FK_hypo_high);
#9;7;1;11
intersect(MK_hyper_high, FK_hyper_high)
intersect(MK_hypo_high, FK_hypo_high)
##0 ##0
intersect(MK_hyper_high, FK_hypo_high)
intersect(FK_hyper_high, MK_hypo_high)

#plot for correlated data 
father_data_analysis <-read.table(file="D:/huahua/3.correlated_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_father.txt",header=T,sep="\t")
head(father_data_analysis)
father_data_analysis_sig<-father_data_analysis[which(father_data_analysis$group == "candidate_region"),]
mother_data_analysis <-read.table(file="D:/huahua/3.correlated_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_mother.txt",header=T,sep="\t")
head(mother_data_analysis);
mother_data_analysis_sig<-mother_data_analysis[which(mother_data_analysis$group == "candidate_region"),]

#FK_hyper_high
FK_hyper_high_selected<- father_data_analysis_sig[which(father_data_analysis_sig$DMR_region %in% FK_hyper_high),]
#FK_hypo_high
FK_hypo_high_selected<- father_data_analysis_sig[which(father_data_analysis_sig$DMR_region %in% FK_hypo_high),]
#MK_hyper_high
MK_hyper_high_selected<- mother_data_analysis_sig[which(mother_data_analysis_sig$DMR_region %in% MK_hyper_high),]
#MK_hypo_high
MK_hypo_high_selected<- mother_data_analysis_sig[which(mother_data_analysis_sig$DMR_region %in% MK_hypo_high),]
plot_cordata<-rbind(FK_hyper_high_selected,FK_hypo_high_selected,MK_hyper_high_selected,MK_hypo_high_selected)
dim(plot_cordata)
plot_cordata$group<-c(rep("FK_hyper",nrow(FK_hyper_high_selected)),rep("FK_hypo",nrow(FK_hypo_high_selected)),
                      rep("MK_hyper",nrow(MK_hyper_high_selected)),rep("MK_hypo",nrow(MK_hypo_high_selected)))

kids_coldata<-kids_alldata[,c("Row.names","SYMBOL")]
colnames(kids_coldata)<-c("DMR_region","gene")
head(kids_coldata);head(plot_cordata)

plot_cordata2<-merge(plot_cordata,kids_coldata,by="DMR_region")
plot_cordata2$DMR_gene<-paste(plot_cordata2$DMR_region,plot_cordata2$gene,sep="_")
dim(plot_cordata);dim(plot_cordata2)
head(plot_cordata2)
plot_cordata2[which(plot_cordata2$DMR_region %in% FK_hyper_high),]
plot_cordata2$group<-factor(plot_cordata2$group,levels =c("FK_hyper","FK_hypo","MK_hyper","MK_hypo"))
plot_cordata2 <- plot_cordata2[order(plot_cordata2$group,decreasing = F),]
plot_cordata2$DMR_gene<-factor(plot_cordata2$DMR_gene,levels =as.character(unique((plot_cordata2$DMR_gene))))

head(plot_cordata2)
#formal plotting 
cor_plot0<-ggplot(data=plot_cordata2,aes(x=kids, y=parent))+geom_point()+
  stat_smooth(method="lm",se=FALSE,aes(colour=DMR_gene))+
  stat_cor(data=plot_cordata2, method = "spearman")+
  facet_wrap(~ DMR_gene, scales = "free",ncol = 7)+ NoLegend() 
cor_plot0
ggsave(file="D:/huahua/1.DMR_called/seperated_DMRs/cor_plot_for_seleced_both_related_and_cor_kids_AMA_DMR_seprated.pdf",cor_plot0,width = 80*7, height =70*ceiling(length(unique(plot_cordata$DMR_region))/7), units = "mm")

cor_plot1<-ggplot(data=plot_cordata2,aes(x=kids, y=parent))+geom_point()+
  stat_smooth(method="lm",se=FALSE)+
  stat_cor(data=plot_cordata2, method = "spearman")+
  facet_wrap(~ DMR_gene, scales = "free",ncol = 7)+ NoLegend() 
cor_plot1
cor_plot3<-ggplot(data=plot_cordata2,aes(x=kids, y=parent))+geom_point()+
  theme_bw()+
  stat_smooth(method="lm",se=FALSE)+
  stat_cor(data=plot_cordata2, method = "spearman")+
  facet_wrap(~ DMR_gene, scales = "free",ncol = 4)+ NoLegend() 
cor_plot3
ggsave(file="D:/huahua/1.DMR_called/seperated_DMRs/cor_plot_for_seleced_both_related_and_cor_kids_AMA_DMR_seprated1.pdf",cor_plot1,width = 80*7, height =70*ceiling(length(unique(plot_cordata$DMR_region))/7), units = "mm")
ggsave(file="D:/huahua/1.DMR_called/seperated_DMRs/cor_plot_for_seleced_both_related_and_cor_kids_AMA_DMR_seprated_trans_2.pdf",cor_plot3, height= 70*7, width =70*ceiling(length(unique(plot_cordata$DMR_region))/7), units = "mm")

cor_plot2<-ggscatter(plot_cordata2,x = "kids", y = "parent",
                     add = "reg.line",conf.int = TRUE,                          # Add confidence interval
                     color = "DMR_gene",# palette =ppCor_all2,           # Color by groups "cyl"
                     #facet.by = "DMR_gene",
                     #)+stat_cor(aes(color = gene_DMRs), label.x = 3)           # Add correlation coefficient
)+stat_cor(aes(color = DMR_gene), method = "spearman", label.x = 3)           # Add correlation coefficient
cor_plot2
cor_plot4<-cor_plot2+ facet_wrap(~ DMR_gene, scales = "free",ncol=7)+
  theme(legend.title = element_text(color = "black", size = 7),legend.position = "right")
cor_plot4
ggsave(file="D:/huahua/1.DMR_called/seperated_DMRs/cor_plot_for_seleced_both_related_and_cor_kids_AMA_DMR_seprated2.pdf",cor_plot4,width = 98*7, height =80*ceiling(length(unique(plot_cordata2$DMR_region))/7), units = "mm")
ggsave(file="D:/huahua/1.DMR_called/seperated_DMRs/cor_plot_for_seleced_both_related_and_cor_kids_AMA_DMR_seprated3.pdf",cor_plot2,width =400, height =450, units = "mm")


##analysis and enrichment analysis for  total ten list of DMRs
length(hyper_DMRs_kids_pure);length(hypo_DMRs_kids_pure)
length(MK_related_hyper_pure);length(MK_related_hypo_pure);
length(FK_related_hyper_pure);length(FK_related_hypo_pure);
length(MK_cor_DMR_sig_hyper);length(MK_cor_DMR_sig_hypo);
length(FK_cor_DMR_sig_hyper);length(FK_cor_DMR_sig_hypo)

length(FMK_common_hyper_no_cor);length(FMK_common_hypo_no_cor)
length(MF_cor_DMRs_hyper);length(MF_cor_DMRs_hypo)
unique(kids_alldata[which(kids_alldata$Row.names %in% FMK_common_hyper_no_cor),]$SYMBOL)
#"MIR4472-1"
unique(kids_alldata[which(kids_alldata$Row.names %in% FMK_common_hypo_no_cor),]$SYMBOL)
#"MIR6800"   "LINC02572" "RBMS3"     "FUBP3"   
unique(kids_alldata[which(kids_alldata$Row.names %in% MF_cor_DMRs_hyper),]$SYMBOL)
#"USP36"
unique(kids_alldata[which(kids_alldata$Row.names %in% MF_cor_DMRs_hypo),]$SYMBOL)
# "FAAP24"       "PCMT1"        "LOC101927450"



head(data_region_merge_DMR_all_big)
DMR_related_DEGs_list<-c(list(unique(kids_alldata[which(kids_alldata$Row.names %in% hyper_DMRs_kids_pure),]$SYMBOL)),
                   list(unique(kids_alldata[which(kids_alldata$Row.names %in% hypo_DMRs_kids_pure),]$SYMBOL)),
                   list(unique(kids_alldata[which(kids_alldata$Row.names %in% MK_related_hyper_pure),]$SYMBOL)),
                   list(unique(kids_alldata[which(kids_alldata$Row.names %in% MK_related_hypo_pure),]$SYMBOL)),
                   list(unique(kids_alldata[which(kids_alldata$Row.names %in% FK_related_hyper_pure),]$SYMBOL)),
                   list(unique(kids_alldata[which(kids_alldata$Row.names %in% FK_related_hypo_pure),]$SYMBOL)),
                   list(unique(kids_alldata[which(kids_alldata$Row.names %in% MK_cor_DMR_sig_hyper),]$SYMBOL)),
                   list(unique(kids_alldata[which(kids_alldata$Row.names %in% MK_cor_DMR_sig_hypo),]$SYMBOL)),
                   list(unique(kids_alldata[which(kids_alldata$Row.names %in% FK_cor_DMR_sig_hyper),]$SYMBOL)),
                   list(unique(kids_alldata[which(kids_alldata$Row.names %in% FK_cor_DMR_sig_hypo),]$SYMBOL)))
names(DMR_related_DEGs_list)<-c("K-S_hyper_DMRs","K-S_hypo_DMRs","M_related_hyper_DMRs","M_related_hypo_DMRs",
                                "F_related_hyper_DMRs","F_related_hypo_DMRs","M_cor_hyper_DMRs","M_cor_hypo_DMRs","F_cor_hyper_DMRs","F_cor_hypo_DMRs")
saveRDS(DMR_related_DEGs_list, file = "D:/huahua/1.DMR_called/seperated_DMRs/ten_list_DEGs.rds")

#enrichment
for ( DEGs_name in names(DMR_related_DEGs_list)){
  #DEGs_name <- "F_related_hyper_DMRs"  ##test line
  DEGs_list<-DMR_related_DEGs_list[[DEGs_name]]
  print(DEGs_name)
  DEGs_list3<- bitr(DEGs_list, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) #translate into other types ID
  BP <- enrichGO(DEGs_list3$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2, readable=T)  
  MF <- enrichGO(DEGs_list3$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, minGSSize = 2, readable=T) 
  CC <- enrichGO(DEGs_list3$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2,  readable=T) 
  head(summary(BP))
  head(as.data.frame(BP_simp@result))
  kk <- clusterProfiler::enrichKEGG(gene = DEGs_list3$ENTREZID,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
  kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
  head(as.data.frame(kk@result))
  ## Reactome pathway enrichment analysis
  rpea <- enrichPathway(gene= DEGs_list3$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
  head(as.data.frame(rpea@result))
  write.csv(as.data.frame(CC@result), file=paste0("D:/huahua/1.DMR_called/seperated_DMRs/",DEGs_name,"_15q005_CC_GO.csv")) 
  write.csv(as.data.frame(BP@result),file=paste0("D:/huahua/1.DMR_called/seperated_DMRs/",DEGs_name,"_15q005_BP_GO.csv")) 
  write.csv(as.data.frame(MF@result),file=paste0("D:/huahua/1.DMR_called/seperated_DMRs/",DEGs_name,"_15q005_MF_GO.csv")) 
  write.csv((as.data.frame(kk@result)), file=paste0("D:/huahua/1.DMR_called/seperated_DMRs/",DEGs_name,"_15q005_KEGG.csv")) 
  write.csv((as.data.frame(rpea@result)),file=paste0("D:/huahua/1.DMR_called/seperated_DMRs/",DEGs_name,"_15q005_RPEA.csv"))  
  gene_GO_RPEA_erichment_results<-list(CC,BP,MF,kk,rpea)
  names(gene_GO_RPEA_erichment_results)<-c("CC","BP","MF","kkeg","rpea")
  saveRDS(gene_GO_RPEA_erichment_results, file = paste("D:/huahua/1.DMR_called/seperated_DMRs/",DEGs_name,"_15q005_GO_RPEA_result.rds",sep = ""))
  }  

##plot enrichment 
#plot for changed proteins
#construct dataframe for all GO terms
sample_names<-c("K-S_hyper_DMRs","K-S_hypo_DMRs","M_related_hyper_DMRs","M_related_hypo_DMRs", "F_related_hyper_DMRs",
                "F_related_hypo_DMRs","M_cor_hyper_DMRs","M_cor_hypo_DMRs","F_cor_hyper_DMRs","F_cor_hypo_DMRs")
GO_lists<-list()
for ( DEGs_name in sample_names){
  #DEGs_name<-"F_related_hypo_DMRs"
  print(DEGs_name)
  GO_table<-read.csv(paste0("D:/huahua/1.DMR_called/seperated_DMRs/",DEGs_name,"_15q005_BP_GO.csv")) 
  # head(GO_table);dim(GO_table)
  GO_lists<-c(GO_lists,list(GO_table))
}
length(sample_names);length(GO_lists)
names(GO_lists)<-sample_names


dim(GO_lists[[1]][which(GO_lists[[1]]$pvalue <0.05 & GO_lists[[1]]$Count>1),])#64 
dim(GO_lists[[2]][which(GO_lists[[2]]$pvalue <0.05 & GO_lists[[2]]$Count>1),])#64 
dim(GO_lists[[3]][which(GO_lists[[3]]$pvalue <0.05 & GO_lists[[3]]$Count>1),])#28 
dim(GO_lists[[4]][which(GO_lists[[4]]$pvalue <0.05 & GO_lists[[4]]$Count>1),])#28 
dim(GO_lists[[5]][which(GO_lists[[5]]$pvalue <0.05 & GO_lists[[5]]$Count>1),])#6 
dim(GO_lists[[6]][which(GO_lists[[6]]$pvalue <0.05 & GO_lists[[6]]$Count>1),])#10 
dim(GO_lists[[7]][which(GO_lists[[7]]$pvalue <0.05 & GO_lists[[7]]$Count>1),])#4
dim(GO_lists[[8]][which(GO_lists[[8]]$pvalue <0.05 & GO_lists[[8]]$Count>1),])#73
dim(GO_lists[[9]][which(GO_lists[[9]]$pvalue <0.05 & GO_lists[[9]]$Count>1),])#44
dim(GO_lists[[10]][which(GO_lists[[10]]$pvalue <0.05 & GO_lists[[10]]$Count>1),])#15

Enrich_all_plot<-as.data.frame(rbind(GO_lists[[1]],GO_lists[[2]],GO_lists[[3]],GO_lists[[4]],GO_lists[[5]],
                                     GO_lists[[6]],GO_lists[[7]],GO_lists[[8]],GO_lists[[9]],GO_lists[[10]]))
Enrich_all_plot$group<-c(rep("K-S_hyper_DMRs",nrow(GO_lists[[1]])),rep("K-S_hypo_DMRs",nrow(GO_lists[[2]])),
                         rep("M_related_hyper_DMRs",nrow(GO_lists[[3]])),rep("M_related_hypo_DMRs",nrow(GO_lists[[4]])),
                         rep("F_related_hyper_DMRs",nrow(GO_lists[[5]])),rep("F_related_hypo_DMRs",nrow(GO_lists[[6]])),
                         rep("M_cor_hyper_DMRs",nrow(GO_lists[[7]])),rep("M_cor_hypo_DMRs",nrow(GO_lists[[8]])),
                         rep("F_cor_hyper_DMRs",nrow(GO_lists[[9]])),rep("F_cor_hypo_DMRs",nrow(GO_lists[[10]]))
)

dim(Enrich_all_plot)#6184   11
Enrich_all_plot$Log_Padjust<- c(-log10(Enrich_all_plot$p.adjust))
range(Enrich_all_plot$Log_Padjust)# 0.08283877 2.60992824
Enrich_all_plot$Log_pvalue<- c(-log10(Enrich_all_plot$pvalue))
range(Enrich_all_plot$Log_pvalue)# 0.08283877 5.23317753
range(Enrich_all_plot$Count)#1 6

Enrich_all_plot_sig_p005_more_than2<-Enrich_all_plot[which(Enrich_all_plot$pvalue<0.05 & Enrich_all_plot$Count>1),]
dim(Enrich_all_plot_sig_p005_more_than2)#336  13
Enrich_all_number <- data.frame(table(as.character(Enrich_all_plot_sig_p005_more_than2$Description)))
Enrich_all_number2 <- Enrich_all_number[which(Enrich_all_number$Freq>1),]
dim(Enrich_all_number2)# 28    2
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot2<-merge(Enrich_all_plot_sig_p005_more_than2,Enrich_all_number)
dim(Enrich_all_plot2)#336  14
write.table(Enrich_all_plot2, file="D:/huahua/1.DMR_called/seperated_DMRs/GO_BP_merge.txt",row.names=T, col.names=T) 

#for no selected GO terms of change protein
Enrich_all_plot_noselect<-read.table("D:/huahua/1.DMR_called/seperated_DMRs/GO_BP_merge.txt",header =T) 
head(Enrich_all_plot_noselect)
length(unique(Enrich_all_plot_noselect$Description))#307 
length(unique(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>1),]$Description))#307 
Enrich_all_plot_noselect<-Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$freq,decreasing = F),]
Enrich_all_plot_noselect$group<-factor(Enrich_all_plot_noselect$group,levels =sample_names)
Enrich_all_plot_noselect <- Enrich_all_plot_noselect[order(Enrich_all_plot_noselect$group),]
Enrich_all_plot_noselect$Description<-factor(Enrich_all_plot_noselect$Description,levels = c(unique(as.character(Enrich_all_plot_noselect$Description))))

range(Enrich_all_plot_noselect$Log_pvalue) #1.304112 5.233178
range(Enrich_all_plot_noselect$Log_Padjust) #0.4537519 2.6099282
range(Enrich_all_plot_noselect$Count) #2 6
range(Enrich_all_plot_noselect$freq) #1 3

#data_for_plot<-Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>2),]
data_for_plot<-Enrich_all_plot_noselect
data_for_plot$group<-factor(data_for_plot$group,levels =sample_names)
data_for_plot <- data_for_plot[order(data_for_plot$group),]
data_for_plot <- Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$freq>1),]
data_for_plot<-data_for_plot[order(data_for_plot$freq,decreasing = F),]
head(data_for_plot)
data_for_plot$Description<-factor(data_for_plot$Description,levels = c(unique(as.character(data_for_plot$Description))))

ggplot(data_for_plot,aes(x=group,y=Description,size=Count,colour=Log_pvalue))+
  geom_point(alpha =0.8,na.rm = TRUE)+
  scale_size(breaks = c(2,3,4,5),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(5,'RdYlBu')),breaks=c(0,1.5,2,2.5,3,3.5),name='-log10(p_value)')+  
  scale_y_discrete(labels=function(x) str_wrap(x, width=150))+
  theme_classic()+labs(x="",y="GO terms",title="Kids AMA related change protein BP enrichment")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black',vjust=0,hjust=1),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 10),
        legend.position ="right",legend.direction = "vertical")
write.table(data_for_plot, file="D:/huahua/1.DMR_called/seperated_DMRs/GO_BP_merge_appear_more_than_two.txt",row.names=T, col.names=T) 

data_for_plot<-Enrich_all_plot_noselect
data_for_plot$group<-factor(data_for_plot$group,levels =sample_names)
data_for_plot <- data_for_plot[order(data_for_plot$group),]
data_for_plot<-data_for_plot[order(data_for_plot$freq,decreasing = F),]
head(data_for_plot)
data_for_plot$Description<-factor(data_for_plot$Description,levels = c(unique(as.character(data_for_plot$Description))))

ggplot(data_for_plot,aes(x=group,y=Description,size=Count,colour=Log_pvalue))+
  geom_point(alpha =0.8,na.rm = TRUE)+
  scale_size(breaks = c(2,3,4,5),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(5,'RdYlBu')),breaks=c(0,1.5,3,4.5,6),name='-log10(p_value)')+  
  scale_y_discrete(labels=function(x) str_wrap(x, width=150))+
  theme_classic()+labs(x="",y="GO terms",title="Kids AMA related change protein BP enrichment")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 6,colour = 'black',vjust=0.5,hjust=1),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 10),
        legend.position ="right",legend.direction = "vertical")
write.table(data_for_plot, file="D:/huahua/1.DMR_called/seperated_DMRs/GO_BP_merge_appear_all_wait_selected.txt",row.names=T, col.names=T) 

