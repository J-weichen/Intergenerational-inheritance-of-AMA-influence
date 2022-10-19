rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)#org.Hs.eg.db 的数据类型以及使用简介：可用于数据类型转换
keytypes(org.Hs.eg.db)#查看org.Hs.eg.db数据对象里面包含着各大主流数据库的数据
library(RColorBrewer)

Mom_kid_common_DEGsdata<-read.table("/mnt/data/chenwei/huahua/5.valificat_result/AMA_common_DEG_between_kids_and_mother_trend_add.txt",head = T,sep="\t")
head(Mom_kid_common_DEGsdata)
Mom_kid_common_up<-na.omit(unique(Mom_kid_common_DEGsdata[which(Mom_kid_common_DEGsdata$trend == "common_UP_DEGs"),]$MK_DEGs))
Mom_kid_common_down<-na.omit(unique(Mom_kid_common_DEGsdata[which(Mom_kid_common_DEGsdata$trend == "common_Down_DEGs"),]$MK_DEGs))
gene_list<-c(list(Mom_kid_common_up),list(Mom_kid_common_down))
names(gene_list)<-c("Mom_kid_common_up","Mom_kid_common_down")
#list.files命令将input文件夹下所有文件名输入
for ( fname in c("Mom_kid_common_up","Mom_kid_common_down")){
  # fname<-"Mom_kid_common_up"     
  gene_list0<-gene_list[[fname]]
  print(fname)
  gene_list1<- bitr(gene_list0, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) 
  BP <- enrichGO(gene_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2, readable=T)  
  MF <- enrichGO(gene_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, minGSSize = 2, readable=T) 
  CC <- enrichGO(gene_list1$ENTREZID,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 2,  readable=T) 
  head(summary(BP))
  CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
  BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
  MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
  head(as.data.frame(BP_simp@result))
  kk <- clusterProfiler::enrichKEGG(gene = gene_list1$ENTREZID,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
  kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
  head(as.data.frame(kk@result))
  ## Reactome pathway enrichment analysis
  rpea <- enrichPathway(gene= gene_list1$ENTREZID,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
  head(as.data.frame(rpea@result))
  
  write.table(as.data.frame(CC@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/",fname,"_unselect_CC_GO.txt",sep = ""), row.names=T, col.names=T) 
  write.table(as.data.frame(BP@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/",fname,"_unselect_BP_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(MF@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/",fname,"_unselect_MF_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(CC_simp@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/",fname,"_unselect_CC_GO_simp.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(BP_simp@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/",fname,"_unselect_BP_GO_simp.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(MF_simp@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/",fname,"_unselect_MF_GO_simp.txt",sep = ""), row.names=T, col.names=T) 
  write.table((as.data.frame(kk@result)), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/",fname,"_unselect_KEGG.txt",sep = ""),row.names=T, col.names=T) 
  write.table((as.data.frame(rpea@result)), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/",fname,"_unselect_RPEA.txt",sep = ""),row.names=T, col.names=T) 
  gene_GO_RPEA_erichment_results<-list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
  names(gene_GO_RPEA_erichment_results)<-c("CC","BP","MF","CC_simp","BP_simp","MF_simp","kkeg","rpea")
  saveRDS(gene_GO_RPEA_erichment_results, file = paste("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/",fname,"_q005_unselect_GO_RPEA_result.rds",sep = ""))
}

##plot enrichment
#construct dataframe for all GO terms 
MK_common_up_BP <-read.table("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/Mom_kid_common_up_unselect_BP_GO.txt",head = T)
MK_common_down_BP <-read.table("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/Mom_kid_common_down_unselect_BP_GO.txt",head = T)

dim(MK_common_up_BP[which(MK_common_up_BP$Count>=3),]);dim(MK_common_down_BP[which(MK_common_down_BP$Count>=3),])#145  9  126  9
Enrich_all_plot<-as.data.frame(rbind(MK_common_up_BP,MK_common_down_BP))
Enrich_all_plot$group<-c(rep("MK_common_up",nrow(MK_common_up_BP)),rep("MK_common_down",nrow(MK_common_down_BP)))
dim(Enrich_all_plot)#2644   10
head(Enrich_all_plot)

Enrich_all_plot2<-Enrich_all_plot[which(Enrich_all_plot$pvalue<0.05 & Enrich_all_plot$Count>=2),];dim(Enrich_all_plot2)#340  10
Enrich_all_plot2$Log_pvalue<- c(-log10(Enrich_all_plot2$pvalue))
range(Enrich_all_plot2$Log_pvalue)# 1.302200 5.068468
range(Enrich_all_plot2$Count)# 2 8

Enrich_all_number<-data.frame(table(as.character(Enrich_all_plot2$Description)))
#Enrich_all_number2<-Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot2<-merge(Enrich_all_plot2,Enrich_all_number)
write.table(Enrich_all_plot2, file="/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/pvalue_0.05_Count_2_GO_BP_merge.txt",row.names=T, col.names=T) 

#no selected 
Enrich_all_plot_noselect<-read.table("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/pvalue_0.05_Count_2_GO_BP_merge.txt",header =T) 
head(Enrich_all_plot_noselect)
length(unique(as.character(Enrich_all_plot_noselect$Description)))#318
length(unique(as.character(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>=3),]$Description)))#138
range(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>=3),]$Count)# 3 8
Enrich_all_plot_noselect2<-Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count >= 3),]
Enrich_all_number <-data.frame(table(as.character(Enrich_all_plot_noselect2$Description)))
colnames(Enrich_all_number)<-c("Description","freq2")
Enrich_all_plot2<-merge(Enrich_all_plot_noselect2,Enrich_all_number)
head(Enrich_all_plot2)
length(unique(as.character(Enrich_all_plot2$Description)))#138
write.table(Enrich_all_plot2, file="/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/GO_BP_merge_count_more_than_three.txt",row.names=T, col.names=T) 


Enrich_all_plot_noselect$group<-factor(Enrich_all_plot_noselect$group,levels = c("MK_common_up","MK_common_down"))
Enrich_all_plot_noselect2<-Enrich_all_plot_noselect2[order(Enrich_all_plot_noselect2$freq,Enrich_all_plot_noselect2$group,decreasing = T ),]
Enrich_all_plot_noselect2$Description<-factor(Enrich_all_plot_noselect2$Description,levels = rev(c(unique(as.character(Enrich_all_plot_noselect2$Description)))))
length(unique(as.character(Enrich_all_plot_noselect2[which(Enrich_all_plot_noselect2$freq>1),]$Description)))#16 

plot_BP<-ggplot(Enrich_all_plot_noselect2, aes(x=group,y=Description,size=Count,colour=Log_pvalue))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,2,3,4,5,6,7,8),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,1,2,3,4,5),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="AMA related mom and kids commmon DEGs: BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/Enrich_all_BP_noselect.pdf",plot_BP,width = 12, height =30)
write.table(Enrich_all_plot_noselect2, file="/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/GO_BP_merge_count_more_than_two_wait_selected.txt",row.names=T, col.names=T) 

#for selected  
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/common_MK_DEGs_GO_BP_merge_count_more_than_two_selected_one.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
Enrich_all_plot_select$group<-factor(Enrich_all_plot_select$group,levels = c("MK_common_up","MK_common_down"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$freq,Enrich_all_plot_select$group,decreasing = T ),]
Enrich_all_plot_select$Description<-factor(Enrich_all_plot_select$Description,levels = rev(c(unique(as.character(Enrich_all_plot_select$Description)))))
length(unique(as.character(Enrich_all_plot_select[which(Enrich_all_plot_select$freq>1),]$Description)))#6

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=group,y=Description,size=Count,colour=Log_pvalue))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,2,3,4,5,6,7,8),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,1,2,3,4,5),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="AMA related mom and kids commmon DEGs: BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 8,colour = 'black',vjust=1,hjust=1,angle = 90),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/common_MK_DEGs_Enrich_all_BP_select_one.pdf",plot_BP,width = 10, height =10)
write.table(Enrich_all_plot_select, file="/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/common_MK_DEGs_GO_BP_merge_count_more_than_two_selected_one_wait.txt",row.names=T, col.names=T) 

#for selected  final 
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/common_MK_DEGs_GO_BP_merge_count_more_than_three_selected_two.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
Enrich_all_plot_select$group<-factor(Enrich_all_plot_select$group,levels = c("MK_common_up","MK_common_down"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$freq,Enrich_all_plot_select$group,decreasing = T ),]
Enrich_all_plot_select$Description<-factor(Enrich_all_plot_select$Description,levels = rev(c(unique(as.character(Enrich_all_plot_select$Description)))))
length(unique(as.character(Enrich_all_plot_select[which(Enrich_all_plot_select$freq>1),]$Description)))#5

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=group,y=Description,size=Count,colour=Log_pvalue))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,2,3,4,5,6,7,8),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,1,2,3,4,5),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))+
  theme_classic()+labs(x="",y="GO terms",title="AMA related mom and kids commmon DEGs: BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 8,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/common_MK_DEGs_Enrich_all_BP_select_final.pdf",plot_BP,width = 6, height =6)
write.table(Enrich_all_plot_select, file="/mnt/data/chenwei/huahua/5.valificat_result/common_DEGs_enrichment/common_MK_DEGs_GO_BP_merge_count_more_than_two_selected_final.txt",row.names=T, col.names=T) 
