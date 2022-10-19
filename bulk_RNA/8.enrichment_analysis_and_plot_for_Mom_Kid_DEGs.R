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


#FOR Mom
tag_name<-"Mom"; sampleA <-"AMA";sampleB <-"Young"
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file13_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
Mom_AMA_DEGs<-read.table(file)
Mom_AMA_DEGs_up_1.5<-Mom_AMA_DEGs[which(Mom_AMA_DEGs$log2FoldChange>0),]$ID
Mom_AMA_DEGs_down_1.5<-Mom_AMA_DEGs[which(Mom_AMA_DEGs$log2FoldChange<0),]$ID
#FOR kids
#read expression data 
tag_name<-"Kid"; sampleA <-"AMA";sampleB <-"Young" 
file <- paste0("/mnt/data/chenwei/huahua/5.valificat_result/file13_all_filter_",tag_name,"_verification_",sampleA,"_vs_",sampleB,".DEG_information_pvalue005_FC1.5.txt")
kids_AMA_DEGs<-read.table(file)
kids_AMA_DEGs_up_1.5<-kids_AMA_DEGs[which(kids_AMA_DEGs$log2FoldChange>0),]$ID
kids_AMA_DEGs_down_1.5<-kids_AMA_DEGs[which(kids_AMA_DEGs$log2FoldChange<0),]$ID
length(Mom_AMA_DEGs_up_1.5);length(Mom_AMA_DEGs_down_1.5);length(kids_AMA_DEGs_up_1.5);length(kids_AMA_DEGs_down_1.5)

gene_list<-c(list(Mom_AMA_DEGs_up_1.5),list(Mom_AMA_DEGs_down_1.5),list(kids_AMA_DEGs_up_1.5),list(kids_AMA_DEGs_down_1.5))
names(gene_list)<-c("Mom_up","Mom_down","Kid_up","Kid_down")
#list.files命令将input文件夹下所有文件名输入
for ( fname in c("Mom_up","Mom_down","Kid_up","Kid_down")){
  # fname<-"Mom_up"     
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
  
  write.table(as.data.frame(CC@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/",fname,"_unselect_CC_GO.txt",sep = ""), row.names=T, col.names=T) 
  write.table(as.data.frame(BP@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/",fname,"_unselect_BP_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(MF@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/",fname,"_unselect_MF_GO.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(CC_simp@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/",fname,"_unselect_CC_GO_simp.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(BP_simp@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/",fname,"_unselect_BP_GO_simp.txt",sep = ""),row.names=T, col.names=T) 
  write.table(as.data.frame(MF_simp@result), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/",fname,"_unselect_MF_GO_simp.txt",sep = ""), row.names=T, col.names=T) 
  write.table((as.data.frame(kk@result)), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/",fname,"_unselect_KEGG.txt",sep = ""),row.names=T, col.names=T) 
  write.table((as.data.frame(rpea@result)), file=paste("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/",fname,"_unselect_RPEA.txt",sep = ""),row.names=T, col.names=T) 
  gene_GO_RPEA_erichment_results<-list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
  names(gene_GO_RPEA_erichment_results)<-c("CC","BP","MF","CC_simp","BP_simp","MF_simp","kkeg","rpea")
  saveRDS(gene_GO_RPEA_erichment_results, file = paste("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/",fname,"_q005_unselect_GO_RPEA_result.rds",sep = ""))
}


##plot enrichment
#construct dataframe for all GO terms 
Mom_up_BP <-read.table("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_up_unselect_BP_GO_simp.txt",head = T)
Mom_down_BP <-read.table("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_down_unselect_BP_GO_simp.txt",head = T)
Kid_up_BP <-read.table("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Kid_up_unselect_BP_GO_simp.txt",head = T)
Kid_down_BP <-read.table("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Kid_down_unselect_BP_GO_simp.txt",head = T)

dim(Mom_up_BP[which(Mom_up_BP$Count>=3),]);dim(Mom_down_BP[which(Mom_down_BP$Count>=3),])#37  9  47  9
dim(Kid_up_BP[which(Kid_up_BP$Count>=3),]);dim(Kid_down_BP[which(Kid_down_BP$Count>=3),])#86  9  141  9

Enrich_all_plot<-as.data.frame(rbind(Mom_up_BP,Mom_down_BP,Kid_up_BP,Kid_down_BP))
Enrich_all_plot$group<-c(rep("Mom_up",nrow(Mom_up_BP)),rep("Mom_down",nrow(Mom_down_BP)),rep("Kid_up",nrow(Kid_up_BP)),rep("Kid_down",nrow(Kid_down_BP)))
dim(Enrich_all_plot)#94 10
head(Enrich_all_plot)

Enrich_all_plot2<-Enrich_all_plot[which(Enrich_all_plot$pvalue<0.05 & Enrich_all_plot$Count>=2),];dim(Enrich_all_plot2)#329  10
Enrich_all_plot2$Log_pvalue<- c(-log10(Enrich_all_plot2$pvalue))
range(Enrich_all_plot2$Log_pvalue)# 2.54490 74.11978
range(Enrich_all_plot2$Count)#  2 107

Enrich_all_number<-data.frame(table(as.character(Enrich_all_plot2$Description)))
#Enrich_all_number2<-Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot2<-merge(Enrich_all_plot2,Enrich_all_number)
write.table(Enrich_all_plot2, file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_2_GO_BP_merge.txt",row.names=T, col.names=T) 

#no selected 
Enrich_all_plot_noselect<-read.table("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_2_GO_BP_merge.txt",header =T) 
head(Enrich_all_plot_noselect)
length(unique(as.character(Enrich_all_plot_noselect$Description)))#311
length(unique(as.character(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>=3),]$Description)))#293
range(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>=3),]$Count)#  3 107
Enrich_all_plot_noselect2<-Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count >= 3),]
Enrich_all_number <-data.frame(table(as.character(Enrich_all_plot_noselect2$Description)))
colnames(Enrich_all_number)<-c("Description","freq2")
Enrich_all_plot2<-merge(Enrich_all_plot_noselect2,Enrich_all_number)
head(Enrich_all_plot2)
length(unique(as.character(Enrich_all_plot2$Description)))#293
write.table(Enrich_all_plot2, file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_GO_BP_merge.txt",row.names=T, col.names=T) 


Enrich_all_plot_noselect$group<-factor(Enrich_all_plot_noselect$group,levels = c("Kid_down","Kid_up","Mom_down","Mom_up"))
Enrich_all_plot_noselect2<-Enrich_all_plot_noselect2[order(Enrich_all_plot_noselect2$freq,Enrich_all_plot_noselect2$group,decreasing = T ),]
Enrich_all_plot_noselect2$Description<-factor(Enrich_all_plot_noselect2$Description,levels = rev(c(unique(as.character(Enrich_all_plot_noselect2$Description)))))
length(unique(as.character(Enrich_all_plot_noselect2[which(Enrich_all_plot_noselect2$freq>1),]$Description))) # 11

range(Enrich_all_plot_noselect2$Log_pvalue)#2.54490 74.11978
range(Enrich_all_plot_noselect2$Count)     # 3 107
Enrich_all_plot_noselect2$Log_pvalue2<-Enrich_all_plot_noselect2$Log_pvalue
Enrich_all_plot_noselect2[which(Enrich_all_plot_noselect2$Log_pvalue2 >10),]$Log_pvalue2<-10
plot_BP<-ggplot(Enrich_all_plot_noselect2, aes(x=group,y=Description,size=Count,colour=Log_pvalue2))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(0,10,20,40,80,120),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,2,4,6,8,10),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="AMA related DEGs in mom or  kids: BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_noselect.pdf",plot_BP,width = 12, height =35)
write.table(Enrich_all_plot_noselect2, file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_wait_selected.txt",row.names=T, col.names=T) 


#for selected  
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_selected_one.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
length(unique(Enrich_all_plot_select$Description))#196

Enrich_all_plot_select$group<-factor(Enrich_all_plot_select$group,levels = c("Kid_up","Kid_down","Mom_down","Mom_up"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$freq,Enrich_all_plot_select$group,decreasing = T ),]
Enrich_all_plot_select$Description<-factor(Enrich_all_plot_select$Description,levels = rev(c(unique(as.character(Enrich_all_plot_select$Description)))))
length(unique(as.character(Enrich_all_plot_select[which(Enrich_all_plot_select$freq>1),]$Description)))#8

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=group,y=Description,size=Count,colour=Log_pvalue2))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(3,6,9,12,15,18),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,4,8,12),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="Abortion related ICP gene  BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_select_one.pdf",plot_BP,width = 12, height =30)
write.table(Enrich_all_plot_select, file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_select_one_wait.txt",row.names=T, col.names=T) 

#for selected  two
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_select_two.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
length(unique(Enrich_all_plot_select$Description))#196

Enrich_all_plot_select$group<-factor(Enrich_all_plot_select$group,levels = c("Kid_up","Kid_down","Mom_down","Mom_up"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$freq,Enrich_all_plot_select$group,decreasing = T ),]
Enrich_all_plot_select$Description<-factor(Enrich_all_plot_select$Description,levels = rev(c(unique(as.character(Enrich_all_plot_select$Description)))))
length(unique(as.character(Enrich_all_plot_select[which(Enrich_all_plot_select$freq>1),]$Description)))#8

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=group,y=Description,size=Count,colour=Log_pvalue2))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(3,6,9,12,15,18),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,3,6,9,12),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="Abortion related ICP gene  BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=1,hjust=1,angle = 60),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_select_two.pdf",plot_BP,width = 12, height =24)
write.table(Enrich_all_plot_select, file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_select_two_wait.txt",row.names=T, col.names=T) 

#for selected  final
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_three_Enrich_all_BP_select_final.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
length(unique(Enrich_all_plot_select$Description))#98

Enrich_all_plot_select$group<-factor(Enrich_all_plot_select$group,levels = c("Kid_up","Kid_down","Mom_up","Mom_down"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$freq,Enrich_all_plot_select$group,decreasing = T ),]
Enrich_all_plot_select$Description<-factor(Enrich_all_plot_select$Description,levels = rev(c(unique(as.character(Enrich_all_plot_select$Description)))))
length(unique(as.character(Enrich_all_plot_select[which(Enrich_all_plot_select$freq>1),]$Description)))#8

range(Enrich_all_plot_select$Log_pvalue2)# 2.5449 10.0000
range(Enrich_all_plot_select$Count)# 3 107

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=group,y=Description,size=Count,colour=Log_pvalue2))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(5,10,20,40,60,80,100,120),range = c(1,8),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,2,4,6,8,10),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
  theme_classic()+labs(x="",y="GO terms",title="Kid or mom related DEGs:: BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_final.pdf",plot_BP,width = 9, height =15)
write.table(Enrich_all_plot_select, file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_select_final.txt",row.names=T, col.names=T) 

#for selected  final 2
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_select_final2.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
length(unique(Enrich_all_plot_select$Description))#76

Enrich_all_plot_select$group<-factor(Enrich_all_plot_select$group,levels = c("Kid_up","Kid_down","Mom_up","Mom_down"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$freq,Enrich_all_plot_select$group,decreasing = T ),]
Enrich_all_plot_select$Description<-factor(Enrich_all_plot_select$Description,levels = rev(c(unique(as.character(Enrich_all_plot_select$Description)))))
length(unique(as.character(Enrich_all_plot_select[which(Enrich_all_plot_select$freq>1),]$Description)))#8

range(Enrich_all_plot_select$Log_pvalue2)# 2.639675 10.0000
range(Enrich_all_plot_select$Count)# 3 107

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=group,y=Description,size=Count,colour=Log_pvalue2))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(5,10,20,40,60,80,100,120),range = c(1,8),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,2,4,6,8,10),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))+
  theme_classic()+labs(x="",y="GO terms",title="Kid or mom related DEGs:: BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 12,colour = 'black',vjust=0.5,hjust=1,angle = 90),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_final2.pdf",plot_BP,width = 8, height =12)
write.table(Enrich_all_plot_select, file="/mnt/data/chenwei/huahua/5.valificat_result/mom_kid_DEG_enrichment/Mom_or_kid_pvalue_0.05_Count_more_than_two_Enrich_all_BP_final2.txt",row.names=T, col.names=T) 
