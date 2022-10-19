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

##plot enrichment for DMRs of kids

Kid_hyper_BP <-read.csv("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/kids_AMA_vs_Young_hyper_CpGs_200bp_15q005_unselect_BP_GO.csv",head = T)
Kid_hypo_BP <-read.csv("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/kids_AMA_vs_Young_hypo_CpGs_200bp_15q005_unselect_BP_GO.csv",head = T)

dim(Kid_hyper_BP[which(Kid_hyper_BP$Count>=3),]);dim(Kid_hypo_BP[which(Kid_hypo_BP$Count>=3),])#145  9  126  9
Enrich_all_plot<-as.data.frame(rbind(Kid_hyper_BP,Kid_hypo_BP))
Enrich_all_plot$group<-c(rep("Kid_hyper",nrow(Kid_hyper_BP)),rep("Kid_hypo",nrow(Kid_hypo_BP)))
dim(Enrich_all_plot)#4212   11
head(Enrich_all_plot)

Enrich_all_plot2<-Enrich_all_plot[which(Enrich_all_plot$pvalue<0.05 & Enrich_all_plot$Count>=2),];dim(Enrich_all_plot2)#270  11
Enrich_all_plot2$Log_pvalue<- c(-log10(Enrich_all_plot2$pvalue))
range(Enrich_all_plot2$Log_pvalue)# 1.302895 4.096528
range(Enrich_all_plot2$Count)# 2 11

Enrich_all_number<-data.frame(table(as.character(Enrich_all_plot2$Description)))
#Enrich_all_number2<-Enrich_all_number[which(Enrich_all_number$Freq>1),]
colnames(Enrich_all_number)<-c("Description","freq")
Enrich_all_plot2<-merge(Enrich_all_plot2,Enrich_all_number)
write.table(Enrich_all_plot2, file="/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_pvalue_0.05_Count_2_GO_BP_merge.txt",row.names=T, col.names=T) 

#no selected 
Enrich_all_plot_noselect<-read.table("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_pvalue_0.05_Count_2_GO_BP_merge.txt",header =T) 
head(Enrich_all_plot_noselect)
length(unique(as.character(Enrich_all_plot_noselect$Description)))#265
length(unique(as.character(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>=3),]$Description)))# 160
range(Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count>=3),]$Count)# 3 11
Enrich_all_plot_noselect2<-Enrich_all_plot_noselect[which(Enrich_all_plot_noselect$Count >= 3),]
Enrich_all_number <-data.frame(table(as.character(Enrich_all_plot_noselect2$Description)))
colnames(Enrich_all_number)<-c("Description","freq2")
Enrich_all_plot2<-merge(Enrich_all_plot_noselect2,Enrich_all_number)
head(Enrich_all_plot2)
length(unique(as.character(Enrich_all_plot2$Description)))#160
write.table(Enrich_all_plot2, file="/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_GO_BP_merge_count_more_than_two.txt",row.names=T, col.names=T) 


Enrich_all_plot_noselect$group<-factor(Enrich_all_plot_noselect$group,levels = c("Kid_hyper","Kid_hypo"))
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
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_Enrich_DMRs_all_BP_noselect.pdf",plot_BP,width = 12, height =30)
write.table(Enrich_all_plot_noselect2, file="/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_DMRs_GO_BP_merge_count_more_than_two_wait_selected.txt",row.names=T, col.names=T) 

#for selected  
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_DMRs_GO_BP_merge_count_more_than_two_selected_one.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
Enrich_all_plot_select$group<-factor(Enrich_all_plot_select$group,levels = c("Kid_hyper","Kid_hypo"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$freq,Enrich_all_plot_select$group,decreasing = T ),]
Enrich_all_plot_select$Description<-factor(Enrich_all_plot_select$Description,levels = rev(c(unique(as.character(Enrich_all_plot_select$Description)))))
length(unique(as.character(Enrich_all_plot_select[which(Enrich_all_plot_select$freq>1),]$Description)))#2
range(Enrich_all_plot_select$Count)#3 11
range(Enrich_all_plot_select$Log_pvalue)#1.308557 2.735235 

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=group,y=Description,size=Count,colour=Log_pvalue))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,3,6,9,12),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,1,1.5,2,2.5,3),name='-log10(pvalue)')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))+
  theme_classic()+labs(x="",y="GO terms",title="AMA related kids DMR related genes: BP enrichment")+
  theme(legend.text = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 8,colour = 'black',vjust=1,hjust=1,angle = 90),
        legend.title = element_text(size = 12),
        legend.position ="right",legend.direction = "vertical")
plot_BP
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_DMR_related_DEGs_enrich_all_BP_selected_one.pdf",plot_BP,width = 10, height =10)
write.table(Enrich_all_plot_select, file="/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_DMRs_related_genes_GO_BP_merge_count_more_than_two_selected_one_wait.txt",row.names=T, col.names=T) 

#for selected  final 
Enrich_all_plot_select<-read.table("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_DMRs 2_related_genes_GO_BP_merge_count_more_than_two_selected_two.txt",sep="\t",header =T) 
head(Enrich_all_plot_select)
Enrich_all_plot_select$group<-factor(Enrich_all_plot_select$group,levels = c("Kid_hyper","Kid_hypo"))
Enrich_all_plot_select<-Enrich_all_plot_select[order(Enrich_all_plot_select$freq,Enrich_all_plot_select$group,decreasing = T ),]
Enrich_all_plot_select$Description<-factor(Enrich_all_plot_select$Description,levels = rev(c(unique(as.character(Enrich_all_plot_select$Description)))))
length(unique(as.character(Enrich_all_plot_select[which(Enrich_all_plot_select$freq>1),]$Description)))#5
range(Enrich_all_plot_select$Log_pvalue)#1.308557 2.735235 

plot_BP<-ggplot(Enrich_all_plot_select, aes(x=group,y=Description,size=Count,colour=Log_pvalue))+
  geom_point(alpha =0.5,na.rm = TRUE)+
  scale_size(breaks = c(1,3,6,9,12),range = c(1,6),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,1,1.5,2,2.5,3),name='-log10(pvalue)')+ 
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
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_DMRs 2_related_genes_GO_BP_merge_count_selected_final.pdf",plot_BP,width = 6, height =6)
write.table(Enrich_all_plot_select, file="/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/Kid_DMRs 2_related_genes_GO_BP_merge_count_selected_final.txt",row.names=T, col.names=T) 
