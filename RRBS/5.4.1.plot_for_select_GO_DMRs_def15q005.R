rm(list = ls())
#载入各种包
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

##plot for selected GO
GO_simple_plot_hyper <-read.csv("/media/data2/lucunlin/enrich_plot/anotation_and_enrich_all_select/In_Vivo-In_vitro_hyper_work.csv",header =T,stringsAsFactors = FALSE,sep = ",")
GO_simple_plot_hypo <-read.csv("/media/data2/lucunlin/enrich_plot/anotation_and_enrich_all_select/In_Vivo-In_vitro_hypo_work.csv",header =T,stringsAsFactors = FALSE,sep = ",")
head(GO_simple_plot_hyper);head(GO_simple_plot_hypo)
str(GO_simple_plot_hyper);str(GO_simple_plot_hypo)


GO_simple_plot<-as.data.frame(rbind(GO_simple_plot_hyper,GO_simple_plot_hypo))
GO_simple_plot$group<-c(rep("hyper",nrow(GO_simple_plot_hyper)),rep("hypo",nrow(GO_simple_plot_hypo)))
str(GO_simple_plot)
GO_simple_plot$qvalue_log10<-c(-log(GO_simple_plot$qvalue,10))
GO_simple_plot$pvalue_log10<-c(-log(GO_simple_plot$pvalue,10))

#BP_plot$group<- factor(x =BP_plot$group,levels =c("Trobt_Vims","EVTs","VCTs","STBs") )
GO_simple_plot$group<- factor(x =GO_simple_plot$group,levels =c("hyper","hypo") )
levels(GO_simple_plot$group)

GO_simple_plot1<-GO_simple_plot[,c("Description","group","qvalue_log10","Count")]
GO_simple_plot1$Description<-as.character(GO_simple_plot1$Description)
GO_simple_plot1<-GO_simple_plot1[order(GO_simple_plot1$group,GO_simple_plot1$qvalue_log10),]
str(GO_simple_plot1)
head(GO_simple_plot1)

PP_plot<-GO_simple_plot1
head(PP_plot)
#pdf("In_Vivo-In_vitro-DMRs_200bp_50p5_hypo-hyper_GO_simple_all_select.pdf")

p <- ggplot(PP_plot,aes(qvalue_log10 ,reorder(Description,qvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(p.adjust)",y="Description",title="GO_simple enrichment")+
  scale_size_continuous(range=c(5,14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=55))+
  geom_text(aes(label=sprintf("%.3f",qvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")

p


GO_simple_plot2<-GO_simple_plot[,c("Description","group","pvalue_log10","Count")]
GO_simple_plot2$Description<-as.character(GO_simple_plot2$Description)
GO_simple_plot2<-GO_simple_plot2[order(GO_simple_plot2$group,GO_simple_plot2$pvalue_log10),]
str(GO_simple_plot2)
head(GO_simple_plot2)

PP_plot<-GO_simple_plot2
head(PP_plot)
ggplot(PP_plot,aes(pvalue_log10 ,reorder(Description,pvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(pvalue)",y="Description",title="GO enrichment")+
  scale_x_continuous(breaks=seq(0, 14, 1)) +
  scale_size_continuous(range=c(1,9))+
  #scale_y_discrete(labels=function(x) str_wrap(x, width=40))+
  #geom_text(aes(label=sprintf("%.3f",pvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=14),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10),legend.position = "right")
