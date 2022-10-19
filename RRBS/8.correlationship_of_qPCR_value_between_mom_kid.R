rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(scales)
library(ggsci)
library(ggpubr)

#set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
show_col(ppCor)

data_analysis <-read.csv(file="/mnt/data/chenwei/huahua/manuscript/qPCR_value.csv")
head(data_analysis)
levels(data_analysis$gene_name)
# "CD24"       "HTRA3"   "IFIM1"   "MED12L"  "MS4A6A"  "NRPA2"   "SCD"     "SLC28A3" ##"DDX11"   "GMFR" 
data_analysis2 <- dcast(data_analysis,sample_name+gene_name ~ group , value.var = "Q_value")
head(data_analysis2)

#calculation the correlationship for kid and mother 
plot_cordata<-data_analysis2[,c("gene_name","kid","mom","sample_name")]
colnames(plot_cordata)<-c("Gene","kid","mom","Family")
plot_cordata$Age<-ifelse(grepl("A",plot_cordata$Family),"AMA","Young")

plot_cordata$Gene <-as.character(plot_cordata$Gene)
plot_cordata2<-plot_cordata[which(!(plot_cordata$Gene %in% c("DDX11","GMFR"))),]
#plot_cordata<-data_analysis_sig[which(data_analysis_sig$Gene %in% c("SCD","CD24")),]
ceiling(length(unique(plot_cordata2$Gene))/2)
cor_plot0<-ggplot(data=plot_cordata2, aes(x=kid, y=mom))+
  stat_smooth(method="lm",se=TRUE,colour="black",fill="cornflowerblue", size=1, alpha = 0.5)+
  geom_point(aes(colour = Age))+
  scale_color_manual(values=c(ppCor[c(1:2)]))+
  stat_cor(method = "spearman")+  
  facet_wrap(~ Gene, scales = "free",ncol = 3)+theme_bw()
cor_plot0
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure5/R_spearman_correlationship_QPCR_level_of_selected_AMA_common_DEGs_between_kids_and_mother_seprated.pdf",cor_plot0,width = 15, height =15)
##pearson
cor_plot2<-ggscatter(plot_cordata2, x = "kid", y = "mom",
          add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, cor.coef = TRUE)+ facet_wrap(~ Gene, scales = "free",ncol=3)
cor_plot2
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure5/R_pearson_correlationship_QPCR_level_of_selected_AMA_common_DEGs_between_kids_and_mother_seprated.pdf",cor_plot2,width = 15, height =15)

##R2 of pearson 
cor_plot4<-ggscatter(plot_cordata2, x = "kid", y = "mom", add = "reg.line",conf.int = TRUE,color = "Gene")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.5) + facet_wrap(~ Gene, scales = "free",ncol=3)+
  theme(legend.title = element_text(color = "black", size = 7),legend.position = "right")
cor_plot4
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure5/R2_pearson_correlationship_QPCR_level_of_selected_AMA_common_DEGs_between_kids_and_mother_seprated.pdf",cor_plot4,width = 15, height =15)

##R2 of spearman 
cor_plot5<-ggscatter(plot_cordata2, x = "kid", y = "mom", add = "reg.line",conf.int = TRUE,color = "Gene")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),   method = "spearman", label.x = 0.5) + facet_wrap(~ Gene, scales = "free",ncol=3)+
  theme(legend.title = element_text(color = "black", size = 7),legend.position = "right")
cor_plot5
ggsave(file="/mnt/data/chenwei/huahua/manuscript/figure/figure5/R2_spearman_correlationship_QPCR_level_of_selected_AMA_common_DEGs_between_kids_and_mother_seprated.pdf",cor_plot5,width = 15, height =15)
