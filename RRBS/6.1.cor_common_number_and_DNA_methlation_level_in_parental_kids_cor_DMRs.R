rm(list = ls())
library(ggpubr)
library(scales)
library(ggsci)
library(stringr)
library(UpSetR)


pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)
#甲基化对应基因列表的读取
#data preparation
#15% df
compare_name1<-"kids_AMA_vs_Young";split_region1<- "200bp"
#data preparation
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region1,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
head(data_region_merge_DMR_all_big)

DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
DMR_hyper1<-unique(as.character(DMR_hyper$Row.names))
DMR_hypo1<-unique(as.character(DMR_hypo$Row.names))
length(DMR_hyper1);length(DMR_hypo1)#182 #269

compare_name2<-"mother_AMA_vs_Young";split_region2<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name2,"_",split_region2,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
DMR_hyper3<-unique(as.character(DMR_hyper$Row.names))
DMR_hypo3<-unique(as.character(DMR_hypo$Row.names))
length(DMR_hyper3);length(DMR_hypo3)# 236 #280

compare_name2<-"father_AMA_vs_Young";split_region2<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name2,"_",split_region2,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
DMR_hyper4<-unique(as.character(DMR_hyper$Row.names))
DMR_hypo4<-unique(as.character(DMR_hypo$Row.names))
length(DMR_hyper4);length(DMR_hypo4)#206 #258

listinput<-c(list(DMR_hyper1),list(DMR_hyper3),list(DMR_hyper4),list(DMR_hypo1),list(DMR_hypo3),list(DMR_hypo4))
names(listinput)<-c("Kids_hyper","Mother_hyper","Father_hyper","Kids_hypo","Mother_hypo","Father_hypo" )

common_hyper<-Reduce(intersect,listinput[c("Kids_hyper","Mother_hyper","Father_hyper")])
common_hypo<-Reduce(intersect,listinput[c("Kids_hypo","Mother_hypo","Father_hypo")])

commonlist_kid_mom_hyper<-Reduce(intersect,listinput[c("Kids_hyper","Mother_hyper")])
commonlist_kid_mom_hypo<-Reduce(intersect,listinput[c("Kids_hypo","Mother_hypo")])
mother_SF_DMRs_hyper<- commonlist_kid_mom_hyper[!(commonlist_kid_mom_hyper %in% common_hyper)]
mother_SF_DMRs_hypo<- commonlist_kid_mom_hypo[!(commonlist_kid_mom_hypo %in% common_hypo)]

commonlist_kid_dad_hyper<-Reduce(intersect,listinput[c("Kids_hyper","Father_hyper")])
commonlist_kid_dad_hypo<-Reduce(intersect,listinput[c("Kids_hypo","Father_hypo")])
father_SF_DMRs_hyper<- commonlist_kid_dad_hyper[!(commonlist_kid_dad_hyper %in% common_hyper)]

common_mix<-Reduce(intersect,listinput[c("Kids_hypo","Mother_hyper","Father_hypo")])
father_SF_DMRs_hypo<- commonlist_kid_dad_hypo[!(commonlist_kid_dad_hypo %in% c(common_hypo,common_mix))]

length(common_hyper);length(common_hypo);length(mother_SF_DMRs_hyper);length(mother_SF_DMRs_hypo);length(father_SF_DMRs_hyper);length(father_SF_DMRs_hypo)
#4 9 31 34 22 36

#read correlated DMRs
#for mother 
mother_cordata <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlated_pvalue_for_kid_AMA_DMR_between_kids_and_mother.txt",header=T,sep="\t")
MF_cor_DMR_sig<-as.character(mother_cordata[which(mother_cordata$group == "candidate_region"),]$DMR_region)
father_cordata <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlated_pvalue_for_kid_AMA_DMR_between_kids_and_father.txt",header=T,sep="\t")
FF_cor_DMR_sig<-as.character(father_cordata[which(father_cordata$group == "candidate_region"),]$DMR_region)

length(MF_cor_DMR_sig);length(FF_cor_DMR_sig)#66 58

##call cor ratio in commmon DMRs
kid_mom_common_DMRs<-c(commonlist_kid_mom_hyper,commonlist_kid_mom_hypo)
MF_cor_DMR_freq<-data.frame(table(kid_mom_common_DMRs %in% MF_cor_DMR_sig))
colnames(MF_cor_DMR_freq)<-c("state","MF_freq")
kid_dad_common_DMRs<-c(commonlist_kid_dad_hyper,commonlist_kid_dad_hypo)
DF_cor_DMR_freq<-data.frame(table(kid_dad_common_DMRs %in% FF_cor_DMR_sig))
colnames(DF_cor_DMR_freq)<-c("state","DF_freq")
common_corDMR_freq<-merge(MF_cor_DMR_freq,DF_cor_DMR_freq)
common_corDMR_freq

trans_MF_cor_common <-read.table(file="/mnt/data/chenwei/huahua/5.valificat_result/3.correlated_DEGs/correlated_pvalue_of_AMA_common_DEG_between_kids_and_mother.txt",header=T,sep="\t")
table(trans_MF_cor_common$group)
MF_cor_DEG_freq<-data.frame(table(trans_MF_cor_common$group))
colnames(MF_cor_DEG_freq)<-c("state","MF_trans_freq")
MF_cor_DEG_freq$state<-c("TRUE","FALSE")
common_three_freq<-merge(common_corDMR_freq,MF_cor_DEG_freq)
plot_data <- melt(common_three_freq,variable.name="group",value.name = "number",id.vars = c("state"))
plot_data$state<-factor(plot_data$state,levels = c("TRUE","FALSE"))
plot_data$group<-factor(plot_data$group,levels = c("MF_trans_freq","MF_freq","DF_freq"))

ggplot(plot_data,aes(group,number,fill=state))+geom_bar(stat="identity",position="fill")

plot_num_decrease<-ggplot(data=plot_data, mapping=aes(x= group,y=number,fill=state))+
  geom_bar(stat="identity",width=0.8,position= 'stack')+scale_fill_manual(name="Types",values=ppCor[c(3,6)])+
  theme_classic()+labs(x="group",y="DMRs/DEGs number",title="")+
  scale_y_continuous(breaks = seq(0,100,by=20))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),axis.text.y  = element_text(size = 8,colour = 'black'),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_num_decrease
plot_num_decrease2<-plot_num_decrease+geom_text(aes(label=number),position=position_stack(vjust=0.5))#3 vjust参数用来调整标签的为重，vjust=0.5将标签放到对应部位的中部
plot_num_decrease2
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/barplot_number_of_cor_parent_kids_DMRs_DEGs_in_commom_region.pdf",plot_num_decrease2,width = 6, height =6)

#show genes for cor_DMRs
MF_cor_common_DMR<-kid_mom_common_DMRs[which(kid_mom_common_DMRs %in% MF_cor_DMR_sig)]
DF_cor_common_DMR<-kid_dad_common_DMRs[which(kid_dad_common_DMRs %in% FF_cor_DMR_sig)]
DMRs_unique<-unique(c(MF_cor_common_DMR,DF_cor_common_DMR))

#data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region"),]
data_region_merge_DMR_all_big<-read.table("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/kids_AMA_vs_Young_200bp_DMR_myDiff15q005_merge_data.txt",header=T,sep="\t")
cor_DMR_gene<-data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$Row.names %in% DMRs_unique),c("bin_region","distanceToTSS","annotation","SYMBOL")]
colnames(cor_DMR_gene)<-c("DMR_region","distanceToTSS","annotation","SYMBOL")
write.table(cor_DMR_gene, file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/all_cor_DMR_gene_parent_kids_common_AMA_vs_Young_200bp_DMR_myDiff15q005.txt"),quote=F,row.names=F,col.names=T,sep = "\t")
#for mother 
data_analysis <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_mother.txt",header=T,sep="\t")
data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region"),]
data_analysis_sig$DMR_region <-as.character(data_analysis_sig$DMR_region)
plot_cordata<-merge(cor_DMR_gene,data_analysis_sig)
plot_cordata$region_gene<-paste0(plot_cordata$DMR_region,"_",plot_cordata$SYMBOL)
length(unique(plot_cordata$DMR_region));ceiling(length(unique(plot_cordata$DMR_region))/5)
cor_plot0<-ggplot(data=plot_cordata, aes(x=kids, y=parent))+geom_point()+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=plot_cordata, method = "spearman")+  
  facet_wrap(~ region_gene, scales = "free",ncol = 5)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/common_DMRs_correlationship_methlation_level_of__AMA_DMR_between_kids_and_mother_seprated.pdf",cor_plot0,width = 60*5, height =60*ceiling(length(unique(plot_cordata$DMR_region))/5), units = "mm")
#for father 
data_analysis <-read.table(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/correlationship_methlation_level_of_AMA_DMR_between_kids_and_father.txt",header=T,sep="\t")
data_analysis_sig<-data_analysis[which(data_analysis$group == "candidate_region"),]
data_analysis_sig$DMR_region <-as.character(data_analysis_sig$DMR_region)
plot_cordata<-merge(cor_DMR_gene,data_analysis_sig)
plot_cordata$region_gene<-paste0(plot_cordata$DMR_region,"_",plot_cordata$SYMBOL)
length(unique(plot_cordata$DMR_region));ceiling(length(unique(plot_cordata$DMR_region))/5)
cor_plot1<-ggplot(data=plot_cordata, aes(x=kids, y=parent))+geom_point()+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=plot_cordata, method = "spearman")+  
  facet_wrap(~ region_gene, scales = "free",ncol = 5)
ggsave(file="/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/common_DMRs_correlationship_methlation_level_of__AMA_DMR_between_kids_and_mother_seprated.pdf",cor_plot1,width = 60*5, height =60*ceiling(length(unique(plot_cordata$DMR_region))/5), units = "mm")
