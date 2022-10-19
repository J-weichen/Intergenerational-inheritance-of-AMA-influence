rm(list = ls())
#载入各种包
library(scales)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(grid)
library(ggsci)
library(data.table)
library(gridExtra)
library(cowplot)
#SET colors
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#step1:read  DMRs bed files
#readPeakFile读入peaks文件
#def 15%
#selection for kids DMRs

#FOR mother
#read meth data 
#read meth data 
compare_name<-"mother_AMA_vs_Young";split_region<- "200bp";depth<-"6X"
#data preparation
mother_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
head(mother_DMR_all_big)
myDiff15q005_all<-mother_DMR_all_big[,c("seqnames","start","end")]
myDiff15q005_hyper<-mother_DMR_all_big[which(mother_DMR_all_big$meth.diff> 0),c("seqnames","start","end")]
myDiff15q005_hypo<-mother_DMR_all_big[which(mother_DMR_all_big$meth.diff< 0),c("seqnames","start","end")]
colnames(myDiff15q005_all)<-colnames(myDiff15q005_hyper)<-colnames(myDiff15q005_hypo)<-c("chr","start","end")
dim(myDiff15q005_all);dim(myDiff15q005_hyper);dim(myDiff15q005_hypo)

#selected DMRs
target_DMRs<-myDiff15q005_all;DMR_type<-"All_DMR"
#target_DMRs<-myDiff15q005_hyper;DMR_type<-"hyper_DMR"
#target_DMRs<-myDiff15q005_hypo;DMR_type<-"hypo_DMR"

#target_DMRs2<-target_DMRs[,1:3]
head(target_DMRs);dim(target_DMRs)
#target_DMRs2$start<-as.numeric(as.character(target_DMRs2$start));target_DMRs2$end<-as.numeric(as.character(target_DMRs2$end))
table(target_DMRs$chr)
#all DMRs
#chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
#  43    31    30    20    18    18    16    22    31     8    23    34    14    12    30    17    16    26    24    27    22    34 

##hyper_DMR
##chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
##15    12    15    10     9     9    10    10    11     2    12    16     4     8    12     6     2    13    11    15    13    21
#hypo_DMR
#chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
#28    19    15    10     9     9     6    12    20     6    11    18    10     4    18    11    14    13    13    12     9    13

#step2: read region files
files_region_hg38<- list.files("/mnt/data/chenwei/qinmen_BR/00.ref_data/1-30regions/chenwei_file")
Element_percentage<-list()
cnane_region<-c()

for ( f_name in files_region_hg38[1:length(files_region_hg38)]){
  #f_name<-"MEs-23chr.bed"   ##test line 
  f_name1<- unlist(lapply(strsplit(f_name,"[.]"), function(x) x[1]))
  f_name2<- unlist(lapply(strsplit(f_name1,"-"), function(x) x[1]))
  print(as.character(f_name2))
  
  regionnames<-c(f_name2,paste0("no_",f_name2))
  Anno_region <- as.data.frame(read.table(paste("/mnt/data/chenwei/qinmen_BR/00.ref_data/1-30regions/chenwei_file/",f_name,sep=""),header = F,sep = "\t"))
  Anno_region<-Anno_region[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",Anno_region$V1,invert=TRUE),]
  dim(Anno_region);head(Anno_region)
  #target_region<-Anno_region[,1:4]
  #colnames(target_region)<-c("chr","start","end","region_name")
  target_region<-Anno_region[,1:3]
  colnames(target_region)<-c("chr","start","end")
  #calculate intersection
  x = data.table(target_DMRs)
  y = data.table(target_region)
  x;y;dim(x);dim(y)
  ####key work 
  setkey(y, chr, start, end)
  all_data<-foverlaps(x, y, type="any")
  inter_data<-foverlaps(x, y, type="any", nomatch=NULL)
  all_data$DMR_name<-paste0(all_data$chr,"_",all_data$i.start,"_",all_data$i.end)
  inter_data$DMR_name<-paste0(inter_data$chr,"_",inter_data$i.start,"_",inter_data$i.end)
  numlist<-c(length(unique(inter_data$DMR_name)),length(unique(all_data$DMR_name))-length(unique(inter_data$DMR_name)))
  numlist1<-numlist
  names(numlist)<-regionnames
  data.frame(numlist)
  ##data for plot pie 
  df <- data.frame(value = numlist1,Group = regionnames,Group2 = c("Overlap","No_overlap")) %>%
    mutate(Group = factor(Group, levels = rev(regionnames)),
           Group2=  factor(Group2, levels = c("No_overlap","Overlap")),
           #cumulative = length(unique(all_data$DMR_name)),
           # midpoint = cumulative - value / 2,
           label1 =paste0(round(value / sum(value) * 100, 1),"%"),
           label2 = paste0(Group,":"," (",round(value / sum(value) * 100, 1),"%)"),
           label3 = paste0(Group,":"," (n=",value,")"),
           label4 = paste0(Group2,":"," (n=",value,")"))
  df$label3<-factor(df$label3, levels = rev(df$label3))
  df$label4<-factor(df$label4, levels = rev(df$label4))
  df$fraction = df$value / sum(df$value)
  df$ymax = cumsum(df$fraction)
  df$ymin = c(0, head(df$ymax, n = -1))
  
  Element_percentage<-c(Element_percentage,list(df))
  cnane_region<-c(cnane_region,f_name2)
  print(df)
}

length(Element_percentage);length(cnane_region)
names(Element_percentage)<-cnane_region
head(Element_percentage$LTR)
head(Element_percentage[["LTR"]])

saveRDS(Element_percentage, file = paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/location_annotation_Element_percentage_",compare_name,"_",DMR_type,".rds"))

##step3: plot for percentage
#Element_percentage<-readRDS(file = paste0("D:/3.秦萌高龄甲基化/Arrest_project/DMR_file/region_anotation/Element_percentage_",DMR_type,".rds"))
names(Element_percentage)
target_region_list<-c("LTR","LINE","SINE","CpGislands","ICRs_67","MEs","HCP","ICP","LCP","Low_complexity","Retroposon_SVA","Transposon","Satellite","Micro_satellites")
length(target_region_list)

#target_region_list<-c("Low_complexity","LTR","LINE","SINE","Retroposon_SVA","Transposon",
#                      "Satellite","Micro_satellites","CpGislands","ICRs_380","ICRs_67","MEs",
#                      "Centromeres","Intragenic_region_KSM","Intergenic_region_KSM",
#                      "HCP","ICP","LCP","UTR5","UTR3","Exons","Introns","All_Genes","Intergenic_1kb")
p_pie_list<-list()
p_ring_list<-list()
cnane_region2<-c()

for ( region_name in target_region_list){
  #  region_name<-"LTR"  ##test line 
  print(as.character(region_name))
  plot_percent<-Element_percentage[[which(names(Element_percentage) == region_name)]]
  
  p_pie<-ggplot(plot_percent, aes(x = 1, weight = value, fill =  label4)) +
    geom_bar(width = 1, position = "stack") +
    coord_polar(theta = "y") +
    scale_fill_manual(values=ppCor[c(9,10)])+
    theme_bw() + 
    #geom_text(aes(x = 1.3, y = midpoint, label = label1))+
    geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), x = 1, label = label1), size = 5)+   ## 在图中加上百分比：x 调节标签到圆心的距离, y 调节标签的左右位置  
    labs(x = "", y = "", title = "") +
    theme(axis.ticks = element_blank()) +  
    theme(legend.title = element_text(color = "black", size = 7),legend.position = "right")+
    theme(axis.text.x =element_blank(),axis.text.y=element_blank())  +
    theme(panel.grid=element_blank()) + 
    theme(panel.border=element_blank())+guides(fill = guide_legend(title = region_name ))
  # ggtitle("Satellite")+  theme(plot.title =element_text(hjust = 0.5,vjust = 0.5))
  
  ##plot rings
  p_ring<- ggplot(data = plot_percent, aes(fill = label4, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
    geom_rect(colour = "grey30", show.legend = T) +
    coord_polar(theta = "y") +
    scale_fill_manual(values=ppCor[c(9,10)])+
    labs(x = "", y = "", title = "") + 
    xlim(c(0, 4)) +
    theme_bw() +
    theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.border=element_blank()) + 
    theme(legend.title = element_text(color = "black", size = 7), legend.position = "right")+
    geom_text(aes(x = 3.5, y = ((ymin+ymax)/2), label = label1)) +guides(fill = guide_legend(title = region_name))
  #ggtitle("Satellite")+  theme(plot.title =element_text(hjust = 0.5,vjust = 0.5))
  # ggsave(p_pie,file=paste0("/media/data2/lucunlin/result/","In_Vivo-In_vitro_Pie_200bp_def50p05_DMRs_in_",region_name,"_pie_",DMR_type,".pdf"))   ##ggsave批量生成name变量为名字的pdf文件
  #ggsave(p_pie,file=paste0("/media/data2/lucunlin/result/","In_Vivo-In_vitro_Pie_200bp_def50p05_DMRs_in_",region_name,"_pie_",DMR_type,".pdf"))   ##ggsave批量生成name变量为名字的pdf文件
  
  pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/location_annotation_AMA_Pie_200bp_def15q005_DMRs_in_",region_name,"_",compare_name,"_",DMR_type,".pdf"))
  Pie_ring_plot<- ggdraw()+ draw_plot(p_pie, x=0, y=0, width=0.5, height = 1)+draw_plot(p_ring, x=0.5, y=0, width = 0.5, height = 1)
  p<-annotate_figure(Pie_ring_plot, top=text_grob(paste0("AMA_Pie_200bp_def15p05_DMRs_in_",region_name), color = "black",face = "bold", size=12))
  print(p)
  dev.off()
  p_pie_list<-c(p_pie_list,list(p_pie))
  p_ring_list<-c(p_ring_list,list(p_ring))
  cnane_region2<-c(cnane_region2,region_name)
}

length(p_pie_list);length(p_ring_list);length(cnane_region2)
names(p_pie_list)<-cnane_region2
names(p_ring_list)<-cnane_region2

pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/location_annotation_AMA_Pie_rings_200bp_def15q005_DMRs_in_",compare_name,"_",DMR_type,".pdf"))

for (region_name in cnane_region2){
  #region_name<-"Intergenic_1kb"  ##test line 
  
  print(as.character(region_name))
  Pie_ring_plot<- ggdraw()+ draw_plot(p_pie_list[[region_name]], x=0, y=0, width=0.5, height = 1)+draw_plot(p_ring_list[[region_name]], x=0.5, y=0, width = 0.5, height = 1)
  p<-annotate_figure(Pie_ring_plot, top=text_grob(paste0("Arrest_abortion_Pie_200bp_def15q005_DMRs_in_",region_name,"_",DMR_type), color = "black",face = "bold", size=12))
  print(p)
}

dev.off()

#https://www.cnblogs.com/triple-y/p/11635154.html
#https://www.cnblogs.com/triple-y/p/11635154.html
p_pie_group <- ggpubr::ggarrange(p_pie_list[[1]],p_pie_list[[2]],p_pie_list[[3]],p_pie_list[[4]],p_pie_list[[5]],
                                 p_pie_list[[6]],p_pie_list[[7]],p_pie_list[[8]],p_pie_list[[9]],p_pie_list[[10]],
                                 p_pie_list[[11]],p_pie_list[[12]],p_pie_list[[13]],p_pie_list[[14]],
                                 #labels = c('A', 'B', 'C', 'D'),  font.label = list(color = 'black'),
                                 nrow = 3, ncol = 5)
annotate_figure(p_pie_group, top=text_grob(paste0("AMA_6X_CpG_200bp_DMRs_myDiff15p005_",DMR_type), color = "black",face = "bold", size=12))
#In_Vivo-In_vitro_Df50q5_200bp_DMR_percent_pie_merge_plot

p_ring_group <- ggpubr::ggarrange(p_ring_list[[1]],p_ring_list[[2]],p_ring_list[[3]],p_ring_list[[4]], p_ring_list[[5]],
                                  p_ring_list[[6]],p_ring_list[[7]],p_ring_list[[8]],p_ring_list[[9]],p_ring_list[[10]],
                                  p_ring_list[[11]],p_ring_list[[12]],p_ring_list[[13]],p_ring_list[[14]],
                                  #labels = c('A', 'B', 'C', 'D'),  font.label = list(color = 'black'),
                                  nrow = 3, ncol = 5)
annotate_figure(p_ring_group, top=text_grob(paste0("AMA_6X_CpG_200bp_DMRs_myDiff15p005_",DMR_type), color = "black",face = "bold", size=12))

#pdf(paste0("D:/huahua/1.DMR_called/DMR_local_annotation/","all_Arrest_15_Pie_ring_merge_6X_CpG_200bp_def15p005_DMRs_in_",DMR_type,".pdf"))
p1<-annotate_figure(p_pie_group,  top=text_grob(paste0("AMA_6X_CpG_200bp_DMRs_myDiff15p005_",DMR_type), color = "black",face = "bold", size=12))
p2<-annotate_figure(p_ring_group, top=text_grob(paste0("AMA_6X_CpG_200bp_DMRs_myDiff15p005_",DMR_type), color = "black",face = "bold", size=12))

#dev.off()
ggsave(file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/location_annotation_AMA_14_Pie_ring_merge_6X_CpG_200bp_def15p005_DMRs_in_",compare_name,"_",DMR_type,"_1.pdf"),p1,width = 600, height =360, units = "mm")
ggsave(file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/location_annotation_AMA_14_Pie_ring_merge_6X_CpG_200bp_def15p005_DMRs_in_",compare_name,"_",DMR_type,"_2.pdf"),p2,width = 600, height =360, units = "mm")
