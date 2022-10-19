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


#upset plot for all eight lists of DMRs
listinput_new1<-c(list(DMR_hyper1),list(DMR_hyper3),list(DMR_hyper4),list(DMR_hypo1),list(DMR_hypo3),list(DMR_hypo4),list(MF_cor_DMR_sig),list(FF_cor_DMR_sig))
names(listinput_new1)<-c("Kids_hyper","Mother_hyper","Father_hyper","Kids_hypo","Mother_hypo","Father_hypo","MF_cor_DMR_sig","FF_cor_DMR_sig")
upset_dataframe<-as.data.frame(fromList(listinput_new1))
dim(upset_dataframe)# 1235     8
upset_dataframe[1:10,1:5]
upset_dataframe_rowSums<-rowSums(upset_dataframe)
range(upset_dataframe_rowSums)#1:5
upset_dataframe_colSums<-colSums(upset_dataframe)
range(upset_dataframe_colSums)#58 280
upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]
upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400)

list2_q1<-list(query = intersects, params = list("mother_SF_DMRs_hyper","MF_cor_DMR_sig"), color = "red", active = T,query.name = "MS_DMRs")
list2_q2<-list(query = intersects, params = list("mother_SF_DMRs_hypo","MF_cor_DMR_sig"), color = "red", active = T,query.name = "MS_DMRs")
list2_q3<-list(query = intersects, params = list("father_SF_DMRs_hyper","FF_cor_DMR_sig"), color = "purple", active = T,query.name = "FS_DMRs")
list2_q4<-list(query = intersects, params = list("father_SF_DMRs_hypo","FF_cor_DMR_sig"), color = "purple", active = T,query.name = "FS_DMRs")

#list2_q2<-list(query = intersects, params = list("mother_SF_DMRs_hypo","MF_cor_DMR_sig"), color = "navy", active = T,query.name = "common hypo DMRs")
#list2_q3<-list(query = intersects, params = list("Kids_hyper", "Mother_hyper"), color = "orange", active = T,query.name = "mother_SF_DMRs")
#list2_q4<-list(query = intersects, params = list("Kids_hypo","Mother_hypo"), color = "orange", active = T,query.name = "mother_SF_DMRs")
#list2_q6<-list(query = intersects, params = list("Kids_hypo","Father_hypo"), color = "purple", active = T,query.name = "father_SF_DMRs")
pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/upsetplot_cor_DMRs_and_six_lists_DMRs_AMA_vs_Young_200bin_q005_dif15.pdf",width = 10,height=8)

upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400,
      #  sets=c("FMK_common_hyper","mother_SF_DMRs_hyper","father_SF_DMRs_hyper",
      #         "FMK_common_hypo","mother_SF_DMRs_hypo","father_SF_DMRs_hypo","MF_cor_DMR_sig","FF_cor_DMR_sig"),
      sets=c("Kids_hyper","Kids_hypo","Mother_hyper","Mother_hypo","Father_hyper","Father_hypo","MF_cor_DMR_sig","FF_cor_DMR_sig"),
      keep.order = TRUE,
      query.legend = "top",
      sets.bar.color = "brown",
      #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
      shade.color="pink",
      #      matrix.color="purple",
      order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
      point.size = 3,line.size = 1.3,
      mainbar.y.label = "DMRs number Intersections", sets.x.label = "DMRs number per subset",
      mb.ratio = c(0.60, 0.40),
      text.scale = c(1.5, 1.5,1.2,1.5,1.5,1),
      show.numbers = 'yes'#,
      #queries = list(list2_q1,list2_q2,list2_q3,list2_q4)
)
dev.off()

#for hyper kids DMRs
listinput_new1<-c(list(DMR_hyper1),list(DMR_hyper3),list(DMR_hyper4),list(DMR_hypo1),list(DMR_hypo3),list(DMR_hypo4),list(MF_cor_DMR_sig),list(FF_cor_DMR_sig))
names(listinput_new1)<-c("Kids_hyper","Mother_hyper","Father_hyper","Kids_hypo","Mother_hypo","Father_hypo","MF_cor_DMR_sig","FF_cor_DMR_sig")

venn_Age_DMR1<-venn.diagram(list(Kids_hyper=DMR_hyper1,MF_cor_DMR_sig=MF_cor_DMR_sig,FF_cor_DMR_sig=FF_cor_DMR_sig,
                                        Mother_hyper=DMR_hyper3,Father_hyper=DMR_hyper4),
                                   alpha=c(0.8,0.8,0.8,0.8,0.8),lwd=1,lty=1,  col="black" , fill=ppCor[c(1,3,5:7)], 
                                   cex = 1.5, cat.col=ppCor[c(1,3,5:7)],cat.fontface=4,cat.cex = 1.5,
                                   main.cex = 2, main.fontface = 2, main.fontfamily = 3,  filename = NULL)
venn_Age_DMR2<-venn.diagram(list(Kids_hypo=DMR_hypo1,MF_cor_DMR_sig=MF_cor_DMR_sig,FF_cor_DMR_sig=FF_cor_DMR_sig,
                                        Mother_hypo=DMR_hypo3,Father_hypo=DMR_hypo4),
                                   alpha=c(0.8,0.8,0.8,0.8,0.8),lwd=1,lty=1,  col="black" , fill=ppCor[c(2,3,5:7)], 
                                   cex = 1.5, cat.col=ppCor[c(2,3,5:7)],cat.fontface=4,cat.cex = 1.5,
                                   main.cex = 2, main.fontface = 2, main.fontfamily = 3,  filename = NULL)
pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/venn_five_lists_DMRs_AMA_vs_Young_200bin_q005_dif15.pdf",width = 8,height=8)
grid.newpage(); 
grid.draw(venn_Age_DMR1)
grid.newpage();
grid.draw(venn_Age_DMR2)
dev.off()

#upset plot for all eight lists of DMRs for parent kids common
listinput_new<-c(list(common_hyper),list(common_hypo),list(mother_SF_DMRs_hyper),list(mother_SF_DMRs_hypo),
                 list(father_SF_DMRs_hyper),list(father_SF_DMRs_hypo),list(MF_cor_DMR_sig),list(FF_cor_DMR_sig))
names(listinput_new)<-c("FMK_common_hyper","FMK_common_hypo","mother_SF_DMRs_hyper","mother_SF_DMRs_hypo",
                        "father_SF_DMRs_hyper","father_SF_DMRs_hypo","MF_cor_DMR_sig","FF_cor_DMR_sig")

upset_dataframe<-as.data.frame(fromList(listinput_new))
dim(upset_dataframe)#213   8
upset_dataframe[1:10,1:5]
upset_dataframe_rowSums<-rowSums(upset_dataframe)
range(upset_dataframe_rowSums)#1:3
upset_dataframe_colSums<-colSums(upset_dataframe)
range(upset_dataframe_colSums)#2 54
upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]
upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400)

list2_q1<-list(query = intersects, params = list("mother_SF_DMRs_hyper","MF_cor_DMR_sig"), color = "red", active = T,query.name = "MS_DMRs")
list2_q2<-list(query = intersects, params = list("mother_SF_DMRs_hypo","MF_cor_DMR_sig"), color = "red", active = T,query.name = "MS_DMRs")
list2_q3<-list(query = intersects, params = list("father_SF_DMRs_hyper","FF_cor_DMR_sig"), color = "purple", active = T,query.name = "FS_DMRs")
list2_q4<-list(query = intersects, params = list("father_SF_DMRs_hypo","FF_cor_DMR_sig"), color = "purple", active = T,query.name = "FS_DMRs")

#list2_q2<-list(query = intersects, params = list("mother_SF_DMRs_hypo","MF_cor_DMR_sig"), color = "navy", active = T,query.name = "common hypo DMRs")
#list2_q3<-list(query = intersects, params = list("Kids_hyper", "Mother_hyper"), color = "orange", active = T,query.name = "mother_SF_DMRs")
#list2_q4<-list(query = intersects, params = list("Kids_hypo","Mother_hypo"), color = "orange", active = T,query.name = "mother_SF_DMRs")
#list2_q6<-list(query = intersects, params = list("Kids_hypo","Father_hypo"), color = "purple", active = T,query.name = "father_SF_DMRs")
pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/upsetplot_cor_DMRs_and_fine_common_parent_kids_DMRs_AMA_vs_Young_200bin_q005_dif15.pdf",width = 8,height=6)

upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400,
      sets=c("FMK_common_hyper","mother_SF_DMRs_hyper","father_SF_DMRs_hyper",
             "FMK_common_hypo","mother_SF_DMRs_hypo","father_SF_DMRs_hypo", "MF_cor_DMR_sig","FF_cor_DMR_sig"),
      keep.order = TRUE,
      query.legend = "top",
      sets.bar.color = "brown",
      #      sets.bar.color=c(rep("black",times=5),ppCor[9:11]),
      shade.color="pink",
      #      matrix.color="purple",
      order.by = c("freq", "degree"),decreasing = c(TRUE,FALSE),
      point.size = 3,line.size = 1.3,
      mainbar.y.label = "DMRs number Intersections", sets.x.label = "DMRs number per subset",
      mb.ratio = c(0.60, 0.40),
      text.scale = c(1.5, 1.5,1.2,1.5,1.5,1),
      show.numbers = 'yes',
      queries = list(list2_q1,list2_q2,list2_q3,list2_q4)
)

dev.off()

###ratio for DMRs with qvalue < 0.05 in parents samples
#for father
father_data_merge_meth_all<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/father_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv",row.names = 1,check.names=F)
head(father_data_merge_meth_all)
father_data_merge_meth_all$bin_region<-unlist(lapply(strsplit(rownames(father_data_merge_meth_all),"[.]"), function(x) paste(x[1],x[2],x[3],sep = "_")))
father_DMRs_meth<-father_data_merge_meth_all[which(father_data_merge_meth_all$bin_region %in% FF_cor_DMR_sig),]

head(father_DMRs_meth,n=3)
dim(father_DMRs_meth)#58 27
table(father_data_merge_meth_all$qvalue <0.05)
#FALSE   TRUE 
#490197 211740 
table(father_DMRs_meth$qvalue <0.05)
#FALSE   TRUE 
# 6    52 

#卡方比较
data_kf1<-as.data.frame(table(father_data_merge_meth_all$qvalue <0.05))
data_kf2<-as.data.frame(table(father_DMRs_meth$qvalue <0.05))
data_kf<-matrix(c(rev(data_kf1$Freq),rev(data_kf2$Freq)),nrow = 2,ncol=2, byrow = T)
como <- chisq.test(data_kf)
p_HG<- ifelse(min(como$expected)>=5, round(chisq.test(data_kf,correct = F)$p.value,20),round(fisher.test(data_kf)$p.value,20))
p_HG#0
data_kf3<-as.data.frame(data_kf)
colnames(data_kf3)<-c("sig","no_sig");rownames(data_kf3)<-c(paste0("All[n=",nrow(father_data_merge_meth_all),"]"),paste0("Common[n=",nrow(father_DMRs_meth),"]"))
data_kf3$type<-rownames(data_kf3)
#cell_pvalue<-data.frame(Var1=cell_name,p_value=p_value)
#cell_number_merge3<-merge(cell_number_merge2,cell_pvalue,by="Var1")
cell_number_merge4 <- melt(data_kf3,variable.name="type2",value.name = "Num",id.vars = c("type") )
cell_number_merge4$p_value<-p_HG

Pos_rate_plot<-ggplot(data=cell_number_merge4, mapping=aes(x= type,y=Num,fill=type2))+
  geom_bar(stat="identity",width=0.8,position='fill')+
  geom_text(aes(x= 1.5, y=1.1, label = round(p_value,20)), data =cell_number_merge4,vjust = 1,color="black", size=5 )+
  scale_fill_manual(name="Types",values=ppCor[c(7,5)])+
  theme_classic()+labs(x="",y="AMA significant rate",title=paste0("significant rate for father"))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=1,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
ggsave("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/Group_compasion_for_significant_bin_rate_for_father_all_and_FF_cor_DMR_sig.pdf",Pos_rate_plot,width=4, height=3)


#差值分布
father_DMRs_meth$group1<-as.factor(ifelse(father_DMRs_meth$meth.diff>0, "hyper", "hypo"))
head(father_DMRs_meth)
plot_distribute<-ggplot(father_DMRs_meth, aes(x=meth.diff)) +
  geom_histogram(aes(y=..density..,fill=group1), binwidth=1,colour="gray") + # 这一步很重要,使用density代替y轴
  scale_fill_manual(name="group1",values=ppCor[1:2])+
  geom_density(alpha=.2, fill="#FF6666")+ # 重叠部分采用透明设置
  geom_vline(xintercept=c(-50,-25,-15,15,25,50),linetype="dashed",colour=c("blue","blue","blue","red","red","red"))+
  theme_bw()+scale_x_continuous(breaks=seq(-55, 55, 5))+
  theme(axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5,  hjust = 0.5))+
  theme(axis.title.y = element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5))+labs(title = "methykit difference: AMA vs Young")

ggsave("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/father_meth_defference_distribution_for_FF_cor_DMR_sig.pdf",plot_distribute,width=8, height=6)

#plot heatmap
father_DMRs_meth<-father_DMRs_meth[order(father_DMRs_meth$meth.diff,decreasing = T),]
father_sample_list<-c("10-2-1D","16-2-4C","17-2-2E","18-2","19-2-2H","20F-4G","6-2","7F","8-2","E16F",
                      "11-2","12-2-3H","13-2-1H","14-2","15-2","1F","2-2","3F","4-2","5-2")
father_heat_plot<-father_DMRs_meth[,father_sample_list]
head(father_heat_plot)
huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
rownames(huahua_meta)<-huahua_meta$library_code
huahua_meta2<-huahua_meta[father_sample_list,]
head(huahua_meta2)
colnames(father_heat_plot)<-huahua_meta2$analysis_name
as.character(huahua_meta2$analysis_name[order(huahua_meta2$analysis_name,decreasing = F)])
sample_list2 <-c(paste0("AMA_F_",1:10),paste0("YOUNG_F_",1:10))
father_heat_plot<-father_heat_plot[,sample_list2]
head(father_heat_plot)
#设定 annotations
annotation_col<-data.frame(Age_group = factor(c(rep("AMA", 10),rep("Young",10))))
rownames(annotation_col) = colnames(father_heat_plot)
annotation_row = data.frame(Class = factor(rep(c("UP_DMRs","Down_DMRs"), c(length(which(father_DMRs_meth$meth.diff>0)), length(which(father_DMRs_meth$meth.diff<0))))))
rownames(annotation_row) = rownames(father_DMRs_meth)
ann_colors = list(Age_group=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(UP_DMRs=ppCor[3],Down_DMRs=ppCor[4]))
p1<-pheatmap(father_heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = annotation_col, annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             gaps_row = length(which(father_DMRs_meth$meth.diff>0)),
             gaps_col =10,
             cutree_col = 2,treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","orange","firebrick3"))(50),
             main ="DNA methylation level of FF_corDMRs: original value",angle_col ="90")
#scale by row
p2<-pheatmap(father_heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = annotation_col,annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             gaps_row = length(which(father_DMRs_meth$meth.diff>0)),
             gaps_col =10,cutree_col = 2,
             treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             scale ="row", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(50),
             main ="DNA methylation level of FF_corDMRs:scale by row",angle_col ="90")

pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/father_meth_heatmap_for_FF_cor_DMR_sig_noscale.pdf")
print(p1)
dev.off()

pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/father_meth_heatmap_for_FF_cor_DMR_sig_rowscale.pdf")
print(p2)
dev.off()


###ratio for DMRs with qvalue < 0.05 in parents samples
#for mother
mother_data_merge_meth_all<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/mother_20_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth_PVALUE.csv",row.names = 1,check.names=F)
head(mother_data_merge_meth_all)
mother_data_merge_meth_all$bin_region<-unlist(lapply(strsplit(rownames(mother_data_merge_meth_all),"[.]"), function(x) paste(x[1],x[2],x[3],sep = "_")))
mother_DMRs_meth<-mother_data_merge_meth_all[which(mother_data_merge_meth_all$bin_region %in% MF_cor_DMR_sig),]

head(mother_DMRs_meth,n=3)
dim(mother_DMRs_meth)#66 27
table(mother_data_merge_meth_all$qvalue <0.05)
#FALSE   TRUE 
#554914 147023  

table(mother_DMRs_meth$qvalue <0.05)
#FALSE   TRUE 
# 15    51 
#卡方比较
data_kf1<-as.data.frame(table(mother_data_merge_meth_all$qvalue <0.05))
data_kf2<-as.data.frame(table(mother_DMRs_meth$qvalue <0.05))
data_kf<-matrix(c(rev(data_kf1$Freq),rev(data_kf2$Freq)),nrow = 2,ncol=2, byrow = T)
como <- chisq.test(data_kf)
p_HG<- ifelse(min(como$expected)>=5, round(chisq.test(data_kf,correct = F)$p.value,30),round(fisher.test(data_kf)$p.value,30))
p_HG#2.5e-29
data_kf3<-as.data.frame(data_kf)
colnames(data_kf3)<-c("sig","no_sig");rownames(data_kf3)<-c(paste0("All[n=",nrow(mother_data_merge_meth_all),"]"),paste0("Common[n=",nrow(mother_DMRs_meth),"]"))
data_kf3$type<-rownames(data_kf3)
#cell_pvalue<-data.frame(Var1=cell_name,p_value=p_value)
#cell_number_merge3<-merge(cell_number_merge2,cell_pvalue,by="Var1")
cell_number_merge4 <- melt(data_kf3,variable.name="type2",value.name = "Num",id.vars = c("type") )
cell_number_merge4$p_value<-p_HG

Pos_rate_plot<-ggplot(data=cell_number_merge4, mapping=aes(x= type,y=Num,fill=type2))+
  geom_bar(stat="identity",width=0.8,position='fill')+
  geom_text(aes(x= 1.5, y=1.1, label = round(p_value,20)), data =cell_number_merge4,vjust = 1,color="black", size=5 )+
  scale_fill_manual(name="Types",values=ppCor[c(7,5)])+
  theme_classic()+labs(x="",y="AMA significant rate",title=paste0("significant rate for mother"))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=1,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
ggsave("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/Group_compasion_for_significant_bin_rate_for_mother_all_and_MF_cor_DMR_sig.pdf",Pos_rate_plot,width=4, height=3)


#差值分布
mother_DMRs_meth$group1<-as.factor(ifelse(mother_DMRs_meth$meth.diff>0, "hyper", "hypo"))
head(mother_DMRs_meth)
plot_distribute<-ggplot(mother_DMRs_meth, aes(x=meth.diff)) +
  geom_histogram(aes(y=..density..,fill=group1), binwidth=1,colour="gray") + # 这一步很重要,使用density代替y轴
  scale_fill_manual(name="group1",values=ppCor[1:2])+
  geom_density(alpha=.2, fill="#FF6666")+ # 重叠部分采用透明设置
  geom_vline(xintercept=c(-50,-25,-15,15,25,50),linetype="dashed",colour=c("blue","blue","blue","red","red","red"))+
  theme_bw()+scale_x_continuous(breaks=seq(-55, 55, 5))+
  theme(axis.title.x = element_text(size = 12, face = "bold", vjust = 0.5,  hjust = 0.5))+
  theme(axis.title.y = element_text(size = 12,face = "bold", vjust = 0.5, hjust = 0.5))+labs(title = "methykit difference: AMA vs Young")
plot_distribute
ggsave("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/mother_meth_defference_distribution_for_MF_cor_DMR_sig.pdf",plot_distribute,width=8, height=6)

#plot heatmap
mother_DMRs_meth<-mother_DMRs_meth[order(mother_DMRs_meth$meth.diff,decreasing = T),]
mother_sample_list<-c("10-1-in-9","16-1","17-1-2D","18-1","19-1-2G","20M-4F","6-1","7-1","8-1","E16M",
                      "11-1","12-1-3G","13-1-1F","14-1","15-1","1M","2-1","3M","4-1","5-1")
mother_heat_plot<-mother_DMRs_meth[,mother_sample_list]
head(mother_heat_plot)
huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
rownames(huahua_meta)<-huahua_meta$library_code
huahua_meta2<-huahua_meta[mother_sample_list,]
head(huahua_meta2)
colnames(mother_heat_plot)<-huahua_meta2$analysis_name
as.character(huahua_meta2$analysis_name[order(huahua_meta2$analysis_name,decreasing = F)])
sample_list2 <-c(paste0("AMA_M_",1:10),paste0("YOUNG_M_",1:10))
mother_heat_plot<-mother_heat_plot[,sample_list2]
head(mother_heat_plot)
#设定 annotations
annotation_col<-data.frame(Age_group = factor(c(rep("AMA", 10),rep("Young",10))))
rownames(annotation_col) = colnames(mother_heat_plot)
annotation_row = data.frame(Class = factor(rep(c("UP_DMRs","Down_DMRs"), c(length(which(mother_DMRs_meth$meth.diff>0)), length(which(mother_DMRs_meth$meth.diff<0))))))
rownames(annotation_row) = rownames(mother_DMRs_meth)
ann_colors = list(Age_group=c(Young=ppCor[2],AMA=ppCor[1]),Class=c(UP_DMRs=ppCor[3],Down_DMRs=ppCor[4]))
p1<-pheatmap(mother_heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = annotation_col, annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             gaps_row = length(which(mother_DMRs_meth$meth.diff>0)),
             gaps_col =10,
             cutree_col = 2,treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","orange","firebrick3"))(50),
             main ="DNA methylation level of MF_corDMRs: original value",angle_col ="90")
#scale by row
p2<-pheatmap(mother_heat_plot, cluster_row =FALSE,cluster_col =FALSE,na_col = "gray",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = annotation_col,annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             gaps_row = length(which(mother_DMRs_meth$meth.diff>0)),
             gaps_col =10,cutree_col = 2,
             treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             scale ="row", 
             border=FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(50),
             main ="DNA methylation level of MF_corDMRs:scale by row",angle_col ="90")

pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/mother_meth_heatmap_for_MF_cor_DMR_sig_noscale.pdf")
print(p1)
dev.off()

pdf("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/mother_meth_heatmap_for_MF_cor_DMR_sig_rowscale.pdf")
print(p2)
dev.off()

###evaluation for parents in kids 
#kids
split_region<- "200bp"

data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t",check.names = F)
kids_DMRs<-data_region_merge_DMR_all_big$Row.names
Kids_mom_DMRs_meth<-mother_data_merge_meth_all[which(mother_data_merge_meth_all$bin_region %in% kids_DMRs),]
Kids_dad_DMRs_meth<-father_data_merge_meth_all[which(father_data_merge_meth_all$bin_region %in% kids_DMRs),]

table(mother_data_merge_meth_all$qvalue <0.05)
#FALSE   TRUE 
#554914 147023 

table(father_data_merge_meth_all$qvalue <0.05)
#FALSE   TRUE 
#490197 211740 
table(Kids_mom_DMRs_meth$qvalue <0.05)
#FALSE   TRUE 
# 121   330 
table(Kids_dad_DMRs_meth$qvalue <0.05)
#FALSE  TRUE 
#101   350
table(mother_DMRs_meth$qvalue <0.05)
#FALSE   TRUE 
# 15    51 
table(father_DMRs_meth$qvalue <0.05)
#FALSE  TRUE 
# 6    52 
#卡方比较
data_kf1<-as.data.frame(table(mother_data_merge_meth_all$qvalue <0.05))
data_kf2<-as.data.frame(table(Kids_mom_DMRs_meth$qvalue <0.05))
data_kf<-matrix(c(rev(data_kf1$Freq),rev(data_kf2$Freq)),nrow = 2,ncol=2, byrow = T)
como <- chisq.test(data_kf)
p_HG<- ifelse(min(como$expected)>=5, round(chisq.test(data_kf,correct = F)$p.value,60),round(fisher.test(data_kf)$p.value,60))
p_HG#0
data_kf3<-as.data.frame(data_kf)
colnames(data_kf3)<-c("sig","no_sig");rownames(data_kf3)<-c(paste0("Allmotherbin[n=",nrow(mother_data_merge_meth_all),"]"),paste0("kid_mom[n=",nrow(Kids_mom_DMRs_meth),"]"))
data_kf3$type<-rownames(data_kf3)
cell_number_merge4 <- melt(data_kf3,variable.name="type2",value.name = "Num",id.vars = c("type") )
cell_number_merge4$p_value<-p_HG

Pos_rate_plot<-ggplot(data=cell_number_merge4, mapping=aes(x= type,y=Num,fill=type2))+
  geom_bar(stat="identity",width=0.8,position='fill')+
  geom_text(aes(x= 1.5, y=1.1, label = round(p_value,5)), data =cell_number_merge4,vjust = 1,color="black", size=5 )+
  scale_fill_manual(name="Types",values=ppCor[c(7,5)])+
  theme_classic()+labs(x="",y="AMA significant rate",title=paste0("significant rate for kid_mother"))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=1,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
Pos_rate_plot
ggsave("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/Group_compasion_for_significant_bin_rate_for_mother_in_kid_and_mother_all_sig.pdf",Pos_rate_plot,width=4, height=3)


#卡方比较
data_kf1<-as.data.frame(table(father_data_merge_meth_all$qvalue <0.05))
data_kf2<-as.data.frame(table(Kids_dad_DMRs_meth$qvalue <0.05))
data_kf<-matrix(c(rev(data_kf1$Freq),rev(data_kf2$Freq)),nrow = 2,ncol=2, byrow = T)
como <- chisq.test(data_kf)
p_HG<- ifelse(min(como$expected)>=5, round(chisq.test(data_kf,correct = F)$p.value,60),round(fisher.test(data_kf)$p.value,60))
p_HG#0
data_kf3<-as.data.frame(data_kf)
colnames(data_kf3)<-c("sig","no_sig");rownames(data_kf3)<-c(paste0("Allfatherbin[n=",nrow(father_data_merge_meth_all),"]"),paste0("kid_dad[n=",nrow(Kids_dad_DMRs_meth),"]"))
data_kf3$type<-rownames(data_kf3)
cell_number_merge4 <- melt(data_kf3,variable.name="type2",value.name = "Num",id.vars = c("type") )
cell_number_merge4$p_value<-p_HG

Pos_rate_plot<-ggplot(data=cell_number_merge4, mapping=aes(x= type,y=Num,fill=type2))+
  geom_bar(stat="identity",width=0.8,position='fill')+
  geom_text(aes(x= 1.5, y=1.1, label = round(p_value,5)), data =cell_number_merge4,vjust = 1,color="black", size=5 )+
  scale_fill_manual(name="Types",values=ppCor[c(7,5)])+
  theme_classic()+labs(x="",y="AMA significant rate",title=paste0("significant rate for kid_father"))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=1,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
Pos_rate_plot
ggsave("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/Group_compasion_for_significant_bin_rate_for_father_in_kid_and_father_all_sig.pdf",Pos_rate_plot,width=4, height=3)


#卡方比较
data_kf1<-as.data.frame(table(Kids_dad_DMRs_meth$qvalue <0.05))
data_kf2<-as.data.frame(table(father_DMRs_meth$qvalue <0.05))
data_kf<-matrix(c(rev(data_kf1$Freq),rev(data_kf2$Freq)),nrow = 2,ncol=2, byrow = T)
como <- chisq.test(data_kf)
p_HG<- ifelse(min(como$expected)>=5, round(chisq.test(data_kf,correct = F)$p.value,30),round(fisher.test(data_kf)$p.value,30))
p_HG#0.0340052
data_kf3<-as.data.frame(data_kf)
colnames(data_kf3)<-c("sig","no_sig");rownames(data_kf3)<-c(paste0("Allkids[n=",nrow(Kids_dad_DMRs_meth),"]"),paste0("Common[n=",nrow(father_DMRs_meth),"]"))
data_kf3$type<-rownames(data_kf3)
cell_number_merge4 <- melt(data_kf3,variable.name="type2",value.name = "Num",id.vars = c("type") )
cell_number_merge4$p_value<-p_HG

Pos_rate_plot<-ggplot(data=cell_number_merge4, mapping=aes(x= type,y=Num,fill=type2))+
  geom_bar(stat="identity",width=0.8,position='fill')+
  geom_text(aes(x= 1.5, y=1.1, label = round(p_value,5)), data =cell_number_merge4,vjust = 1,color="black", size=5 )+
  scale_fill_manual(name="Types",values=ppCor[c(7,5)])+
  theme_classic()+labs(x="",y="AMA significant rate",title=paste0("significant rate for father"))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=1,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
Pos_rate_plot
ggsave("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/Group_compasion_for_significant_bin_rate_for_father_in_kid_and_FF_cor_DMR_sig_2.pdf",Pos_rate_plot,width=4, height=3)


#卡方比较
data_kf1<-as.data.frame(table(Kids_mom_DMRs_meth$qvalue <0.05))
data_kf2<-as.data.frame(table(mother_DMRs_meth$qvalue <0.05))
data_kf<-matrix(c(rev(data_kf1$Freq),rev(data_kf2$Freq)),nrow = 2,ncol=2, byrow = T)
como <- chisq.test(data_kf)
p_HG<- ifelse(min(como$expected)>=5, round(chisq.test(data_kf,correct = F)$p.value,30),round(fisher.test(data_kf)$p.value,30))
p_HG#0.4796177
data_kf3<-as.data.frame(data_kf)
colnames(data_kf3)<-c("sig","no_sig");rownames(data_kf3)<-c(paste0("Allkids[n=",nrow(Kids_mom_DMRs_meth),"]"),paste0("Common[n=",nrow(mother_DMRs_meth),"]"))
data_kf3$type<-rownames(data_kf3)
cell_number_merge4 <- melt(data_kf3,variable.name="type2",value.name = "Num",id.vars = c("type") )
cell_number_merge4$p_value<-p_HG

Pos_rate_plot<-ggplot(data=cell_number_merge4, mapping=aes(x= type,y=Num,fill=type2))+
  geom_bar(stat="identity",width=0.8,position='fill')+
  geom_text(aes(x= 1.5, y=1.1, label = round(p_value,5)), data =cell_number_merge4,vjust = 1,color="black", size=5 )+
  scale_fill_manual(name="Types",values=ppCor[c(7,5)])+
  theme_classic()+labs(x="",y="AMA significant rate",title=paste0("significant rate for mother"))+
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=1,hjust=0.5,angle = 0),legend.title = element_text(size = 9))
Pos_rate_plot
ggsave("/mnt/data/chenwei/huahua/4.methy_result/3.cor_DMRs/Group_compasion_for_significant_bin_rate_for_mother_in_kid_and_MF_cor_DMR_sig_2.pdf",Pos_rate_plot,width=4, height=3)