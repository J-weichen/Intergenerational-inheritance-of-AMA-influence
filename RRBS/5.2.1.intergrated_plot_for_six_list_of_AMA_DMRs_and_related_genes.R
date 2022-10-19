rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ggpubr)
library(scales)
library(ggsci)
library(stringr)
library(UpSetR)
library(pheatmap)
library(VennDiagram)

pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#甲基化对应基因列表的读取
#data preparation
#15% df
compare_name1<-"kids_AMA_vs_Young"
split_region<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
hyper_gene1<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
hypo_gene1<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))

DMR_hyper1<-unique(as.character(DMR_hyper$Row.names))
DMR_hypo1<-unique(as.character(DMR_hypo$Row.names))
length(hyper_gene1);length(hypo_gene1)#170 236
length(DMR_hyper1);length(DMR_hypo1)#182 #269

compare_name2<-"mother_AMA_vs_Young"
split_region2<- "200bp"
mom_data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name2,"_",split_region2,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper<-mom_data_region_merge_DMR_all_big[mom_data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-mom_data_region_merge_DMR_all_big[mom_data_region_merge_DMR_all_big$meth.diff <0,]
hyper_gene3<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
hypo_gene3<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))

DMR_hyper3<-unique(as.character(DMR_hyper$Row.names))
DMR_hypo3<-unique(as.character(DMR_hypo$Row.names))
length(hyper_gene3);length(hypo_gene3)#170 253
length(DMR_hyper3);length(DMR_hypo3)# 236 #280

compare_name2<-"father_AMA_vs_Young"
split_region2<- "200bp"
dad_data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name2,"_",split_region2,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper<-dad_data_region_merge_DMR_all_big[dad_data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-dad_data_region_merge_DMR_all_big[dad_data_region_merge_DMR_all_big$meth.diff <0,]
hyper_gene4<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
hypo_gene4<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))

DMR_hyper4<-unique(as.character(DMR_hyper$Row.names))
DMR_hypo4<-unique(as.character(DMR_hypo$Row.names))
length(hyper_gene4);length(hypo_gene4)#185  238
length(DMR_hyper4);length(DMR_hypo4)#206 #258

#upset plot for all six lists of DMRs
listinput<-c(list(DMR_hyper1),list(DMR_hyper3),list(DMR_hyper4),list(DMR_hypo1),list(DMR_hypo3),list(DMR_hypo4))
names(listinput)<-c("Kids_hyper","Mother_hyper","Father_hyper","Kids_hypo","Mother_hypo","Father_hypo" )
upset_dataframe<-as.data.frame(fromList(listinput))
dim(upset_dataframe)
upset_dataframe[1:10,1:5]
upset_dataframe_rowSums<-rowSums(upset_dataframe)
range(upset_dataframe_rowSums)#1:3
upset_dataframe_colSums<-colSums(upset_dataframe)
range(upset_dataframe_colSums)#182 280
upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]
upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400)

list2_q1<-list(query = intersects, params = list("Kids_hyper", "Mother_hyper","Father_hyper"), color = "red", active = T,query.name = "common hyper DMRs")
list2_q2<-list(query = intersects, params = list("Kids_hypo","Mother_hypo","Father_hypo"), color = "navy", active = T,query.name = "common hypo DMRs")
list2_q3<-list(query = intersects, params = list("Kids_hyper", "Mother_hyper"), color = "orange", active = T,query.name = "mother_SF_DMRs")
list2_q4<-list(query = intersects, params = list("Kids_hypo","Mother_hypo"), color = "orange", active = T,query.name = "mother_SF_DMRs")
list2_q5<-list(query = intersects, params = list("Kids_hyper", "Father_hyper"), color = "purple", active = T,query.name = "father_SF_DMRs")
list2_q6<-list(query = intersects, params = list("Kids_hypo","Father_hypo"), color = "purple", active = T,query.name = "father_SF_DMRs")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/upsetplot_six_lists_DMRs_AMA_vs_Young_200bin_q005_dif15.pdf",width = 10,height=8)

upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400,
      sets=c("Kids_hyper","Mother_hyper","Father_hyper","Kids_hypo","Mother_hypo","Father_hypo"),
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
      queries = list(list2_q1,list2_q2,list2_q3,list2_q4,list2_q5,list2_q6)
)
dev.off()

#upset plot for all six lists of gene of DMRs
listinput<-c(list(hyper_gene1),list(hyper_gene3),list(hyper_gene4),list(hypo_gene1),list(hypo_gene3),list(hypo_gene4))
names(listinput)<-c("Kids_hyper","Mother_hyper","Father_hyper","Kids_hypo","Mother_hypo","Father_hypo" )
upset_dataframe<-as.data.frame(fromList(listinput))
dim(upset_dataframe)
upset_dataframe[1:10,1:5]
upset_dataframe_rowSums<-rowSums(upset_dataframe)
range(upset_dataframe_rowSums)#1:4
upset_dataframe_colSums<-colSums(upset_dataframe)
range(upset_dataframe_colSums)#170 253
upset_dataframe_colSums[which(upset_dataframe_colSums== 1)]
upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400)

list2_q1<-list(query = intersects, params = list("Kids_hyper", "Mother_hyper","Father_hyper"), color = "red", active = T,query.name = "common hyper genes")
list2_q2<-list(query = intersects, params = list("Kids_hypo","Mother_hypo","Father_hypo"), color = "navy", active = T,query.name = "common hypo genes")
list2_q3<-list(query = intersects, params = list("Kids_hyper", "Mother_hyper"), color = "orange", active = T,query.name = "mother_SF_genes")
list2_q4<-list(query = intersects, params = list("Kids_hypo","Mother_hypo"), color = "orange", active = T,query.name = "mother_SF_genes")
#list2_q5<-list(query = intersects, params = list("Kids_hyper", "Father_hyper"), color = "purple", active = T,query.name = "father_SF_genes")
list2_q5<-list(query = intersects, params = list("Kids_hypo","Father_hypo"), color = "purple", active = T,query.name = "father_SF_genes")
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/upsetplot_six_lists_genes_near_DMRs_AMA_vs_Young_200bin_q005_dif15.pdf",width = 10,height=8)

upset(upset_dataframe, main.bar.color = "black",nsets = length(colnames(upset_dataframe)), nintersects = 400,
      sets=c("Kids_hyper","Mother_hyper","Father_hyper","Kids_hypo","Mother_hypo","Father_hypo"),
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
      queries = list(list2_q1,list2_q2,list2_q3,list2_q4,list2_q5)
)
dev.off()
#For DMRs venn
##for four list of mom and father 
tag<- "four list of mom and dad DMRs"
mainname<-tag
venn2 <-venn.diagram(list(Mother_hypo=DMR_hypo3,Mother_hyper=DMR_hyper3,Father_hypo= DMR_hypo4,Father_hyper= DMR_hyper4),
                     alpha=c(0.5,0.5,0.5,0.5),
                     lwd=1,lty=1,col="white",fill=ppCor[c(5,8,6,10)], 
                     cex = 1.5,cat.col=ppCor[c(5,8,6,10)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=mainname, 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)
grid.newpage();grid.draw(venn2)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_four_lists_mom_data_DMRs_AMA_vs_Young_200bin_q005_dif15.pdf",width = 8,height=8)
grid.draw(venn2)
dev.off()
simple_DMR_hyper_common<-setdiff(Reduce(intersect,list(DMR_hyper3,DMR_hyper4)),c(DMR_hypo3,DMR_hypo4))
simple_DMR_hypo_common<-setdiff(Reduce(intersect,list(DMR_hypo3,DMR_hypo4)),c(DMR_hyper3,DMR_hyper4))
dad_data_region_merge_DMR_all_big[which(dad_data_region_merge_DMR_all_big$Row.names %in% simple_DMR_hyper_common),c("annotation","SYMBOL")]
dad_data_region_merge_DMR_all_big[which(dad_data_region_merge_DMR_all_big$Row.names %in% simple_DMR_hypo_common),c("annotation","SYMBOL")]


#for all four group 
tag<- "Three list of hyper DMRs"
mainname<-tag
venn1 <-venn.diagram(list(Mother=DMR_hyper3,Father= DMR_hyper4,Kids=DMR_hyper1),
                     alpha=c(0.5,0.5,0.5),
                     lwd=1,lty=1,col="white",fill=ppCor[c(5:6,10)], 
                     cex = 1.5,cat.col=ppCor[c(5:6,10)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=mainname, 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)

tag<- "Three list of hypo DMRs"
mainname<-tag
venn2 <-venn.diagram(list(Mother=DMR_hypo3,Father= DMR_hypo4,Kids=DMR_hypo1),
                     alpha=c(0.5,0.5,0.5),
                     lwd=1,lty=1,col="white",fill=ppCor[c(5:6,10)], 
                     cex = 1.5,cat.col=ppCor[c(5:6,10)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=mainname, 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_three_lists_DMRs_AMA_vs_Young_200bin_q005_dif15.pdf",width = 8,height=8)
grid.newpage(); 
grid.draw(venn1)
grid.newpage();
grid.draw(venn2)
dev.off()

DMR_hyper_arrest_common<-Reduce(intersect,list(DMR_hyper1,DMR_hyper3,DMR_hyper4))
DMR_hypo_arrest_common<-Reduce(intersect,list(DMR_hypo1,DMR_hypo3,DMR_hypo4))

#For genes of DMRs venn
##DMR_hyper
tag<- "Three list of genes of hyper DMRs"
mainname<-tag
grid.newpage(); #清空画板，开始画新图
venn3 <-venn.diagram(list(Mother=hyper_gene3,Father= hyper_gene4,Kids=hyper_gene1),
                     alpha=c(0.5,0.5,0.5),
                     lwd=1,lty=1,col="white",fill=ppCor[c(3:4,9)], 
                     cex = 1.5,cat.col=ppCor[c(3:4,9)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=mainname, 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)

##DMR_hypo
tag<- "Three list of genes of hypo DMRs"
mainname<-tag
grid.newpage(); #清空画板，开始画新图
venn4 <-venn.diagram(list(Mother=hypo_gene3,Father= hypo_gene4,Kids=hypo_gene1),
                     alpha=c(0.5,0.5,0.5),
                     lwd=1,lty=1,col="white",fill=ppCor[c(3:4,9)], 
                     cex = 1.5,cat.col=ppCor[c(3:4,9)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                     cat.fontface=4, cat.cex = 1.5, main=mainname, 
                     main.cex = 2, main.fontface = 2, main.fontfamily = 3, filename = NULL)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_three_lists_genes_associated_DMRs_AMA_vs_Young_200bin_q005_dif15.pdf",width = 8,height=8)

grid.newpage()
grid.draw(venn3)
grid.newpage()
grid.draw(venn4)
dev.off()

gene_DMR_hyper_arrest_common_three<-Reduce(intersect,list(hyper_gene1,hyper_gene3,hyper_gene4))
gene_DMR_hypo_arrest_common_three<-Reduce(intersect,list(hypo_gene1,hypo_gene3,hypo_gene4))

##以下未跑
#enrichment for gene of common DMRs
DMR_hyper_arrest_common<-Reduce(intersect,list(DMR_hyper1,DMR_hyper3,DMR_hyper4))
DMR_hypo_arrest_common<-Reduce(intersect,list(DMR_hypo1,DMR_hypo3,DMR_hypo4))
data_DMR_hyper_arrest_common<-data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$bin_region %in% DMR_hyper_arrest_common),]
data_DMR_hypo_arrest_common<-data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$bin_region %in% DMR_hypo_arrest_common),]
write.table(data_DMR_hyper_arrest_common, file="/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/three_DMR_hyper_arrest_common_annotation.txt",row.names=T, col.names=T) 
write.table(data_DMR_hypo_arrest_common, file="/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/three_DMR_hypo_arrest_common_annotation.txt",row.names=T, col.names=T) 

split_region<-"200bp"
#for hyper
target<-unique(data_DMR_hyper_arrest_common$geneId)
str(target)
tender<-"hyper"
gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
head(as.data.frame(BP@result))
CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.1,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))

write.table(as.data.frame(CC@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_CC_GO.txt", row.names=T, col.names=T)) 
write.table(as.data.frame(BP@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_BP_GO.txt",row.names=T, col.names=T)) 
write.table(as.data.frame(MF@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_MF_GO.txt",row.names=T, col.names=T)) 
write.table(as.data.frame(CC_simp@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_CC_GO_simp.txt",row.names=T, col.names=T)) 
write.table(as.data.frame(BP_simp@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_BP_GO_simp.txt",row.names=T, col.names=T)) 
write.table(as.data.frame(MF_simp@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_MF_GO_simp.txt", row.names=T, col.names=T)) 
head(as.data.frame(BP_simp@result))
write.table((as.data.frame(kk@result)), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_KEGG.txt",row.names=T, col.names=T)) 
write.table((as.data.frame(rpea@result)), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_RPEA.txt",row.names=T, col.names=T)) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_GO_RPEA_result.rds"))

#for hypo
target<-unique(data_DMR_hypo_arrest_common$geneId)
str(target)
tender<-"hypo"

gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T) 
head(as.data.frame(BP@result));head(as.data.frame(MF@result));head(as.data.frame(CC@result))

CC_simp <-  clusterProfiler::simplify(CC, cutoff=0.7,by="p.adjust",select_fun=min) 
BP_simp <-  clusterProfiler::simplify(BP, cutoff=0.7,by="p.adjust",select_fun=min) 
MF_simp <-  clusterProfiler::simplify(MF, cutoff=0.7,by="p.adjust",select_fun=min)
head(as.data.frame(BP_simp@result));head(as.data.frame(MF_simp@result));head(as.data.frame(CC_simp@result))
kk <- clusterProfiler::enrichKEGG(gene = target,organism ='hsa',pvalueCutoff = 0.05, qvalueCutoff = 0.05,minGSSize = 2,use_internal_data =TRUE)
kk<-clusterProfiler::setReadable(kk,org.Hs.eg.db, keyType="ENTREZID")
head(as.data.frame(kk@result))
## Reactome pathway enrichment analysis
rpea <- enrichPathway(gene=target,pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05, readable=T,minGSSize =2,organism = "human")
head(as.data.frame(rpea@result))

write.table(as.data.frame(CC@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_CC_GO.txt", row.names=T, col.names=T)) 
write.table(as.data.frame(BP@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_BP_GO.txt",row.names=T, col.names=T)) 
write.table(as.data.frame(MF@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_MF_GO.txt",row.names=T, col.names=T)) 
write.table(as.data.frame(CC_simp@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_CC_GO_simp.txt", row.names=T, col.names=T)) 
write.table(as.data.frame(BP_simp@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_BP_GO_simp.txt",row.names=T, col.names=T)) 
write.table(as.data.frame(MF_simp@result), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_MF_GO_simp.txt",row.names=T, col.names=T)) 
head(as.data.frame(BP_simp@result))
write.table((as.data.frame(kk@result)), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_KEGG.txt", row.names=T, col.names=T)) 
write.table((as.data.frame(rpea@result)), file=paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_RPEA.txt",row.names=T, col.names=T)) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_GO_RPEA_result.rds"))

#plot 
gene_erichment_results_hyper<-readRDS(file = paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hyper_",split_region,"_15q005_GO_RPEA_result.rds"))
gene_erichment_results_hypo<-readRDS(file = paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_hypo_",split_region,"_15q005_GO_RPEA_result.rds"))

KEGG_plot_hyper<-gene_erichment_results_hyper[[7]] 
KEGG_plot_hypo<-gene_erichment_results_hypo[[7]]
dim(KEGG_plot_hyper@result);dim(KEGG_plot_hypo@result)
KEGG_plot<-as.data.frame(rbind(KEGG_plot_hyper@result,KEGG_plot_hypo@result))
KEGG_plot$group<-c(rep("hyper",nrow(KEGG_plot_hyper@result)),rep("hypo",nrow(KEGG_plot_hypo@result)))
str(KEGG_plot)
KEGG_plot$pvalue_log10<-c(-log(KEGG_plot$pvalue,10))

KEGG_plot$group<- factor(x =KEGG_plot$group,levels =c("hyper","hypo") )
levels(KEGG_plot$group)
KEGG_plot<-KEGG_plot[,c("Description","group","pvalue_log10","Count")]
KEGG_plot$Description<-as.character(KEGG_plot$Description)
KEGG_plot<-KEGG_plot[order(KEGG_plot$group,KEGG_plot$pvalue_log10),]
str(KEGG_plot);head(KEGG_plot)

KEGG_plot2<-KEGG_plot[which(KEGG_plot$Count >=2),]
head(KEGG_plot2)

pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_",split_region,"_DMRs_15q005_hypo_hyper_KEGG.pdf"))
ggplot(KEGG_plot2,aes(x=group,y=Description,size=Count,colour=pvalue_log10))+
  geom_point(alpha =0.8,na.rm = TRUE)+
  scale_size(breaks = c(1,2,3,5),range = c(3,7),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,1,2,3,4,5),
                        name='-log10(p_value)')+  
  scale_y_discrete(labels=function(x) str_wrap(x, width=150))+
  theme_classic()+labs(x="",y="GO terms",title="Arrest DMRs ：：KEGG enrichment")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),
        legend.position ="right",legend.direction = "vertical")

dev.off()

GO_BP_plot_hyper<-gene_erichment_results_hyper[[2]] 
GO_BP_plot_hypo<-gene_erichment_results_hypo[[2]]
dim(GO_BP_plot_hyper@result);dim(GO_BP_plot_hypo@result)
GO_BP_plot<-as.data.frame(rbind(GO_BP_plot_hyper@result,GO_BP_plot_hypo@result))
GO_BP_plot$group<-c(rep("hyper",nrow(GO_BP_plot_hyper@result)),rep("hypo",nrow(GO_BP_plot_hypo@result)))
str(GO_BP_plot)
GO_BP_plot$pvalue_log10<-c(-log(GO_BP_plot$pvalue,10))

GO_BP_plot$group<- factor(x =GO_BP_plot$group,levels =c("hyper","hypo") )
levels(GO_BP_plot$group)
GO_BP_plot<-GO_BP_plot[,c("Description","group","pvalue_log10","Count")]
GO_BP_plot$Description<-as.character(GO_BP_plot$Description)
GO_BP_plot<-GO_BP_plot[order(GO_BP_plot$group,GO_BP_plot$pvalue_log10),]
str(GO_BP_plot);head(GO_BP_plot)

GO_BP_plot2<-GO_BP_plot[which(GO_BP_plot$Count >2),]
head(GO_BP_plot2);dim(GO_BP_plot2)
range(GO_BP_plot$Count)
pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/enrichment/","Three_common_DMRs_",split_region,"_DMRs_15q005_hypo_hyper_GO_BP_3gene.pdf"))
ggplot(GO_BP_plot2,aes(x=group,y=Description,size=Count,colour=pvalue_log10))+
  geom_point(alpha =0.8,na.rm = TRUE)+
  scale_size(breaks = c(1,3,5,7,9),range = c(3,8),name='Genes number')+
  scale_color_gradientn(colors = rev(brewer.pal(7,'RdYlBu')),breaks=c(0,1,2,3,4),
                        name='-log10(p_value)')+  
  scale_y_discrete(labels=function(x) str_wrap(x, width=150))+
  theme_classic()+labs(x="",y="GO terms",title="Arrest DMRs ：：GO BP enrichment")+
  theme(legend.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10,colour = 'black'),
        panel.grid.major.y = element_line(color="grey", size= 0.1),
        axis.text.x = element_text(size = 10,colour = 'black',vjust=0.5,hjust=0.5,angle = 0),
        legend.title = element_text(size = 10),
        legend.position ="right",legend.direction = "vertical")

dev.off()
#以下未跑
#plot heatmap  for selected DMRs or genes
compare_name1<-"All_arrest_vs_abortion"
split_region1<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name1,"_",split_region1,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
data_region_merge_DMR_all_big[1:3,]
#for common DMRs
data_DMR_hyper_arrest_common<-data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$bin_region %in% DMR_hyper_arrest_common),]
data_DMR_hypo_arrest_common<-data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$bin_region %in% DMR_hypo_arrest_common),]
data_DMR_arrest_common<-rbind(data_DMR_hyper_arrest_common,data_DMR_hypo_arrest_common)

data_DMR_arrest_common$SYMBOL[duplicated(data_DMR_arrest_common$SYMBOL)]#28
length(as.character(data_DMR_hyper_arrest_common$SYMBOL[duplicated(data_DMR_hyper_arrest_common$SYMBOL)]))#8
#"ANKRD9"  "PLEKHH3" "PLEKHH3" "PLEKHH3" "PLEKHH3" "HMG20B"  "PRAG1"   "ERI1" 
length(as.character(data_DMR_hypo_arrest_common$SYMBOL[duplicated(data_DMR_hypo_arrest_common$SYMBOL)]))#20
# [1] "UBE4B"        "STK40"        "STK40"        "LINC00592"    "PIBF1"        "LOC101927911"
# [7] "SUZ12P1"      "TK1"          "SPPL2B"       "LAMA5"        "LAMA5"        "PKNOX1"      
# [13] "TBC1D10A"     "SH3BP2"       "UBE2K"        "ATP10D"       "COQ2"         "CCND3"       
# [19] "SUN1"         "DENND4C"

Young_Abortion_sample<-c("P11","P12","P13","P14","P15","P31","P32","P33","P34","P35","P36","P41")
Young_Arrest_sample<-c("A1","A3","A4","A22","A25","A31","A32","A33")
AMA_Abortion_sample<-c("P16","P17","P18","P19","P20","P21","P22","P24","P25","P38")
AMA_Arrest_sample<-c("A6","A8","A9","A10","A11","A36")
rownames(data_DMR_arrest_common)<-data_DMR_arrest_common$Row.names
plot_DMRs2<-data_DMR_arrest_common[,c(Young_Abortion_sample,Young_Arrest_sample,AMA_Abortion_sample,AMA_Arrest_sample)]
pheatmap(plot_DMRs2,scale ="row",cluster_row =FALSE,cluster_col =FALSE,na_col = "grey90", clustering_method ="complete")

#设定 annotations
# 生成行 列的注释
#for DNA methylation
annotation_col<-data.frame(Treatment = factor(c(rep("Young_Abortion", 12),rep("Young_Arrest", 8),rep("AMA_Abortion", 10),rep("AMA_Arrest", 6))))
rownames(annotation_col) = colnames(plot_DMRs2)

annotation_row = data.frame( Class = data_DMR_arrest_common$SYMBOL)
rownames(annotation_row) = rownames(plot_DMRs2)

## 自定义分组的颜色
ann_colors = list(Treatment=c(Young_Abortion=ppCor[6],Young_Arrest=ppCor[4],AMA_Abortion=ppCor[2],AMA_Arrest=ppCor[8]))

grid.newpage(); #清空画板，开始画新图
#scale by row
p1<-pheatmap(plot_DMRs2,cluster_row =FALSE,cluster_col =FALSE,na_col = "grey",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = annotation_col,#annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             # gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
             treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             # color = colorRampPalette(c("navy","pink","orange","firebrick3"))(50),
             border=FALSE,
             scale = "row",
             #   color = colorRampPalette(c("purple","grey","orange"))(40),
             #color = colorRampPalette(c("purple","black","gold"))(40),
             color = colorRampPalette(c("navy","white","firebrick3"))(50),
             #  filename = "/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/DMR_methlevel_genes_over_DEGs_single_heatmap_scale_by_row4.pdf",
             main ="DNA methylation level of selected  DMRs :scale by row",angle_col ="90")

#no scale by row
p2<-pheatmap(plot_DMRs2,cluster_row =FALSE,cluster_col =FALSE,na_col = "black",
             clustering_distance_rows ="euclidean",#correlation
             show_rownames = F,show_colnames = T,
             annotation_col = annotation_col,#annotation_row=annotation_row,
             annotation_colors = ann_colors, 
             # gaps_row = length(which(data_region_merge_DMR_all_simple$meth.diff>0)),gaps_col =c(3),cutree_col = 2,
             treeheight_col = 20, #treeheight_row = 30, 
             # labels_row = labels_row,
             #border_color ="red", 
             #color = colorRampPalette(c("navy","pink","orange","firebrick3"))(50),
             border=FALSE,
             #    color = colorRampPalette(c("purple","grey","orange"))(40),
             color = colorRampPalette(c("purple","darkgrey","gold"))(40),
             #    color = colorRampPalette(c("navy","white","firebrick3"))(50),
             #  filename = "/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/DMR_methlevel_genes_over_DEGs_single_heatmap_scale_by_row4.pdf",
             main ="DNA methylation level of selected  DMRs :no scale by row",angle_col ="90")
pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/four_group_common_DMR_6X_CpGs_200bp_DMRs_myDiff15q005_heatmap1.pdf"))
print(p1)
dev.off()
pdf(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/four_group_common_DMR_6X_CpGs_200bp_DMRs_myDiff15q005_heatmap2.pdf"))
print(p2)
dev.off()
dev.new()
