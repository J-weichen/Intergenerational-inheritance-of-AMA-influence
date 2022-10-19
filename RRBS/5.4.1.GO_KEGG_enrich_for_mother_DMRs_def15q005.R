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

#for kids

#def 15%
compare_name<-"mother_AMA_vs_Young"
split_region<- "200bp"
depth<-"6X"

# no select
annotation_DMR_all  <-read.csv(paste0("D:/huahua/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_all_myDiff15q005.csv"))
annotation_DMR_hyper<-read.csv(paste0("D:/huahua/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hyper_myDiff15q005.csv"))
annotation_DMR_hypo <-read.csv(paste0("D:/huahua/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hypo_myDiff15q005.csv"))

#within 3kb
annotation_DMR_all_3kb<-read.csv(paste0("D:/huahua/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_all_myDiff15q005_within_TSS_3kb.csv"))
annotation_DMR_hyper_3kb<-read.csv(paste0("D:/huahua/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hyper_myDiff15q005_within_TSS_3kb.csv"))
annotation_DMR_hypo_3kb<-read.csv(paste0("D:/huahua/1.DMR_called/",compare_name,"_annotation_for_",split_region,"_DMRs_hypo_myDiff15q005_within_TSS_3kb.csv"))

##make gene list 
#gene within 3kb
gene_3kb<-list(DMRs_all=annotation_DMR_all_3kb$geneId,DMRs_hyper=annotation_DMR_hyper_3kb$geneId,DMRs_hypo=annotation_DMR_hypo_3kb$geneId)
#unselected
gene_no_select <- list(DMRs_all=annotation_DMR_all$geneId,DMRs_hyper=annotation_DMR_hyper$geneId,DMRs_hypo=annotation_DMR_hypo$geneId)

##enrichment analysis for special list
##for genes with 3kb
gene<-gene_3kb #gene_no_select
#for hyper
target<-unique(gene$DMRs_hyper)
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

write.csv(as.data.frame(CC@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_CC_GO.csv"))
write.csv(as.data.frame(BP@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_BP_GO.csv")) 
write.csv(as.data.frame(MF@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_MF_GO.csv")) 
write.csv(as.data.frame(CC_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_CC_GO_simp.csv"))
write.csv(as.data.frame(BP_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_BP_GO_simp.csv")) 
write.csv(as.data.frame(MF_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_MF_GO_simp.csv")) 
head(as.data.frame(BP_simp@result))
write.csv((as.data.frame(kk@result)), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_KEGG.csv") )
write.csv((as.data.frame(rpea@result)), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_RPEA.csv")) 

gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_GO_RPEA_result.rds"))

#for hypo
target<-unique(gene$DMRs_hypo)
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

write.csv(as.data.frame(CC@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_CC_GO.csv"))
write.csv(as.data.frame(BP@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_BP_GO.csv")) 
write.csv(as.data.frame(MF@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_MF_GO.csv")) 
write.csv(as.data.frame(CC_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_CC_GO_simp.csv")) 
write.csv(as.data.frame(BP_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_BP_GO_simp.csv")) 
write.csv(as.data.frame(MF_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_MF_GO_simp.csv")) 
head(as.data.frame(BP_simp@result))
write.csv((as.data.frame(kk@result)), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_KEGG.csv")) 
write.csv((as.data.frame(rpea@result)), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_RPEA.csv")) 

gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_GO_RPEA_result.rds"))


####for unselect genes
gene<-gene_no_select
#for hyper
target<-unique(gene$DMRs_hyper)
str(target)
tender<-"hyper"
gene_GO_RPEA_erichment_results<-list()
BP <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T)  
MF <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
CC <- enrichGO(target,"org.Hs.eg.db", keyType = "ENTREZID",ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1, readable=T) 
head(summary(BP))
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

write.csv(as.data.frame(CC@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_CC_GO.csv"))
write.csv(as.data.frame(BP@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_BP_GO.csv")) 
write.csv(as.data.frame(MF@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_MF_GO.csv")) 
write.csv(as.data.frame(CC_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_CC_GO_simp.csv")) 
write.csv(as.data.frame(BP_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_BP_GO_simp.csv")) 
write.csv(as.data.frame(MF_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_MF_GO_simp.csv")) 
head(as.data.frame(BP_simp@result))
write.csv((as.data.frame(kk@result)), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_KEGG.csv")) 
write.csv((as.data.frame(rpea@result)), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_RPEA.csv")) 


gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_GO_RPEA_result.rds"))

#for hypo
target<-unique(gene$DMRs_hypo)
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
write.csv(as.data.frame(CC@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_CC_GO.csv"))
write.csv(as.data.frame(BP@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_BP_GO.csv")) 
write.csv(as.data.frame(MF@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_MF_GO.csv")) 
write.csv(as.data.frame(CC_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_CC_GO_simp.csv")) 
write.csv(as.data.frame(BP_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_BP_GO_simp.csv")) 
write.csv(as.data.frame(MF_simp@result), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_MF_GO_simp.csv")) 
head(as.data.frame(BP_simp@result))
write.csv((as.data.frame(kk@result)), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_KEGG.csv")) 
write.csv((as.data.frame(rpea@result)), file=paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_RPEA.csv")) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_GO_RPEA_result.rds"))

#结果展示 https://www.jianshu.com/p/e133ab3169fa
gene_erichment_results_hyper_3kb<-readRDS(file = paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_3kb_GO_RPEA_result.rds"))
gene_erichment_results_hypo_3kb<-readRDS(file = paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_3kb_GO_RPEA_result.rds"))

gene_erichment_results_hyper_unselect<-readRDS(file = paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_15q005_unselect_GO_RPEA_result.rds"))
gene_erichment_results_hypo_unselct<-readRDS(file = paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_15q005_unselect_GO_RPEA_result.rds"))

#for genes in 3kb
#For hyper
gene_erichment_results<-gene_erichment_results_hyper_3kb #gene_erichment_results_hyper_unselect
names(gene_erichment_results)
CC_plot<-gene_erichment_results[[1]]
BP_plot<-gene_erichment_results[[2]]
MF_plot<-gene_erichment_results[[3]]

#结果展示
merge_three<-as.data.frame(rbind(head(BP_plot@result,n=10),head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
merge_three$qvalue_log10<-c(-log(merge_three$qvalue,10))
merge_three$pvalue_log10<-c(-log(merge_three$pvalue,10))
dim(merge_three)
merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))

pdf(paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_DMRs_15q005_enrichment_BP_CC_MF_3kb.pdf"))
ggbarplot(merge_three, x="Description", y="qvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          #          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=90)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55)) #按组排序 x.text.angle=90
#修进：http://blog.sciencenet.cn/blog-3334560-1091714.html
ggbarplot(merge_three, x="Description", y="pvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          #          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=90)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55)) #按组排序 x.text.angle=90
#修进：http://blog.sciencenet.cn/blog-3334560-1091714.html
dev.off()

#For hypo
gene_erichment_results<-gene_erichment_results_hypo_3kb#gene_erichment_results_hypo_unselct
CC_plot<-gene_erichment_results[[1]]
BP_plot<-gene_erichment_results[[2]]
MF_plot<-gene_erichment_results[[3]]

#结果展示
merge_three<-as.data.frame(rbind(head(BP_plot@result,n=10),head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
dim(merge_three)
#merge_three$qlogtran<-c(-log(merge_three$qvalue,10))
#merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))
#merge_two<-as.data.frame(rbind(head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
merge_three$qvalue_log10<-c(-log(merge_three$qvalue,10))
merge_three$pvalue_log10<-c(-log(merge_three$pvalue,10))
merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))
pdf(paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_DMRs_15q005_enrichment_BP_CC_MF_3kb.pdf"))
ggbarplot(merge_three, x="Description", y="qvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          #          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=90)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55)) #按组排序 x.text.angle=90
#修进：http://blog.sciencenet.cn/blog-3334560-1091714.html
ggbarplot(merge_three, x="Description", y="pvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          #          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=60)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55))
dev.off()



KEGG_plot_hyper<-gene_erichment_results_hyper_3kb[[7]] #gene_erichment_results_hyper_unselect
KEGG_plot_hypo<-gene_erichment_results_hypo_3kb[[7]] #gene_erichment_results_hypo_unselect

dim(KEGG_plot_hyper@result)
dim(KEGG_plot_hypo@result)

KEGG_plot<-as.data.frame(rbind(head(KEGG_plot_hyper@result,n=10),head(KEGG_plot_hypo@result,n=10)))
KEGG_plot$group<-c(rep("hyper",10),rep("hypo",10))
#KEGG_plot$group<-c(rep("hyper",10),rep("hypo",9))
str(KEGG_plot)
KEGG_plot$qvalue_log10<-c(-log(KEGG_plot$qvalue,10))
KEGG_plot$pvalue_log10<-c(-log(KEGG_plot$pvalue,10))

#BP_plot$group<- factor(x =BP_plot$group,levels =c("Trobt_Vims","EVTs","VCTs","STBs") )
KEGG_plot$group<- factor(x =KEGG_plot$group,levels =c("hyper","hypo") )
levels(KEGG_plot$group)
KEGG_plot<-KEGG_plot[,c("Description","group","qvalue_log10","Count")]
KEGG_plot$Description<-as.character(KEGG_plot$Description)
KEGG_plot<-KEGG_plot[order(KEGG_plot$group,KEGG_plot$qvalue_log10),]
str(KEGG_plot)
head(KEGG_plot)

PP_plot<-KEGG_plot
head(PP_plot)

pdf(paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_",depth,"_CpGs_",split_region,"_DMRs_15q005_hypo-hyper_KEGG_3kb.pdf"))

p <- ggplot(PP_plot,aes(qvalue_log10 ,reorder(Description,qvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(p.adjust)",y="Description",title="KEGG enrichment")+
  scale_size_continuous(range=c(5,14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=55))+
  geom_text(aes(label=sprintf("%.3f",qvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")

p


KEGG_plot<-as.data.frame(rbind(head(KEGG_plot_hyper@result,n=10),head(KEGG_plot_hypo@result,n=10)))
KEGG_plot$group<-c(rep("hyper",10),rep("hypo",10))
#KEGG_plot$group<-c(rep("hyper",10),rep("hypo",9))
str(KEGG_plot)
KEGG_plot$qvalue_log10<-c(-log(KEGG_plot$qvalue,10))
KEGG_plot$pvalue_log10<-c(-log(KEGG_plot$pvalue,10))

#BP_plot$group<- factor(x =BP_plot$group,levels =c("Trobt_Vims","EVTs","VCTs","STBs") )
KEGG_plot$group<- factor(x =KEGG_plot$group,levels =c("hyper","hypo") )
levels(KEGG_plot$group)

KEGG_plot<-KEGG_plot[,c("Description","group","pvalue_log10","Count")]
KEGG_plot$Description<-as.character(KEGG_plot$Description)
KEGG_plot<-KEGG_plot[order(KEGG_plot$group,KEGG_plot$pvalue_log10),]
str(KEGG_plot)
head(KEGG_plot)

PP_plot<-KEGG_plot
head(PP_plot)
ggplot(PP_plot,aes(pvalue_log10 ,reorder(Description,pvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(pvalue)",y="Description",title="KEGG enrichment")+
  scale_size_continuous(range=c(5,14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=55))+
  geom_text(aes(label=sprintf("%.3f",pvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")


dev.off()


##for genes unselected
#For hyper
gene_erichment_results<-gene_erichment_results_hyper_unselect
BP_plot<-gene_erichment_results[[1]]
MF_plot<-gene_erichment_results[[2]]
CC_plot<-gene_erichment_results[[3]]

#结果展示

merge_three<-as.data.frame(rbind(head(BP_plot@result,n=10),head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
dim(merge_three)
merge_three$qvalue_log10<-c(-log(merge_three$qvalue,10))
merge_three$pvalue_log10<-c(-log(merge_three$pvalue,10))

merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))

pdf(paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hyper_",depth,"_CpGs_",split_region,"_DMRs_15q005_enrichment_BP_CC_MF_unselect.pdf"))
ggbarplot(merge_three, x="Description", y="qvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          #          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=90)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55)) #按组排序 x.text.angle=90
#修进：http://blog.sciencenet.cn/blog-3334560-1091714.html
ggbarplot(merge_three, x="Description", y="pvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
                   rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=90)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55)) #按组排序 x.text.angle=90
#修进：http://blog.sciencenet.cn/blog-3334560-1091714.html
dev.off()

#For hypO
gene_erichment_results<-gene_erichment_results_hypo_unselct
BP_plot<-gene_erichment_results[[2]]
MF_plot<-gene_erichment_results[[3]]
CC_plot<-gene_erichment_results[[1]]

#结果展示
merge_three<-as.data.frame(rbind(head(BP_plot@result,n=10),head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
dim(merge_three)
#merge_three$qlogtran<-c(-log(merge_three$qvalue,10))
#merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))
#merge_two<-as.data.frame(rbind(head(CC_plot@result,n=10),head(MF_plot@result,n=10)))
merge_three$qvalue_log10<-c(-log(merge_three$qvalue,10))
merge_three$pvalue_log10<-c(-log(merge_three$pvalue,10))
merge_three$group<-c(rep("BP",10),rep("CC",10),rep("MF",10))

pdf(paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_hypo_",depth,"_CpGs_",split_region,"_DMRs_15q005_enrichment_CC_MF_unselect.pdf"))
ggbarplot(merge_three, x="Description", y="qvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
          #          rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=90)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55)) #按组排序 x.text.angle=90
#修进：http://blog.sciencenet.cn/blog-3334560-1091714.html
ggbarplot(merge_three, x="Description", y="pvalue_log10", fill = "group", color = "white", 
          palette = "aaas", #杂志Science的配色 "npg"
          #          sort.val = "asc", #上升排序,区别于desc，具体看图演示 "desc", #下降排序
          sort.val = "desc", 
                    rotate=TRUE,
          sort.by.groups=TRUE,x.text.angle=60)+
  scale_x_discrete(labels=function(x) str_wrap(x, width=55))

ggplot(merge_three,aes(pvalue_log10 ,reorder(Description,pvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[3:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(pvalue)",y="Description",title="GO enrichment")+
  scale_size_continuous(range=c(1,10))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))+
  geom_text(aes(label=sprintf("%.3f",qvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")

dev.off()

KEGG_plot_hyper<-gene_erichment_results_hyper_unselect[[7]] 
KEGG_plot_hypo<-gene_erichment_results_hypo_unselct[[7]] 

dim(KEGG_plot_hyper@result)
dim(KEGG_plot_hypo@result)

KEGG_plot<-as.data.frame(rbind(head(KEGG_plot_hyper@result,n=10),head(KEGG_plot_hypo@result,n=10)))
KEGG_plot$group<-c(rep("hyper",10),rep("hypo",10))
str(KEGG_plot)
KEGG_plot$qvalue_log10<-c(-log(KEGG_plot$qvalue,10))
KEGG_plot$pvalue_log10<-c(-log(KEGG_plot$pvalue,10))

#BP_plot$group<- factor(x =BP_plot$group,levels =c("Trobt_Vims","EVTs","VCTs","STBs") )
KEGG_plot$group<- factor(x =KEGG_plot$group,levels =c("hyper","hypo") )
levels(KEGG_plot$group)
KEGG_plot<-KEGG_plot[,c("Description","group","qvalue_log10","Count")]
KEGG_plot$Description<-as.character(KEGG_plot$Description)
KEGG_plot<-KEGG_plot[order(KEGG_plot$group,KEGG_plot$qvalue_log10),]
str(KEGG_plot)
head(KEGG_plot)

PP_plot<-KEGG_plot
head(PP_plot)

pdf(paste0("D:/huahua/5.DMR_function_enrichment/",compare_name,"_",depth,"_CpGs_",split_region,"_DMRs_15q005_hypo-hyper_KEGG_unselect.pdf"))

p <- ggplot(PP_plot,aes(qvalue_log10 ,reorder(Description,qvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(p.adjust)",y="Description",title="KEGG enrichment")+
  scale_size_continuous(range=c(5,14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=55))+
  geom_text(aes(label=sprintf("%.3f",qvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")

p


KEGG_plot<-as.data.frame(rbind(head(KEGG_plot_hyper@result,n=10),head(KEGG_plot_hypo@result,n=10)))
KEGG_plot$group<-c(rep("hyper",10),rep("hypo",10))
str(KEGG_plot)
KEGG_plot$qvalue_log10<-c(-log(KEGG_plot$qvalue,10))
KEGG_plot$pvalue_log10<-c(-log(KEGG_plot$pvalue,10))

#BP_plot$group<- factor(x =BP_plot$group,levels =c("Trobt_Vims","EVTs","VCTs","STBs") )
KEGG_plot$group<- factor(x =KEGG_plot$group,levels =c("hyper","hypo") )
levels(KEGG_plot$group)

KEGG_plot<-KEGG_plot[,c("Description","group","pvalue_log10","Count")]
KEGG_plot$Description<-as.character(KEGG_plot$Description)
KEGG_plot<-KEGG_plot[order(KEGG_plot$group,KEGG_plot$pvalue_log10),]
str(KEGG_plot)
head(KEGG_plot)

PP_plot<-KEGG_plot
head(PP_plot)
ggplot(PP_plot,aes(pvalue_log10 ,reorder(Description,pvalue_log10)))+#让纵轴的Description的显示顺序按GeneRatio_num值排序
  geom_point(aes(size=Count,color=group),alpha =0.5)+# 修改点的大小
  #      scale_fill_manual(values=ppCor[4:1])+
  scale_color_manual(values=ppCor[4:1])+
  #  scale_color_brewer(values = ppCor[4:1])+
  labs(color="Cell types",size="Count",x="-log10(pvalue)",y="Description",title="KEGG enrichment")+
  scale_size_continuous(range=c(5,14))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=55))+
  geom_text(aes(label=sprintf("%.3f",pvalue_log10)), size=4,vjust = 0, nudge_y = 0.25)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title = element_text(size=15),axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=15),legend.position = "right")
dev.off()
