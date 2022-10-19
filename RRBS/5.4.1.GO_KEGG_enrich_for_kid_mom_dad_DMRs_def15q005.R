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

#perform enrichment analysis
for ( fname in c("father_AMA_vs_Young","mother_AMA_vs_Young","kids_AMA_vs_Young")){
#compare_name<-"father_AMA_vs_Young"
compare_name<-fname
print(compare_name)
split_region<- "200bp"
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
hyper_gene4<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
hypo_gene4<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))

##make gene list 
gene <- list(Hyper=hyper_gene4,Hypo=hypo_gene4)
##enrichment analysis for special list
#for hyper
gene_list1<- bitr(unique(gene$Hyper), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) 
target<-gene_list1$ENTREZID
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

write.csv(as.data.frame(CC@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hyper_CpGs_",split_region,"_15q005_unselect_CC_GO.csv"))
write.csv(as.data.frame(BP@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hyper_CpGs_",split_region,"_15q005_unselect_BP_GO.csv")) 
write.csv(as.data.frame(MF@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hyper_CpGs_",split_region,"_15q005_unselect_MF_GO.csv")) 
write.csv(as.data.frame(CC_simp@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hyper_CpGs_",split_region,"_15q005_unselect_CC_GO_simp.csv")) 
write.csv(as.data.frame(BP_simp@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hyper_CpGs_",split_region,"_15q005_unselect_BP_GO_simp.csv")) 
write.csv(as.data.frame(MF_simp@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hyper_CpGs_",split_region,"_15q005_unselect_MF_GO_simp.csv")) 
write.csv((as.data.frame(kk@result)), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hyper_CpGs_",split_region,"_15q005_unselect_KEGG.csv")) 
write.csv((as.data.frame(rpea@result)), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hyper_CpGs_",split_region,"_15q005_unselect_RPEA.csv")) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hyper_CpGs_",split_region,"_15q005_unselect_GO_RPEA_result.rds"))

#for hypo
gene_list1<- bitr(unique(gene$Hypo), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) 
target<-gene_list1$ENTREZID
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

write.csv(as.data.frame(CC@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hypo_CpGs_",split_region,"_15q005_unselect_CC_GO.csv"))
write.csv(as.data.frame(BP@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hypo_CpGs_",split_region,"_15q005_unselect_BP_GO.csv")) 
write.csv(as.data.frame(MF@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hypo_CpGs_",split_region,"_15q005_unselect_MF_GO.csv")) 
write.csv(as.data.frame(CC_simp@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hypo_CpGs_",split_region,"_15q005_unselect_CC_GO_simp.csv")) 
write.csv(as.data.frame(BP_simp@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hypo_CpGs_",split_region,"_15q005_unselect_BP_GO_simp.csv")) 
write.csv(as.data.frame(MF_simp@result), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hypo_CpGs_",split_region,"_15q005_unselect_MF_GO_simp.csv")) 
write.csv((as.data.frame(kk@result)), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hypo_CpGs_",split_region,"_15q005_unselect_KEGG.csv")) 
write.csv((as.data.frame(rpea@result)), file=paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hypo_CpGs_",split_region,"_15q005_unselect_RPEA.csv")) 
gene_GO_RPEA_erichment_results=list(CC,BP,MF,CC_simp,BP_simp,MF_simp,kk,rpea)
saveRDS(gene_GO_RPEA_erichment_results, file = paste0("/mnt/data/chenwei/huahua/manuscript/figure/figure2/enrichment_for_DMRS/",compare_name,"_hypo_CpGs_",split_region,"_15q005_unselect_GO_RPEA_result.rds"))

}
