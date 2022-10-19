rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(VennDiagram)
library(scales)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(UpSetR)
#step0 :set colors
pal <- pal_npg("nrc", alpha=1)(9)
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9)])
Cells_col_raw<-c("#A20056FF","#EE4C97FF","#F39B7FFF","#E18727FF","#FFDC91FF","#00A1D5FF", "#0072B5FF","#20854EFF","#631879FF","#7E6148FF","#DC0000FF")
##extend colors
pal1<-pal_nejm("default",alpha = 1)(8)##(8表示呈现多少个颜色)nejm，共8种
pal2<-pal_jama("default",alpha = 1)(7)##(8表示呈现多少个颜色)nejm，共8种
pal3<- pal_aaas("default",alpha=1)(10)
pal4 <- pal_npg("nrc", alpha=1)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
pal5 <- pal_npg("nrc", alpha=0.5)(10)#nrc是Palette Types，alpha用于调节透明度，共10种
ppCor_all <-c(pal1,pal2,pal3,pal4,pal5)
show_col(ppCor)
#read expression data 

#甲基化对应基因列表的读取
#FOR kids
#read meth data 
compare_name<-"kids_AMA_vs_Young";split_region<- "200bp";depth<-"6X"
#data preparation
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
nrow(DMR_hyper);nrow(DMR_hypo)#182 269
hyper_gene<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
hypo_gene<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))
dim(data_region_merge_DMR_all_big);length(hyper_gene);length(hypo_gene)
#451 #170 #236

#load GenAge database genes
longevity_gene0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/longevity.csv")
longevity_gene<-unique(as.character(unlist(strsplit(as.character(longevity_gene0$Gene.s.),","))))
length(longevity_gene)
Age_gene_lgh0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/human_Age_gene_new_500.csv")
Age_gene_lgh<-unique(as.character(Age_gene_lgh0$Symbol))

genAge_database_extand<-unique(c(Age_gene_lgh,longevity_gene));length(genAge_database_extand) #1223

length(Reduce(intersect,list(hyper_gene,Age_gene_lgh)))#5
length(Reduce(intersect,list(hypo_gene,Age_gene_lgh)))#6

Reduce(intersect,list(hyper_gene,Age_gene_lgh))#"PRDX1" "SOCS2" "RB1"   "TXN"   "SYK"
Reduce(intersect,list(hypo_gene,Age_gene_lgh))#"RPS6KA5" "TRPV1"   "CACNA1A" "PTK2"    "TERF1"   "PSAT1"  
Reduce(intersect,list(hyper_gene,longevity_gene))#"TXNRD1"  "SOCS2"   "CELSR1"  "ST6GAL1" "SYNE1"   "CSMD3"   "DEPTOR"  "SYK" 
Reduce(intersect,list(hypo_gene,longevity_gene))#"LAPTM5","ADARB2","LINC01234","ITGB3","POLRMT","CELSR1","PRKN","TERF1","LMX1B","COL5A1

grid.newpage(); 
venn_Age_longevity <-venn.diagram(list(Age_gene=Age_gene_lgh,longevity_gene=longevity_gene,AMA_DMR_genes_all=unique(c(hyper_gene,hypo_gene))),
                                  alpha=c(0.9,0.9,0.9),lwd=1,lty=1,col="black" , 
                                  fill=ppCor_all[c(3,6,13)],cex = 1.5,cat.col=ppCor_all[c(3,6,13)],cat.fontface=4, cat.cex = 1.5,    
                                  main = "Database:: longevity gene and two list of Age gene",
                                  main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                  filename = NULL)
grid.draw(venn_Age_longevity)
grid.newpage();
venn_Age_DEGs_genAge<-venn.diagram(list(genAge_database=genAge_database_extand,DMR_hyper_gene=hyper_gene,DMR_hypo_gene=hypo_gene),
                                   alpha=c(1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(6,1:2)], cex = 1.5, 
                                   cat.col=ppCor[c(6,1:2)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)
grid.newpage(); #清空画板，开始画新图

pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_annotation_genAge_longevity_gene_kids_DMRs.pdf",width = 8,height = 8)
grid.draw(venn_Age_longevity)
grid.newpage();
grid.draw(venn_Age_DEGs_genAge)
dev.off()

Reduce(intersect,list(genAge_database_extand,hyper_gene))
#"PRDX1"   "RB1"     "SOCS2"   "SYK"     "TXN"     "TXNRD1"  "CELSR1"  "DEPTOR"  "SYNE1"   "CSMD3"   "ST6GAL1"
Reduce(intersect,list(genAge_database_extand,hypo_gene))
#[1] "PSAT1"     "CACNA1A"   "PTK2"      "RPS6KA5"   "TERF1"     "TRPV1"     "ITGB3"     "ADARB2"    "LAPTM5"   
#[10] "CELSR1"    "COL5A1"    "LINC01234" "PRKN"      "LMX1B"     "POLRMT" 
Reduce(intersect,list(genAge_database_extand,hypo_gene,hyper_gene))# "CELSR1"

## SASP marker genes in DEGs
SASP_gene0 <- read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/SASP_human.csv")
SASP_gene<-unique(as.character(SASP_gene0$V1))
Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),SASP_gene))#"RPS6KA5"
Reduce(intersect,list(hypo_gene,SASP_gene))# "RPS6KA5" 
Reduce(intersect,list(hyper_gene,SASP_gene)) 
head(data_region_merge_DMR_all_big)
data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$SYMBOL ==  "RPS6KA5"),]
data_region_merge_DMR_all_big[which(data_region_merge_DMR_all_big$SYMBOL %in% Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),SASP_gene))),c("Row.names","meth.diff","annotation","SYMBOL")]
# Row.names meth.diff                                      annotation  SYMBOL
# chr14_91033401_91033600 -23.22648 Intron (ENST00000614987.5/9252, intron 1 of 16) RPS6KA5

venn_Age_DMRs_SASP_gene<-venn.diagram(list(SASP_gene=SASP_gene,AMA_DMRs_genes_all=unique(c(hyper_gene,hypo_gene))),
                                   alpha=c(1,1),lwd=1,lty=1,  col="black" , fill=ppCor[c(10,6)], 
                                   cex = 1.5, cat.col=ppCor[c(10,6)],cat.fontface=4,cat.cex = 1.5,
                                   main.cex = 2, main.fontface = 2, main.fontfamily = 3,  filename = NULL)
grid.newpage()
grid.draw(venn_Age_DMRs_SASP_gene)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_annotation_SASP_gene_kids_DMRs.pdf",width = 8,height = 8)
grid.draw(venn_Age_DMRs_SASP_gene)
dev.off()


#read imprinting genes ref:"Genomic Imprinting and Physiological Processes in Mammals"
imprinting_gene0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/human_imprint_gene_single.csv",header=T)
head(imprinting_gene0)
imprinting_genes<-unique(as.character(imprinting_gene0$Gene))
length(imprinting_genes)#228

Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),imprinting_genes))
#"RB1"     "CACNA1A"
grid.newpage(); #清空画板，开始画新图
venn_Age_imprinted_genes1 <-venn.diagram(list(imprinted_genes=imprinting_genes,AMA_DMRs_genes_all=unique(c(hyper_gene,hypo_gene))),
                          alpha=c(0.9,0.9),lwd=1,lty=1,col="black" , 
                          fill=ppCor[3:4],cex = 1.5,cat.col=ppCor[3:4],cat.fontface=4, cat.cex = 1.5,    
                          main = "imprinted_genes and AMA DMRs related genes",
                          main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                          filename = NULL)
grid.draw(venn_Age_imprinted_genes1)
venn_Age_imprinted_genes2<-venn.diagram(list(imprinted_genes=imprinting_genes,DMR_hyper_gene=hyper_gene,DMR_hypo_gene=hypo_gene),
                         alpha=c(1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(6,1:2)], cex = 1.5, 
                         cat.col=ppCor[c(6,1:2)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                         filename = NULL)
grid.newpage();
grid.draw(venn_Age_imprinted_genes2)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_Age_imprinted_genes_kids_DMRs.pdf",width = 8,height = 8)
grid.draw(venn_Age_imprinted_genes1)
grid.newpage()
grid.draw(venn_Age_imprinted_genes2)
dev.off()

##human TFs and TFs_cofactor 
Homo_sapiens_TF <- read.table("/mnt/data/chenwei/qinmen_BR/00.ref_data/Homo_sapiens_TF.txt", sep = "\t", header = T)
Homo_sapiens_coTF <- read.table("/mnt/data/chenwei/qinmen_BR/00.ref_data/Homo_sapiens_TF_cofactors.txt", sep = "\t", header = T)
Hm_TFs<-unique(as.character(Homo_sapiens_TF$Symbol))
Hm_coTFs<-unique(as.character(Homo_sapiens_coTF$Symbol))

Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),Hm_TFs))
#  [1] "FOXN4"  "ZNF536" "ZNF565" "ZNF229" "CENPA"  "TOX2"   "LHX6"   "HIVEP3" "NPAS3"  "FOXL1"  "TGIF1"  "ZNF407"
#  [13] "ZNF468" "TSHZ2"  "NKX1-1" "POU3F2" "GLI3"   "TERF1"  "LMX1B"  
Reduce(intersect,list(hyper_gene,Hm_TFs))
#"FOXN4"  "ZNF536" "ZNF565" "ZNF229" "CENPA"  "TOX2"   "LHX6"  
Reduce(intersect,list(hypo_gene,Hm_TFs))
#"HIVEP3" "NPAS3"  "FOXL1"  "TGIF1"  "ZNF407" "ZNF468" "TSHZ2"  "NKX1-1" "POU3F2" "GLI3"   "TERF1"  "LMX1B" 

Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),Hm_coTFs))
#  [1] "PEX14"  "PINK1"  "MEAF6"  "SUFU"   "RB1"    "BCAS3"  "MED16"  "ZNF638" "VGLL2"  "TXN"    "TAF1L"  "SYK"   
# [13] "CENPF"  "CITED4" "CTDP1"  "RNF168" "HDAC9"
Reduce(intersect,list(hyper_gene,Hm_coTFs))
# [1] "PEX14"  "PINK1"  "MEAF6"  "SUFU"   "RB1"    "BCAS3"  "MED16"  "ZNF638" "VGLL2"  "TXN"    "TAF1L"  "SYK"   
Reduce(intersect,list(hypo_gene,Hm_coTFs))
# "CENPF"  "CITED4" "CTDP1"  "RNF168" "HDAC9" 

venn_Age_DEGs_Hg_TF_co_TFs<-venn.diagram(list(DMR_hyper_gene=hyper_gene,DMR_hypo_gene=hypo_gene,Homo_sapiens_TF=Hm_TFs,Homo_sapiens_coTF=Hm_coTFs),
                                   alpha=c(1,1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1:3,5)],cex = 1.5, cat.col=ppCor[c(1:3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                   cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_Hg_TF_co_TFs_kids_DMRs.pdf",width = 8,height = 8)
grid.newpage()
grid.draw(venn_Age_DEGs_Hg_TF_co_TFs)
dev.off()



#甲基化对应基因列表的读取
#FOR mom
#read meth data 
compare_name<-"mother_AMA_vs_Young";split_region<- "200bp";depth<-"6X"
#data preparation
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
nrow(DMR_hyper);nrow(DMR_hypo)#236 280
hyper_gene<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
hypo_gene<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))
dim(data_region_merge_DMR_all_big);length(hyper_gene);length(hypo_gene)
#516 # 212 #253

#load GenAge database genes
longevity_gene0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/longevity.csv")
longevity_gene<-unique(as.character(unlist(strsplit(as.character(longevity_gene0$Gene.s.),","))))
length(longevity_gene)#884
Age_gene_lgh0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/human_Age_gene_new_500.csv")
Age_gene_lgh<-unique(as.character(Age_gene_lgh0$Symbol))

genAge_database_extand<-unique(c(Age_gene_lgh,longevity_gene));length(genAge_database_extand) #1223

length(Reduce(intersect,list(hyper_gene,Age_gene_lgh)))#2
length(Reduce(intersect,list(hypo_gene,Age_gene_lgh)))#11

Reduce(intersect,list(hyper_gene,Age_gene_lgh))#"SOCS2" "SYK" 
Reduce(intersect,list(hypo_gene,Age_gene_lgh))#"EDARADD" "RPS6KA5" "IGF1R"   "ADCY9"   "TRPV1"   "PLCG1"   "IL2"     "CAMK4"   "YWHAZ"   "EHMT1"   "PSAT1"
Reduce(intersect,list(hyper_gene,longevity_gene))#[1] "AGT"      "ADARB2"   "GRAMD1B"  "NAV2"     "INS-IGF2" "NADSYN1"  "TENM4"    "SOCS2"    "HPCAL1"   "LYPD6"   
#[11] "CELSR1"   "MECOM"    "SYNE1"    "ADAMTSL1" "DOCK8"    "SYK" 
Reduce(intersect,list(hypo_gene,longevity_gene))#[1] "CR1","CCND1","CERS3-AS1","ADAMTS7","IGF1R","IL2","CAMK4","MAPKAP1","LMX1B","COL5A1"   
grid.newpage(); 
venn_Age_longevity <-venn.diagram(list(Age_gene=Age_gene_lgh,longevity_gene=longevity_gene,AMA_DMR_genes_all=unique(c(hyper_gene,hypo_gene))),
                                  alpha=c(0.9,0.9,0.9),lwd=1,lty=1,col="black" , 
                                  fill=ppCor_all[c(3,6,13)],cex = 1.5,cat.col=ppCor_all[c(3,6,13)],cat.fontface=4, cat.cex = 1.5,    
                                  main = "Database:: longevity gene and two list of Age gene",
                                  main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                  filename = NULL)
grid.draw(venn_Age_longevity)
grid.newpage();
venn_Age_DEGs_genAge<-venn.diagram(list(genAge_database=genAge_database_extand,DMR_hyper_gene=hyper_gene,DMR_hypo_gene=hypo_gene),
                                   alpha=c(1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(6,1:2)], cex = 1.5, 
                                   cat.col=ppCor[c(6,1:2)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)

pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_annotation_genAge_longevity_gene_mom_DMRs.pdf",width = 8,height = 8)
grid.draw(venn_Age_longevity)
grid.newpage();
grid.draw(venn_Age_DEGs_genAge)
dev.off()

Reduce(intersect,list(genAge_database_extand,hyper_gene))
#[1] "SOCS2"    "SYK"      "AGT"      "HPCAL1"   "ADARB2"   "LYPD6"    "CELSR1"   "NADSYN1"  "GRAMD1B"  "MECOM"   
#[11] "SYNE1"    "TENM4"    "ADAMTSL1" "DOCK8"    "NAV2"     "INS-IGF2"
Reduce(intersect,list(genAge_database_extand,hypo_gene))
#[1] "PSAT1"     "ADCY9"     "CAMK4"     "EDARADD"   "EHMT1"     "IGF1R"     "IL2"       "PLCG1"     "RPS6KA5"  
#[10] "TRPV1"     "YWHAZ"     "CERS3-AS1" "MAPKAP1"   "COL5A1"    "LMX1B"     "ADAMTS7"   "CCND1"     "CR1"     
Reduce(intersect,list(genAge_database_extand,hypo_gene,hyper_gene))

## SASP marker genes in DEGs
SASP_gene0 <- read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/SASP_human.csv")
SASP_gene<-unique(as.character(SASP_gene0$V1))
Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),SASP_gene))#"ETS2"    "RPS6KA5"
Reduce(intersect,list(hypo_gene,SASP_gene))# "RPS6KA5"
Reduce(intersect,list(hyper_gene,SASP_gene)) #ETS2 
data_region_merge_DMR_all_big[which(mom_data_region_merge_DMR_all_big$SYMBOL %in% Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),SASP_gene))),c("Row.names","meth.diff","annotation","SYMBOL")]
#        Row.names         meth.diff                                      annotation  SYMBOL
# chr14_91033401_91033600 -16.13940 Intron (ENST00000614987.5/9252, intron 1 of 16) RPS6KA5
# chr21_38787401_38787600  19.96044                               Distal Intergenic    ETS2

venn_Age_DMRs_SASP_gene<-venn.diagram(list(SASP_gene=SASP_gene,AMA_DMRs_genes_all=unique(c(hyper_gene,hypo_gene))),
                                      alpha=c(1,1),lwd=1,lty=1,  col="black" , fill=ppCor[c(10,6)], 
                                      cex = 1.5, cat.col=ppCor[c(10,6)],cat.fontface=4,cat.cex = 1.5,
                                      main.cex = 2, main.fontface = 2, main.fontfamily = 3,  filename = NULL)
grid.newpage()
grid.draw(venn_Age_DMRs_SASP_gene)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_annotation_SASP_gene_mom_DMRs.pdf",width = 8,height = 8)
grid.newpage()
grid.draw(venn_Age_DMRs_SASP_gene)
dev.off()


#read imprinting genes ref:"Genomic Imprinting and Physiological Processes in Mammals"
imprinting_gene0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/human_imprint_gene_single.csv",header=T)
head(imprinting_gene0)
imprinting_genes<-unique(as.character(imprinting_gene0$Gene))
length(imprinting_genes)#228

Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),imprinting_genes))#"NAV2"   "DLGAP2" "KCNK9" 
grid.newpage(); #清空画板，开始画新图
venn_Age_imprinted_genes1 <-venn.diagram(list(imprinted_genes=imprinting_genes,AMA_DMRs_genes_all=unique(c(hyper_gene,hypo_gene))),
                                         alpha=c(0.9,0.9),lwd=1,lty=1,col="black" , 
                                         fill=ppCor[3:4],cex = 1.5,cat.col=ppCor[3:4],cat.fontface=4, cat.cex = 1.5,    
                                         main = "imprinted_genes and AMA DMRs related genes",
                                         main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                         filename = NULL)
grid.draw(venn_Age_imprinted_genes1)
venn_Age_imprinted_genes2<-venn.diagram(list(imprinted_genes=imprinting_genes,DMR_hyper_gene=hyper_gene,DMR_hypo_gene=hypo_gene),
                                        alpha=c(1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(6,1:2)], cex = 1.5, 
                                        cat.col=ppCor[c(6,1:2)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                        filename = NULL)
grid.newpage();
grid.draw(venn_Age_imprinted_genes2)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_Age_imprinted_genes_mom_DMRs.pdf",width = 8,height = 8)
grid.newpage()
grid.draw(venn_Age_imprinted_genes1)
grid.newpage()
grid.draw(venn_Age_imprinted_genes2)
dev.off()

##human TFs and TFs_cofactor 
Homo_sapiens_TF <- read.table("/mnt/data/chenwei/qinmen_BR/00.ref_data/Homo_sapiens_TF.txt", sep = "\t", header = T)
Homo_sapiens_coTF <- read.table("/mnt/data/chenwei/qinmen_BR/00.ref_data/Homo_sapiens_TF_cofactors.txt", sep = "\t", header = T)
Hm_TFs<-unique(as.character(Homo_sapiens_TF$Symbol))
Hm_coTFs<-unique(as.character(Homo_sapiens_coTF$Symbol))

Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),Hm_TFs))
#  [1] "LMX1A"    "RCOR2"    "FOXN3"    "ZNF500"   "CENPT"    "ETV4"     "ZNF726"   "ZNF229"   "ZNF765"   "GLI2"    
#  [11] "SMARCAL1" "TOX2"     "ETS2"     "TSC22D2"  "MECOM"    "PRDM1"    "ZNF398"   "NEUROD6"  "LHX6"     "ZNF362"  
#  [21] "HIVEP3"   "SIX4"     "BCL6B"    "LRRFIP1"  "TFAP2C"   "SMAD1"    "NKX1-1"   "ZNF696"   "ZNF596"   "ZFP37"   
#  [31] "LMX1B"  
Reduce(intersect,list(hyper_gene,Hm_TFs))
#[1] "LMX1A"    "RCOR2"    "FOXN3"    "ZNF500"   "CENPT"    "ETV4"     "ZNF726"   "ZNF229"   "ZNF765"   "GLI2"    
#[11] "SMARCAL1" "TOX2"     "ETS2"     "TSC22D2"  "MECOM"    "PRDM1"    "ZNF398"   "NEUROD6"  "LHX6
Reduce(intersect,list(hypo_gene,Hm_TFs))
#[1] "ZNF362"  "HIVEP3"  "SIX4"    "BCL6B"   "LRRFIP1" "TFAP2C"  "SMAD1"   "NKX1-1"  "ZNF696"  "ZNF596"  "ZFP37"  
#[12] "LMX1B"
Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),Hm_coTFs))
#[1] "CHTOP"   "SUFU"    "NOTCH3"  "SLC19A1" "ELP1"    "TAF1L"   "SYK"     "CITED4"  "TCERG1L" "CCND1"   "DPF3"   
# [12] "RASD1"   "HDAC4"   "CHD6"    "LIN54"   "CAMK4"   "EHMT1"   "KDM4C"  
Reduce(intersect,list(hyper_gene,Hm_coTFs))
# [1] "CHTOP"   "SUFU"    "NOTCH3"  "SLC19A1" "ELP1"    "TAF1L"   "SYK"     
Reduce(intersect,list(hypo_gene,Hm_coTFs))
# [1] "CITED4"  "TCERG1L" "CCND1"   "DPF3"    "RASD1"   "HDAC4"   "CHD6"    "LIN54"   "CAMK4"   "EHMT1"   "KDM4C"  

venn_Age_DEGs_Hg_TF_co_TFs<-venn.diagram(list(DMR_hyper_gene=hyper_gene,DMR_hypo_gene=hypo_gene,Homo_sapiens_TF=Hm_TFs,Homo_sapiens_coTF=Hm_coTFs),
                                         alpha=c(1,1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1:3,5)],cex = 1.5, cat.col=ppCor[c(1:3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                         cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                         filename = NULL)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_Hg_TF_co_TFs_mom_DMRs.pdf",width = 8,height = 8)
grid.newpage()
grid.draw(venn_Age_DEGs_Hg_TF_co_TFs)
dev.off()



#甲基化对应基因列表的读取
#FOR dad
#read meth data 
compare_name<-"father_AMA_vs_Young";split_region<- "200bp";depth<-"6X"
#data preparation
data_region_merge_DMR_all_big<-read.table(paste0("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/",compare_name,"_",split_region,"_DMR_myDiff15q005_merge_data.txt"),header=T,sep="\t")
DMR_hyper<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff >0,]
DMR_hypo<-data_region_merge_DMR_all_big[data_region_merge_DMR_all_big$meth.diff <0,]
nrow(DMR_hyper);nrow(DMR_hypo)#206 258
hyper_gene<-na.omit(unique(as.character(DMR_hyper$SYMBOL)))
hypo_gene<-na.omit(unique(as.character(DMR_hypo$SYMBOL)))
dim(data_region_merge_DMR_all_big);length(hyper_gene);length(hypo_gene)
#464 #185 #238

#load GenAge database genes
longevity_gene0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/longevity.csv")
longevity_gene<-unique(as.character(unlist(strsplit(as.character(longevity_gene0$Gene.s.),","))))
length(longevity_gene)#884
Age_gene_lgh0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/human_Age_gene_new_500.csv")
Age_gene_lgh<-unique(as.character(Age_gene_lgh0$Symbol))

genAge_database_extand<-unique(c(Age_gene_lgh,longevity_gene));length(genAge_database_extand) #1223

length(Reduce(intersect,list(hyper_gene,Age_gene_lgh)))#2
length(Reduce(intersect,list(hypo_gene,Age_gene_lgh)))#2

Reduce(intersect,list(hyper_gene,Age_gene_lgh))#"PRDX1"  "NFKBIA"
Reduce(intersect,list(hypo_gene,Age_gene_lgh))#"TOP1"   "POU1F1"
Reduce(intersect,list(hyper_gene,longevity_gene))#[1] "ADARB2"    "LINC01234" "MDGA2"     "ST6GAL1"   "PROM1"     "MAPKAP1"
Reduce(intersect,list(hypo_gene,longevity_gene))# [1] "ST3GAL3" "CCND1"   "CACNA1C" "FARP1"   "CERS3"   "UQCRFS1" "ERBB4"   "FHIT"    "PROM1"   "PRKN"    "LMX1B"  

venn_Age_longevity <-venn.diagram(list(Age_gene=Age_gene_lgh,longevity_gene=longevity_gene,AMA_DMR_genes_all=unique(c(hyper_gene,hypo_gene))),
                                  alpha=c(0.9,0.9,0.9),lwd=1,lty=1,col="black" , 
                                  fill=ppCor_all[c(3,6,13)],cex = 1.5,cat.col=ppCor_all[c(3,6,13)],cat.fontface=4, cat.cex = 1.5,    
                                  main = "Database:: longevity gene and two list of Age gene",
                                  main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                  filename = NULL)
grid.draw(venn_Age_longevity)
grid.newpage();
venn_Age_DEGs_genAge<-venn.diagram(list(genAge_database=genAge_database_extand,DMR_hyper_gene=hyper_gene,DMR_hypo_gene=hypo_gene),
                                   alpha=c(1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(6,1:2)], cex = 1.5, 
                                   cat.col=ppCor[c(6,1:2)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                   filename = NULL)
grid.draw(venn_Age_DEGs_genAge)

pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_annotation_genAge_longevity_gene_dad_DMRs.pdf",width = 8,height = 8)
grid.newpage();
grid.draw(venn_Age_longevity)
grid.newpage();
grid.draw(venn_Age_DEGs_genAge)
dev.off()

Reduce(intersect,list(genAge_database_extand,hyper_gene))
#[1] "NFKBIA"    "PRDX1"     "ADARB2"    "MAPKAP1"   "PROM1"     "LINC01234" "MDGA2"     "ST6GAL1"  

Reduce(intersect,list(genAge_database_extand,hypo_gene))
# [1] "POU1F1"  "TOP1"    "CACNA1C" "CERS3"   "ST3GAL3" "PROM1"   "PRKN"    "LMX1B"   "FHIT"    "FARP1"   "CCND1"  
#[12] "ERBB4"   "UQCRFS1" 
Reduce(intersect,list(genAge_database_extand,hypo_gene,hyper_gene))#"PROM1"

## SASP marker genes in DEGs
SASP_gene0 <- read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/SASP_human.csv")
SASP_gene<-unique(as.character(SASP_gene0$V1))
Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),SASP_gene))
Reduce(intersect,list(hypo_gene,SASP_gene))
Reduce(intersect,list(hyper_gene,SASP_gene)) 

venn_Age_DMRs_SASP_gene<-venn.diagram(list(SASP_gene=SASP_gene,AMA_DMRs_genes_all=unique(c(hyper_gene,hypo_gene))),
                                      alpha=c(1,1),lwd=1,lty=1,  col="black" , fill=ppCor[c(10,6)], 
                                      cex = 1.5, cat.col=ppCor[c(10,6)],cat.fontface=4,cat.cex = 1.5,
                                      main.cex = 2, main.fontface = 2, main.fontfamily = 3,  filename = NULL)
grid.newpage()
grid.draw(venn_Age_DMRs_SASP_gene)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_annotation_SASP_gene_dad_DMRs.pdf",width = 8,height = 8)
grid.newpage()
grid.draw(venn_Age_DMRs_SASP_gene)
dev.off()

#read imprinting genes ref:"Genomic Imprinting and Physiological Processes in Mammals"
imprinting_gene0<-read.csv("/mnt/data/chenwei/qinmen_BR/00.ref_data/human_imprint_gene_single.csv",header=T)
head(imprinting_gene0)
imprinting_genes<-unique(as.character(imprinting_gene0$Gene))
length(imprinting_genes)#228

Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),imprinting_genes))
#[1] "RASGRF1" "NUDT12"  "TRAPPC9" "OPCML"   "CACNA1C" "CYP24A1" "PLG"     "RAPGEF5" "DENND3" 
grid.newpage();
venn_Age_imprinted_genes1 <-venn.diagram(list(imprinted_genes=imprinting_genes,AMA_DMRs_genes_all=unique(c(hyper_gene,hypo_gene))),
                                         alpha=c(0.9,0.9),lwd=1,lty=1,col="black" , 
                                         fill=ppCor[3:4],cex = 1.5,cat.col=ppCor[3:4],cat.fontface=4, cat.cex = 1.5,    
                                         main = "imprinted_genes and AMA DMRs related genes",
                                         main.cex = 1.5, main.fontface = 1.5, main.fontfamily = 1.5, 
                                         filename = NULL)
grid.draw(venn_Age_imprinted_genes1)
venn_Age_imprinted_genes2<-venn.diagram(list(imprinted_genes=imprinting_genes,DMR_hyper_gene=hyper_gene,DMR_hypo_gene=hypo_gene),
                                        alpha=c(1,1,1),lwd=1,lty=1,col="black" ,fill=ppCor[c(6,1:2)], cex = 1.5, 
                                        cat.col=ppCor[c(6,1:2)], cat.fontface=4, cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                        filename = NULL)
grid.newpage();
grid.draw(venn_Age_imprinted_genes2)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_Age_imprinted_genes_dad_DMRs.pdf",width = 8,height = 8)
grid.newpage()
grid.draw(venn_Age_imprinted_genes1)
grid.newpage()
grid.draw(venn_Age_imprinted_genes2)
dev.off()

##human TFs and TFs_cofactor 
Homo_sapiens_TF <- read.table("/mnt/data/chenwei/qinmen_BR/00.ref_data/Homo_sapiens_TF.txt", sep = "\t", header = T)
Homo_sapiens_coTF <- read.table("/mnt/data/chenwei/qinmen_BR/00.ref_data/Homo_sapiens_TF_cofactors.txt", sep = "\t", header = T)
Hm_TFs<-unique(as.character(Homo_sapiens_TF$Symbol))
Hm_coTFs<-unique(as.character(Homo_sapiens_coTF$Symbol))

Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),Hm_TFs))
#  [1] "SMAD3"   "NFATC1"  "ZNF536"  "GLI2"    "TSC22D2" "JARID2"  "RREB1"   "KLF6"    "BARX2"   "ZFHX3"   "TGIF1"  
# [12] "RFX1"    "ZNF99"   "TSHZ2"   "MYT1"    "RUNX1"   "ZXDC"    "SMARCC1" "FOXP1"   "POU1F1"  "POU3F2"  "ZBTB10" 
# [23] "PBX3"    "LMX1B" 
Reduce(intersect,list(hyper_gene,Hm_TFs))
#[1] "SMAD3"   "NFATC1"  "ZNF536"  "GLI2"    "TSC22D2" "JARID2"  "RREB1"  
Reduce(intersect,list(hypo_gene,Hm_TFs))
#[1] "KLF6"    "BARX2"   "ZFHX3"   "TGIF1"   "RFX1"    "ZNF99"   "TSHZ2"   "MYT1"    "RUNX1"   "ZXDC"    "SMARCC1"
#[12] "FOXP1"   "POU1F1"  "POU3F2"  "ZBTB10"  "PBX3"    "LMX1B
Reduce(intersect,list(unique(c(hyper_gene,hypo_gene)),Hm_coTFs))
#[1] "PEX14"   "ATF7IP"  "NFKBIA"  "PLK1"    "RASD1"   "VGLL2"   "ELP3"    "CCND1"   "SNX6"    "RAI1"    "URI1"   
#[12] "WTIP"    "ERBB4"   "TOP1"    "CHD6"    "SLC19A1" "FHIT"    "NOTCH4"  "SMURF1" 
Reduce(intersect,list(hyper_gene,Hm_coTFs))
# [1] "PEX14"  "ATF7IP" "NFKBIA" "PLK1"   "RASD1"  "VGLL2"  "ELP3"  

Reduce(intersect,list(hypo_gene,Hm_coTFs))
#  [1] "PEX14"   "CCND1"   "SNX6"    "RAI1"    "URI1"    "WTIP"    "ERBB4"   "TOP1"    "CHD6"    "SLC19A1" "FHIT"   
# [12] "NOTCH4"  "SMURF1" 
venn_Age_DEGs_Hg_TF_co_TFs<-venn.diagram(list(DMR_hyper_gene=hyper_gene,DMR_hypo_gene=hypo_gene,Homo_sapiens_TF=Hm_TFs,Homo_sapiens_coTF=Hm_coTFs),
                                         alpha=c(1,1,1,1),lwd=1,lty=1, col="black" ,  fill=ppCor[c(1:3,5)],cex = 1.5, cat.col=ppCor[c(1:3,5)],#cat.col表示集合名称的显示颜色。 #分类颜色 
                                         cat.fontface=4,cat.cex = 1.5, main.cex = 2, main.fontface = 2, main.fontfamily = 3, 
                                         filename = NULL)
pdf("/mnt/data/chenwei/huahua/4.methy_result/1.DMR_called/venn_data_Hg_TF_co_TFs_dad_DMRs.pdf",width = 8,height = 8)
grid.newpage()
grid.draw(venn_Age_DEGs_Hg_TF_co_TFs)
dev.off()
