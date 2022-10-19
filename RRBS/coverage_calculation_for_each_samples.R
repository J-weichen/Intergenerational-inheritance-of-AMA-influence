#loading package
rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ggpubr)
library(methylKit)
library(scatterplot3d)
library(pcaMethods)
library(scales)
library(ggsci)
library(grid)
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

grid.newpage(); #清空画板，开始画新图

AMA_father_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Father")
YOUNG_father_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Father")
AMA_Mother_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Mother")
YOUNG_Mother_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Mother")
AMA_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Kids")
YOUNG_Kids_files<-list.files("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Kids")

AMA_father_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Father/",AMA_father_files))
YOUNG_father_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Father/",YOUNG_father_files))
AMA_Mother_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Mother/",AMA_Mother_files))
YOUNG_Mother_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Mother/",YOUNG_Mother_files))
AMA_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/AMA_file/Kids/",AMA_Kids_files))
YOUNG_Kids_files2<-as.list(paste0("/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Kids/",YOUNG_Kids_files))

AMA_father_names<- lapply(strsplit(AMA_father_files,"_"), function(x) x[1])
YOUNG_father_names<- lapply(strsplit(YOUNG_father_files,"_"), function(x) x[1])
AMA_Mother_names<- lapply(strsplit(AMA_Mother_files,"_"), function(x) x[1])
YOUNG_Mother_names<- lapply(strsplit(YOUNG_Mother_files,"_"), function(x) x[1])
AMA_Kids_names<- lapply(strsplit(AMA_Kids_files,"_"), function(x) x[1])
YOUNG_Kids_names<- lapply(strsplit(YOUNG_Kids_files,"_"), function(x) x[1])
#covarage 6
myobj=methRead(c(AMA_Kids_files2,YOUNG_Kids_files2,AMA_Mother_files2,YOUNG_Mother_files2,AMA_father_files2,YOUNG_father_files2),
               sample.id=c(AMA_Kids_names,YOUNG_Kids_names,AMA_Mother_names,YOUNG_Mother_names,AMA_father_names,YOUNG_father_names),
               assembly="hg38",pipeline="bismarkCoverage",#default:'amp'
               treatment=c(rep(0,length(AMA_Kids_names)),rep(1,length(YOUNG_Kids_names)),
                           rep(2,length(AMA_Mother_names)),rep(3,length(YOUNG_Mother_names)),
                           rep(4,length(AMA_father_names)),rep(5,length(YOUNG_father_names))),
               context="CpG", resolution = "base", mincov = 1)
dim(as.data.frame(myobj[1]))

saveRDS(myobj, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample_myobj.rds")

depth<-1
Mobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=NULL)
#remain CpGs exist in every single sample 
meth_1bp=methylKit::unite(Mobj, destrand=FALSE,min.per.group=1L,mc.cores = 48)
write.table(as.data.frame(meth_1bp), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample_covarage.txt",quote=F, row.names=F, col.names=T) 

#将样本聚合在一起
meth<-meth_1bp
head(meth)
meth@sample.ids
unique(meth$chr)
meth2<-meth[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",meth$chr,invert=TRUE),]
unique(meth2$chr)
meth2$chr<-factor(meth2$chr,levels =as.character(unique(meth2$chr)))
dim(meth2); dim(na.omit(meth2))
#[1] 25871217      184
#[1] 6277099     184
saveRDS(meth2, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample_meth2.rds")

#meth2<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample_meth2.rds")

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth2,rowids=T))
head(data_raw_bin)
saveRDS(data_raw_bin, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample_methlevel.rds")

myobj<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample_myobj.rds")



depth<-3
Mobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=NULL)
#remain CpGs exist in every single sample 
meth_1bp=methylKit::unite(Mobj, destrand=FALSE,min.per.group=1L,mc.cores = 48)
write.table(as.data.frame(meth_1bp), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_3X_coverage_bismark_1sample_covarage.txt",quote=F, row.names=F, col.names=T) 

#将样本聚合在一起
meth<-meth_1bp
head(meth)
meth@sample.ids
unique(meth$chr)
meth2<-meth[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",meth$chr,invert=TRUE),]
unique(meth2$chr)
meth2$chr<-factor(meth2$chr,levels =as.character(unique(meth2$chr)))
dim(meth2); dim(na.omit(meth2))
saveRDS(meth2, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_3X_coverage_bismark_1sample_meth2.rds")

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth2,rowids=T))
head(data_raw_bin)
saveRDS(data_raw_bin, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_3X_coverage_bismark_1sample_methlevel.rds")


depth<-5
Mobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=NULL)
#remain CpGs exist in every single sample 
meth_1bp=methylKit::unite(Mobj, destrand=FALSE,min.per.group=1L,mc.cores = 48)
write.table(as.data.frame(meth_1bp), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_5X_coverage_bismark_1sample_covarage.txt",quote=F, row.names=F, col.names=T) 

#将样本聚合在一起
meth<-meth_1bp
head(meth)
meth@sample.ids
unique(meth$chr)
meth2<-meth[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",meth$chr,invert=TRUE),]
unique(meth2$chr)
meth2$chr<-factor(meth2$chr,levels =as.character(unique(meth2$chr)))
dim(meth2); dim(na.omit(meth2))
saveRDS(meth2, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_5X_coverage_bismark_1sample_meth2.rds")

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth2,rowids=T))
head(data_raw_bin)
saveRDS(data_raw_bin, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_5X_coverage_bismark_1sample_methlevel.rds")

depth<-10
Mobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=NULL)
#remain CpGs exist in every single sample 
meth_1bp=methylKit::unite(Mobj, destrand=FALSE,min.per.group=1L,mc.cores = 48)
write.table(as.data.frame(meth_1bp), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_10X_coverage_bismark_1sample_covarage.txt",quote=F, row.names=F, col.names=T) 

#将样本聚合在一起
meth<-meth_1bp
head(meth)
meth@sample.ids
unique(meth$chr)
meth2<-meth[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",meth$chr,invert=TRUE),]
unique(meth2$chr)
meth2$chr<-factor(meth2$chr,levels =as.character(unique(meth2$chr)))
dim(meth2); dim(na.omit(meth2))
saveRDS(meth2, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_10X_coverage_bismark_1sample_meth2.rds")

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth2,rowids=T))
head(data_raw_bin)
saveRDS(data_raw_bin, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_10X_coverage_bismark_1sample_methlevel.rds")


depth<-6
Mobj=filterByCoverage(myobj,lo.count=depth,lo.perc=NULL, hi.count=NULL,hi.perc=NULL)
#remain CpGs exist in every single sample 
meth_1bp=methylKit::unite(Mobj, destrand=FALSE,min.per.group=1L,mc.cores = 48)
write.table(as.data.frame(meth_1bp), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_6X_coverage_bismark_1sample_covarage.txt",quote=F, row.names=F, col.names=T) 

#将样本聚合在一起
meth<-meth_1bp
head(meth)
meth@sample.ids
unique(meth$chr)
meth2<-meth[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",meth$chr,invert=TRUE),]
unique(meth2$chr)
meth2$chr<-factor(meth2$chr,levels =as.character(unique(meth2$chr)))
dim(meth2); dim(na.omit(meth2))
saveRDS(meth2, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_6X_coverage_bismark_1sample_meth2.rds")

#cal methylation level
data_raw_bin<-as.data.frame(percMethylation(meth2,rowids=T))
head(data_raw_bin)
saveRDS(data_raw_bin, file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_6X_coverage_bismark_1sample_methlevel.rds")

#plot count 200bp bin in each samples
huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
head(huahua_meta)
#tiles_200_2[1:6,];meth_200bp_2[1:6,]
rownames(huahua_meta)<-huahua_meta$library_code
#

##reading meth value from target depth

###for 1X
depth<-"1X"
data_raw_bin<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample_methlevel.rds")

#filtering and arrange data
evaluation.table<-data_raw_bin
evaluation.table[!(is.na(evaluation.table))]<- 1
evaluation.table[is.na(evaluation.table)]<- 0
bin1bp_detected_number<-colSums(evaluation.table)
bin1bp_detected_number<-data.frame(bin1bp_detected_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)
#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=1,position=position_stack(1.02),angle=90)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_number_All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_1X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_1X)<-c("sample","coverage_1X")

###for 3X
depth<-"3X"
data_raw_bin<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_3X_coverage_bismark_1sample_methlevel.rds")

#filtering and arrange data
evaluation.table<-data_raw_bin
evaluation.table[!(is.na(evaluation.table))]<- 1
evaluation.table[is.na(evaluation.table)]<- 0
bin1bp_detected_number<-colSums(evaluation.table)
bin1bp_detected_number<-data.frame(bin1bp_detected_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)

#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)

bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=1,position=position_stack(1.02),angle=90)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_number_All_60_samples_meth_level_1bp_CpG_3X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_3X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_3X)<-c("sample","coverage_3X")

###for 5X
depth<-"5X"
data_raw_bin<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_5X_coverage_bismark_1sample_methlevel.rds")

#filtering and arrange data
evaluation.table<-data_raw_bin
evaluation.table[!(is.na(evaluation.table))]<- 1
evaluation.table[is.na(evaluation.table)]<- 0
bin1bp_detected_number<-colSums(evaluation.table)
bin1bp_detected_number<-data.frame(bin1bp_detected_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)

#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)

bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=1,position=position_stack(1.02),angle=90)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_number_All_60_samples_meth_level_1bp_CpG_5X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_5X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_5X)<-c("sample","coverage_5X")


###for 6X
depth<-"6X"
data_raw_bin<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_6X_coverage_bismark_1sample_methlevel.rds")

#filtering and arrange data
evaluation.table<-data_raw_bin
evaluation.table[!(is.na(evaluation.table))]<- 1
evaluation.table[is.na(evaluation.table)]<- 0
bin1bp_detected_number<-colSums(evaluation.table)
bin1bp_detected_number<-data.frame(bin1bp_detected_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)

#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)

bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=1,position=position_stack(1.02),angle=90)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_number_All_60_samples_meth_level_1bp_CpG_6X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_6X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_6X)<-c("sample","coverage_6X")

###for 10X
depth<-"10X"
data_raw_bin<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_10X_coverage_bismark_1sample_methlevel.rds")

#filtering and arrange data
evaluation.table<-data_raw_bin
evaluation.table[!(is.na(evaluation.table))]<- 1
evaluation.table[is.na(evaluation.table)]<- 0
bin1bp_detected_number<-colSums(evaluation.table)
bin1bp_detected_number<-data.frame(bin1bp_detected_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)

#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)

bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=1,position=position_stack(1.02),angle=90)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_number_All_60_samples_meth_level_1bp_CpG_10X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_10X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_10X)<-c("sample","coverage_10X")

CpG_number_data<-merge(CpG_number_1X,CpG_number_3X,by="sample")
CpG_number_data<-merge(CpG_number_data,CpG_number_5X,by="sample")
CpG_number_data<-merge(CpG_number_data,CpG_number_6X,by="sample")
CpG_number_data<-merge(CpG_number_data,CpG_number_10X,by="sample")
head(CpG_number_data)
write.table(as.data.frame(CpG_number_data), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_number_of_1bp_CpG_different_coverage_bismark_1sample_covarage.txt",quote=F, row.names=F, col.names=T) 

##coverage
###for 1X
depth<-"1X"
Coverage_data<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample_meth2.rds")
#Coverage_data2<-getData(Coverage_data)
#filtering and arrange data
evaluation.table<-as.data.frame(getData(Coverage_data))
head(evaluation.table)
evaluation.table$region<-paste(evaluation.table$chr,evaluation.table$start,evaluation.table$end,sep = "_")
evaluation_data<-evaluation.table[,c("region",colnames(evaluation.table)[grepl("coverage",colnames(evaluation.table))])]
rownames(evaluation_data)<-evaluation_data$region;evaluation_data<-evaluation_data[,-c(1)]
colnames(evaluation_data)<-Coverage_data@sample.ids
tail(evaluation_data)
bin1bp_coverage_number<-colMeans(evaluation_data,na.rm = T)
bin1bp_detected_number<-data.frame(bin1bp_coverage_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)
#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
bin1bp_detected_number2$bin1bp_number<-round(bin1bp_detected_number2$bin1bp_number,2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=2,position=position_stack(1.02),angle=0)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_average_coverage_depth_All_60_samples_meth_level_1bp_CpG_1X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_1X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_1X)<-c("sample","coverage_1X")

###for 3X
depth<-"3X"
Coverage_data<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_3X_coverage_bismark_1sample_meth2.rds")
#Coverage_data2<-getData(Coverage_data)
#filtering and arrange data
evaluation.table<-as.data.frame(getData(Coverage_data))
head(evaluation.table)
evaluation.table$region<-paste(evaluation.table$chr,evaluation.table$start,evaluation.table$end,sep = "_")
evaluation_data<-evaluation.table[,c("region",colnames(evaluation.table)[grepl("coverage",colnames(evaluation.table))])]
rownames(evaluation_data)<-evaluation_data$region;evaluation_data<-evaluation_data[,-c(1)]
colnames(evaluation_data)<-Coverage_data@sample.ids
tail(evaluation_data)
bin1bp_coverage_number<-colMeans(evaluation_data,na.rm = T)
bin1bp_detected_number<-data.frame(bin1bp_coverage_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)
#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
bin1bp_detected_number2$bin1bp_number<-round(bin1bp_detected_number2$bin1bp_number,2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=2,position=position_stack(1.02),angle=0)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_average_coverage_depth_All_60_samples_meth_level_1bp_CpG_3X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_3X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_3X)<-c("sample","coverage_3X")

### for 5X ###
depth<-"5X"
Coverage_data<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_5X_coverage_bismark_1sample_meth2.rds")
#Coverage_data2<-getData(Coverage_data)
#filtering and arrange data
evaluation.table<-as.data.frame(getData(Coverage_data))
head(evaluation.table)
evaluation.table$region<-paste(evaluation.table$chr,evaluation.table$start,evaluation.table$end,sep = "_")
evaluation_data<-evaluation.table[,c("region",colnames(evaluation.table)[grepl("coverage",colnames(evaluation.table))])]
rownames(evaluation_data)<-evaluation_data$region;evaluation_data<-evaluation_data[,-c(1)]
colnames(evaluation_data)<-Coverage_data@sample.ids
tail(evaluation_data)
bin1bp_coverage_number<-colMeans(evaluation_data,na.rm = T)
bin1bp_detected_number<-data.frame(bin1bp_coverage_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)
#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
bin1bp_detected_number2$bin1bp_number<-round(bin1bp_detected_number2$bin1bp_number,2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=2,position=position_stack(1.02),angle=0)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_average_coverage_depth_All_60_samples_meth_level_1bp_CpG_5X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_5X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_5X)<-c("sample","coverage_5X")


### for 6X ###
depth<-"6X"
Coverage_data<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_6X_coverage_bismark_1sample_meth2.rds")
#Coverage_data2<-getData(Coverage_data)
#filtering and arrange data
evaluation.table<-as.data.frame(getData(Coverage_data))
head(evaluation.table)
evaluation.table$region<-paste(evaluation.table$chr,evaluation.table$start,evaluation.table$end,sep = "_")
evaluation_data<-evaluation.table[,c("region",colnames(evaluation.table)[grepl("coverage",colnames(evaluation.table))])]
rownames(evaluation_data)<-evaluation_data$region;evaluation_data<-evaluation_data[,-c(1)]
colnames(evaluation_data)<-Coverage_data@sample.ids
tail(evaluation_data)
bin1bp_coverage_number<-colMeans(evaluation_data,na.rm = T)
bin1bp_detected_number<-data.frame(bin1bp_coverage_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)
#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
bin1bp_detected_number2$bin1bp_number<-round(bin1bp_detected_number2$bin1bp_number,2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=2,position=position_stack(1.02),angle=0)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_average_coverage_depth_All_60_samples_meth_level_1bp_CpG_6X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_6X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_6X)<-c("sample","coverage_6X")

### for 10X ###
depth<-"10X"
Coverage_data<-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_1bp_CpG_10X_coverage_bismark_1sample_meth2.rds")
#Coverage_data2<-getData(Coverage_data)
#filtering and arrange data
evaluation.table<-as.data.frame(getData(Coverage_data))
head(evaluation.table)
evaluation.table$region<-paste(evaluation.table$chr,evaluation.table$start,evaluation.table$end,sep = "_")
evaluation_data<-evaluation.table[,c("region",colnames(evaluation.table)[grepl("coverage",colnames(evaluation.table))])]
rownames(evaluation_data)<-evaluation_data$region;evaluation_data<-evaluation_data[,-c(1)]
colnames(evaluation_data)<-Coverage_data@sample.ids
tail(evaluation_data)

bin1bp_coverage_number<-colMeans(evaluation_data,na.rm = T)
bin1bp_detected_number<-data.frame(bin1bp_coverage_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin1bp_detected_number) <- "bin1bp_number"
rownames(bin1bp_detected_number);rownames(huahua_meta)
bin1bp_detected_number2<-merge(bin1bp_detected_number,huahua_meta,by=0)
head(bin1bp_detected_number2);dim(bin1bp_detected_number2)
bin1bp_detected_number2$group2<-paste(bin1bp_detected_number2$Age_group,bin1bp_detected_number2$Sample.types,sep="_")

bin1bp_detected_number2$Age_group <-factor(bin1bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin1bp_detected_number2$Sample.types <-factor(bin1bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin1bp_detected_number2$group2 <-factor(bin1bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin1bp_detected_number2$Gender_RRBS<-factor(bin1bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-as.character(bin1bp_detected_number2$analysis_name)
#bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$group2,bin1bp_detected_number2$bin1bp_number,decreasing = T),]
bin1bp_detected_number2<-bin1bp_detected_number2[order(bin1bp_detected_number2$analysis_name,decreasing = F),]
head(bin1bp_detected_number2)
bin1bp_detected_number2$analysis_name<-factor(bin1bp_detected_number2$analysis_name,levels=as.character(bin1bp_detected_number2$analysis_name))
head(bin1bp_detected_number2)
bin1bp_detected_number2$bin1bp_number<-round(bin1bp_detected_number2$bin1bp_number,2)
#plot
plot_bin1bp_detected_number0<-ggplot(data=bin1bp_detected_number2, mapping=aes(x=analysis_name,y=bin1bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin1bp_detected_number1 <-plot_bin1bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin1bp_number), color="black", size=2,position=position_stack(1.02),angle=0)+
  theme_classic()+labs(x="",y="Number of CpGs",title=paste0("Number of CpGs detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin1bp_detected_number1
ggsave(plot_bin1bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_average_coverage_depth_All_60_samples_meth_level_1bp_CpG_10X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin1bp_detected_number2)
CpG_number_10X<-bin1bp_detected_number2[,c("analysis_name","bin1bp_number")]
colnames(CpG_number_10X)<-c("sample","coverage_10X")
head(bin1bp_detected_number2)

### formational analysis ###
depth<-"200bp_6X_3L_8s"
Coverage_data <-readRDS(file = "/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth2.rds")
#Coverage_data2<-getData(Coverage_data)
#filtering and arrange data
evaluation.table<-as.data.frame(getData(Coverage_data))
head(evaluation.table)
evaluation.table$region<-paste(evaluation.table$chr,evaluation.table$start,evaluation.table$end,sep = "_")
evaluation_data<-evaluation.table[,c("region",colnames(evaluation.table)[grepl("coverage",colnames(evaluation.table))])]
rownames(evaluation_data)<-evaluation_data$region;evaluation_data<-evaluation_data[,-c(1)]
colnames(evaluation_data)<-Coverage_data@sample.ids
tail(evaluation_data)
colMaxs(as.matrix(evaluation_data),na.rm = T)
colMins(as.matrix(evaluation_data),na.rm = T)

bin200bp_coverage_number<-colMeans(evaluation_data,na.rm = T)
bin200bp_detected_number<-data.frame(bin200bp_coverage_number)

head(huahua_meta);dim(huahua_meta)
colnames(bin200bp_detected_number) <- "bin200bp_number"
rownames(bin200bp_detected_number);rownames(huahua_meta)
bin200bp_detected_number2<-merge(bin200bp_detected_number,huahua_meta,by=0)
head(bin200bp_detected_number2);dim(bin200bp_detected_number2)
bin200bp_detected_number2$group2<-paste(bin200bp_detected_number2$Age_group,bin200bp_detected_number2$Sample.types,sep="_")

bin200bp_detected_number2$Age_group <-factor(bin200bp_detected_number2$Age_group,level=c("YOUNG","AMA"))
bin200bp_detected_number2$Sample.types <-factor(bin200bp_detected_number2$Sample.types,level=c("kids","mother","father"))
bin200bp_detected_number2$group2 <-factor(bin200bp_detected_number2$group2,level=c("YOUNG_kids","AMA_kids","YOUNG_mother","AMA_mother","YOUNG_father","AMA_father"))
#bin200bp_detected_number2$Gender_RRBS<-factor(bin200bp_detected_number2$Gender_RRBS,levels = c("Female","Male"))
head(bin200bp_detected_number2)
bin200bp_detected_number2$analysis_name<-as.character(bin200bp_detected_number2$analysis_name)
#bin200bp_detected_number2<-bin200bp_detected_number2[order(bin200bp_detected_number2$group2,bin200bp_detected_number2$bin200bp_number,decreasing = T),]
bin200bp_detected_number2<-bin200bp_detected_number2[order(bin200bp_detected_number2$analysis_name,decreasing = F),]
head(bin200bp_detected_number2)
bin200bp_detected_number2$analysis_name<-factor(bin200bp_detected_number2$analysis_name,levels=as.character(bin200bp_detected_number2$analysis_name))
head(bin200bp_detected_number2)
bin200bp_detected_number2$bin200bp_number<-round(bin200bp_detected_number2$bin200bp_number,2)
#plot
plot_bin200bp_detected_number0<-ggplot(data=bin200bp_detected_number2, mapping=aes(x=analysis_name,y=bin200bp_number,fill=group2))+geom_bar(stat="identity",width=0.8)
plot_bin200bp_detected_number1 <-plot_bin200bp_detected_number0+scale_fill_manual(values=ppCor[c(7:8,2:1,3,5)])+
  geom_text(stat="identity",aes(label=bin200bp_number), color="black", size=2,position=position_stack(1.02),angle=0)+
  theme_classic()+labs(x="",y="Coverage depth",title=paste0("Coverage depth detected_in_covarage_",depth))+guides(fill=guide_legend(ncol=1)) +
  theme(legend.text = element_text(size = 10, colour = "black"),axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),axis.text.y  = element_text(size = 10,colour = 'black'),
        axis.text.x = element_text(size = 9,colour = 'black',vjust=0.5,hjust=1,angle = 90),legend.title = element_text(size = 9))
plot_bin200bp_detected_number1
ggsave(plot_bin200bp_detected_number1,file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_average_coverage_depth_All_60_samples_meth_level_200bp_CpG_10X_coverage_bismark_1sample.pdf",width = 16, height = 6)

head(bin200bp_detected_number2)
CpG_number_200bp<-bin200bp_detected_number2[,c("analysis_name","bin200bp_number")]
colnames(CpG_number_200bp)<-c("sample","200bp_6X_3L_8s")

CpG_number_data<-merge(CpG_number_1X,CpG_number_3X,by="sample")
CpG_number_data<-merge(CpG_number_data,CpG_number_5X,by="sample")
CpG_number_data<-merge(CpG_number_data,CpG_number_6X,by="sample")
CpG_number_data<-merge(CpG_number_data,CpG_number_10X,by="sample")
CpG_number_data<-merge(CpG_number_data,CpG_number_200bp,by="sample")

head(CpG_number_data)
write.table(as.data.frame(CpG_number_data), file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/response_CpGs_average_coverage_depth_All_60_samples_number_of_1bp_CpG_different_coverage_bismark_1sample.txt",quote=F, row.names=F, col.names=T) 

########################
#mean DNA methylation level after filters
colnames(data_raw_bin)
data_raw_bin$AMA_Kids_mean<-apply(data_raw_bin[,1:10],1,function(x) mean(x,na.rm=T))
data_raw_bin$YOUNG_Kids_mean<-apply(data_raw_bin[,11:20],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_Mother_mean<-apply(data_raw_bin[,21:30],1,function(x) mean(x,na.rm=T))
data_raw_bin$YOUNG_Mother_mean<-apply(data_raw_bin[,31:40],1,function(x) mean(x,na.rm=T))
data_raw_bin$AMA_father_mean<-apply(data_raw_bin[,41:50],1,function(x) mean(x,na.rm=T))
data_raw_bin$YOUNG_father_mean<-apply(data_raw_bin[,51:60],1,function(x) mean(x,na.rm=T))
dim(data_raw_bin)
write.table(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.txt",quote=F,row.names=T)
write.csv(data_raw_bin, file="/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_meth.csv",quote=F, row.names=T)
dim(data_raw_bin[,61:66])#701937      6
head(data_raw_bin[,61:66])#701937      6

pdf("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/All_60_samples_meth_level_CpG_6X_meth_bismark_200bp_3L_8sample_similarity_logisty.pdf")
p2<-ggplot(data = na.omit(data_raw_bin[,61:62]), mapping = aes(x = AMA_Kids_mean, y = YOUNG_Kids_mean)) + 
  stat_bin2d(bins = 300) + scale_fill_gradient(low = 'steelblue', high = 'darkred', limits = c(0,100), breaks = c(0,25,50,100))
p_plot1 <- p2+ stat_cor(method = "pearson", label.x = 10, label.y = 90)+xlim(0,100)+ xlab("AMA_Kids_mean") + 
  theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5,  hjust = 0.5))+ylab("YOUNG_Kids_mean") + 
  theme(axis.title.y = element_text(size = 16,face = "bold", vjust = 0.5, hjust = 0.5))+theme_bw()
print(p_plot1)

p4<-ggplot(data = na.omit(data_raw_bin[,63:64]), mapping = aes(x = AMA_Mother_mean, y = YOUNG_Mother_mean)) + 
  stat_bin2d(bins = 300) + scale_fill_gradient(low = 'steelblue', high = 'darkred', limits = c(0,100), breaks = c(0,25,50,100))
p_plot3 <- p4+ stat_cor(method = "pearson", label.x = 10, label.y = 90)+xlim(0,100)+ xlab("AMA_Mother_mean") + 
  theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5,  hjust = 0.5))+ylab("YOUNG_Mother_mean") + 
  theme(axis.title.y = element_text(size = 16,face = "bold", vjust = 0.5, hjust = 0.5))+theme_bw()
print(p_plot3)

p3<-ggplot(data = na.omit(data_raw_bin[,65:66]), mapping = aes(x = AMA_father_mean, y = YOUNG_father_mean)) + 
  stat_bin2d(bins = 300) + scale_fill_gradient(low = 'steelblue', high = 'darkred', limits = c(0,100), breaks = c(0,25,50,100))
p_plot4 <- p3+ stat_cor(method = "pearson", label.x = 10, label.y = 90)+xlim(0,100)+ xlab("AMA_father_mean") + 
  theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5,  hjust = 0.5))+ylab("YOUNG_father_mean") + 
  theme(axis.title.y = element_text(size = 16,face = "bold", vjust = 0.5, hjust = 0.5))+theme_bw()
print(p_plot4)

dev.off()