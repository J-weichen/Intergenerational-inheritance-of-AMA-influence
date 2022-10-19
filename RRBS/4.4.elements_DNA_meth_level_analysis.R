rm(list = ls())
.libPaths("/mnt/data/chenwei/software/R_lib") 
dyn.load('/mnt/data/chenwei/software/glpk/libglpk.so.40')
library(ggplot2)
library(ggpubr) 
library(methylKit)
library(data.table)
library(scales)
library(ggsci)
library(stringr)
#调颜色
pal <- pal_npg("nrc", alpha=1)(9)#nrc是Palette Types，alpha用于调节透明度
nejm<-pal_nejm("default",alpha = 1)(8)##(9表示呈现多少个颜色)nejm
ppCor <-c(nejm,pal[c(2,5,9,1,4,6,7,8)])
show_col(ppCor)

#step2: read and prepare region files
files_region_hg38<- list.files("/mnt/data/chenwei/qinmen_BR/00.ref_data/1-30regions/chenwei_file")
data_list_region<-list();cnane_region<-c()
#data_list_bedGR<-list()

for ( f_name in files_region_hg38[1:length(files_region_hg38)]){
  #f_name<-"ICRs_380-23chr.bed" ##test line 
  FF<-read.table(paste("/mnt/data/chenwei/qinmen_BR/00.ref_data/1-30regions/chenwei_file/",f_name,sep=""), header=F)
  FF_bed<-FF[,c(1:3)]
  head(FF_bed);dim(FF_bed)
  FF_bed<-unique(FF_bed);dim(FF_bed)
  colnames(FF_bed)<-c("chr","start","end")
  FF_bed<-FF_bed[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",FF_bed$chr,invert=TRUE),]
  FF_bed <- data.frame(FF_bed,V4=paste0("tile",c(1:nrow(FF_bed))))
  head(FF_bed);dim(FF_bed)
  FF_bed$chr<-factor(FF_bed$chr,levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                           "chr20","chr21","chr22"))
  data_list_region<-c(data_list_region,list(as.data.frame(FF_bed)))
  print(as.character(f_name))
  print(dim(as.data.frame(FF_bed)))
  f_name1<- unlist(lapply(strsplit(f_name,"[.]"), function(x) x[1]))
  f_name2<- unlist(lapply(strsplit(f_name1,"-"), function(x) x[1]))
  cnane_region<-c(cnane_region,f_name2)
  #cnane_single_father<-c(cnane_region,unlist(lapply(strsplit(f_name2,"-"), function(x) paste("single",x[2],x[3],x[4],sep="-"))))
}
length(data_list_region);length(cnane_region)
names(data_list_region)<-cnane_region
head(data_list_region$LTR)
head(data_list_region[["LTR"]])

##step3: calculate DNA methylation level for selected regions 
target_region_list<-c("Low_complexity","LTR","LINE","SINE","Retroposon_SVA","Transposon",
                      "Satellite","Micro_satellites","CpGislands","ICRs_67","MEs",
                      "HCP","ICP","LCP","UTR5","UTR3","Exons","Introns","All_Genes","Intergenic_region_KSM")

#step3 reading cov files

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

sample_list<-unlist(c(AMA_Kids_files2,YOUNG_Kids_files2,AMA_Mother_files2,YOUNG_Mother_files2,AMA_father_files2,YOUNG_father_files2))

#step 4 calculated methylation level in regions
data_list_region_meth_2<-c()
data_list_region_meth_n2<-c()
data_list_region_meth_common2<-c()
cnane_region_meth<-c()

for (region_name in target_region_list){
  # region_name<-"MEs"   ##test line 
  print(as.character(region_name))
  FF<-data_list_region[[which(names(data_list_region) == region_name)]]
  data_list_region_meth_1<-c()
  region_list<-c()
  interSE<-list()
  
  for(sample in sample_list){
    # sample <- "/mnt/data/chenwei/huahua/3.mapped_data/bismark_cov_file/YOUNG_file/Father/E16F_1_CpG.bismark.cov"
    Sample_Name<-as.character(unlist(lapply(strsplit(sample,"/"), function(x) x[10])))
    print(Sample_Name)
    Sample_Name2<-as.character(unlist(lapply(strsplit(Sample_Name,"_"), function(x) x[1])))
   # data_meth1<-read.table(sample, header=F, stringsAsFactors=F)
    data_meth1 <- fread(sample,header = F,fill=TRUE, na.strings = "",stringsAsFactors=F)
    head(data_meth1);head(FF)
    #for CpGs with no less
    than six corvarage
    data_meth2<-as.data.table(data_meth1[which((data_meth1$V5 + data_meth1$V6)>=6),c(1:4)])
    #str(data_meth2)
    FF2<-as.data.table(FF[,c(1:3)])
    colnames(data_meth2)<-c("chr","start","end","meth");colnames(FF2)<-c("chr","start","end")
    
    # head(data_meth2)
    setkey(FF2, chr, start, end)
    data_meth3<-as.data.frame(foverlaps(data_meth2,FF2,  nomatch=0))
    #dim(data_meth3);head(data_meth3)
    data_meth3<-data_meth3[grep("chrUn_*|*_alt|*random|chrM|chrY|chrX",data_meth3$chr,invert=TRUE),]
    #dim(data_meth3)#3742    7
    data_meth3$ID <- paste(data_meth3$chr,data_meth3$start,data_meth3$end,sep = "_")
    length(unique(data_meth3$ID))# 363
    ##select  region with more than 3 CpGs
    CpG_num<-table(data_meth3$ID)
    data_meth4<-data_meth3[which(data_meth3$ID %in% names(CpG_num[CpG_num>=3])),]
    head(data_meth4);dim(data_meth4)#3696    7
    data_meth5 <-aggregate(data_meth4[,6],list(data_meth4[,7]),FUN=mean, na.rm=TRUE, na.action=NULL)
    colnames(data_meth5)<-c("region",Sample_Name2)
    interSE<-c(interSE,list(as.data.frame(data_meth5)))
    region_num<-nrow(data_meth5)
    data_meth_mean<-mean(data_meth5[,2],na.rm=T)
    # print(data_meth_mean)
    names(data_meth_mean)<-paste0(region_name,"-",Sample_Name2)
    names(region_num)<-paste0(region_name,"-",Sample_Name2)
    data_list_region_meth_1<-c(data_list_region_meth_1,data_meth_mean)
    region_list<-c(region_list,region_num)
  
    } 

saveRDS(interSE,paste0("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/",region_name,"_all_sample_meth_data.rds"))
  
#for common region
commonlist<-Reduce(intersect,list(interSE[[1]]$region,interSE[[2]]$region,interSE[[3]]$region,interSE[[4]]$region,interSE[[5]]$region,interSE[[6]]$region,interSE[[7]]$region,
                                  interSE[[8]]$region,interSE[[9]]$region,interSE[[10]]$region,interSE[[11]]$region,interSE[[12]]$region,interSE[[13]]$region,interSE[[14]]$region,
                                  interSE[[15]]$region,interSE[[16]]$region,interSE[[17]]$region,interSE[[18]]$region,interSE[[19]]$region,interSE[[20]]$region,
                                  interSE[[21]]$region,interSE[[22]]$region,interSE[[23]]$region,interSE[[24]]$region,interSE[[25]]$region,interSE[[26]]$region,interSE[[27]]$region,
                                  interSE[[28]]$region,interSE[[29]]$region,interSE[[30]]$region,interSE[[31]]$region,interSE[[32]]$region,interSE[[33]]$region,interSE[[34]]$region,
                                  interSE[[35]]$region,interSE[[36]]$region,interSE[[37]]$region,interSE[[38]]$region,interSE[[39]]$region,interSE[[40]]$region,interSE[[41]]$region,
                                  interSE[[42]]$region,interSE[[43]]$region,interSE[[44]]$region,interSE[[45]]$region,interSE[[46]]$region,interSE[[47]]$region,
                                  interSE[[48]]$region,interSE[[49]]$region,interSE[[50]]$region,interSE[[51]]$region,interSE[[52]]$region,interSE[[53]]$region,interSE[[54]]$region,
                                  interSE[[55]]$region,interSE[[56]]$region,interSE[[57]]$region,interSE[[58]]$region,interSE[[59]]$region,interSE[[60]]$region
))

all_list<-unique(c(interSE[[1]]$region,interSE[[2]]$region,interSE[[3]]$region,interSE[[4]]$region,interSE[[5]]$region,interSE[[6]]$region,interSE[[7]]$region,
                   interSE[[8]]$region,interSE[[9]]$region,interSE[[10]]$region,interSE[[11]]$region,interSE[[12]]$region,interSE[[13]]$region,interSE[[14]]$region,
                   interSE[[15]]$region,interSE[[16]]$region,interSE[[17]]$region,interSE[[18]]$region,interSE[[19]]$region,interSE[[20]]$region,
                   interSE[[21]]$region,interSE[[22]]$region,interSE[[23]]$region,interSE[[24]]$region,interSE[[25]]$region,interSE[[26]]$region,interSE[[27]]$region,
                   interSE[[28]]$region,interSE[[29]]$region,interSE[[30]]$region,interSE[[31]]$region,interSE[[32]]$region,interSE[[33]]$region,interSE[[34]]$region,
                   interSE[[35]]$region,interSE[[36]]$region,interSE[[37]]$region,interSE[[38]]$region,interSE[[39]]$region,interSE[[40]]$region,interSE[[41]]$region,
                   interSE[[42]]$region,interSE[[43]]$region,interSE[[44]]$region,interSE[[45]]$region,interSE[[46]]$region,interSE[[47]]$region,
                   interSE[[48]]$region,interSE[[49]]$region,interSE[[50]]$region,interSE[[51]]$region,interSE[[52]]$region,interSE[[53]]$region,interSE[[54]]$region,
                   interSE[[55]]$region,interSE[[56]]$region,interSE[[57]]$region,interSE[[58]]$region,interSE[[59]]$region,interSE[[60]]$region
))

common_data_mean<-c(mean(interSE[[1]][which(interSE[[1]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[2]][which(interSE[[2]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[3]][which(interSE[[3]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[4]][which(interSE[[4]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[5]][which(interSE[[5]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[6]][which(interSE[[6]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[7]][which(interSE[[7]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[8]][which(interSE[[8]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[9]][which(interSE[[9]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[10]][which(interSE[[10]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[11]][which(interSE[[11]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[12]][which(interSE[[12]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[13]][which(interSE[[13]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[14]][which(interSE[[14]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[15]][which(interSE[[15]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[16]][which(interSE[[16]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[17]][which(interSE[[17]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[18]][which(interSE[[18]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[19]][which(interSE[[19]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[20]][which(interSE[[20]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[21]][which(interSE[[21]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[22]][which(interSE[[22]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[23]][which(interSE[[23]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[24]][which(interSE[[24]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[25]][which(interSE[[25]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[26]][which(interSE[[26]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[27]][which(interSE[[27]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[28]][which(interSE[[28]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[29]][which(interSE[[29]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[30]][which(interSE[[30]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[31]][which(interSE[[31]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[32]][which(interSE[[32]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[33]][which(interSE[[33]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[34]][which(interSE[[34]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[35]][which(interSE[[35]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[36]][which(interSE[[36]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[37]][which(interSE[[37]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[38]][which(interSE[[38]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[39]][which(interSE[[39]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[40]][which(interSE[[40]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[41]][which(interSE[[41]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[42]][which(interSE[[42]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[43]][which(interSE[[43]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[44]][which(interSE[[44]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[45]][which(interSE[[45]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[46]][which(interSE[[46]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[47]][which(interSE[[47]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[48]][which(interSE[[48]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[49]][which(interSE[[49]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[50]][which(interSE[[50]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[51]][which(interSE[[51]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[52]][which(interSE[[52]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[53]][which(interSE[[53]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[54]][which(interSE[[54]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[55]][which(interSE[[55]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[56]][which(interSE[[56]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[57]][which(interSE[[57]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[58]][which(interSE[[58]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[59]][which(interSE[[59]]$region %in% commonlist ),2],na.rm=T),
                    mean(interSE[[60]][which(interSE[[60]]$region %in% commonlist ),2],na.rm=T))
  # names(common_data_mean)<-names(region_list)
  names(common_data_mean)<- names(region_list)
  
  print(data_list_region_meth_1)
  print(c(region_list,length(commonlist),length(all_list)))
  print(common_data_mean)
  
  data_list_region_meth_2<-c(data_list_region_meth_2,data_list_region_meth_1)
  data_list_region_meth_n2<-c(data_list_region_meth_n2,region_list)
  data_list_region_meth_common2<-c(data_list_region_meth_common2,common_data_mean)
}

#for all regions respective for individual
mean_data<-as.data.frame(data_list_region_meth_2)
colnames(mean_data)<-"meth_mean"
mean_data_common<-as.data.frame(data_list_region_meth_common2)
colnames(mean_data_common)<-"meth_mean"
head(mean_data);head(mean_data_common)
rownames(mean_data)
saveRDS(mean_data,"/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/all_region_mean_data.rds")
saveRDS(mean_data_common,"/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/all_region_mean_data_common.rds")
saveRDS(as.data.frame(data_list_region_meth_n2),"/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/number_of_all_region_mean_data.rds")

########re-orgenization
mean_data$sample_region <- rownames(mean_data)
#mean_data$sample <-unlist(lapply(strsplit(mean_data$sample_region,"-"), function(x) x[2]))
tail(mean_data$sample_region)
mean_data$sample<-gsub("_1_CpG.bismark.cov", "",mean_data$sample_region)
mean_data$sample <-unlist(lapply(strsplit(mean_data$sample,"-"), function(x) paste(x[2],x[3],x[4],x[5],sep = "-")))
mean_data$sample <- gsub("-NA","",mean_data$sample)
mean_data$region <-unlist(lapply(strsplit(mean_data$sample_region,"-"), function(x) x[1]))

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
huahua_meta$sample<-huahua_meta$library_code
head(huahua_meta)
mean_data2<-merge(mean_data,huahua_meta,by="sample")
head(mean_data2)
mean_data2$group_region<-paste(mean_data2$region,mean_data2$Age_group,sep="-")
head(mean_data2,n=7)
dim(mean_data2)
table(mean_data2$region)
write.csv(mean_data2,"/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/AMA_all_60sample_6X_CpG_DNA_methylation_in_target_region_all_bin.csv")

mean_data_common$sample_region <- rownames(mean_data_common)
head(mean_data_common)
#mean_data$sample <-unlist(lapply(strsplit(mean_data$sample_region,"-"), function(x) x[2]))
tail(mean_data_common$sample_region)
mean_data_common$sample<-gsub("_1_CpG.bismark.cov", "",mean_data_common$sample_region)
mean_data_common$sample <-unlist(lapply(strsplit(mean_data_common$sample,"-"), function(x) paste(x[2],x[3],x[4],x[5],sep = "-")))
mean_data_common$sample <- gsub("-NA","",mean_data_common$sample)
mean_data_common$region <-unlist(lapply(strsplit(mean_data_common$sample_region,"-"), function(x) x[1]))

huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
huahua_meta$sample<-huahua_meta$library_code
head(huahua_meta)
mean_data_common2<-merge(mean_data_common,huahua_meta,by="sample")
head(mean_data_common2)
mean_data_common2$group_region<-paste(mean_data_common2$region,mean_data_common2$Age_group,sep="-")
head(mean_data_common2,n=7)
dim(mean_data_common2)
table(mean_data_common2$region)
write.csv(mean_data_common2,"/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/AMA_all_60sample_6X_CpG_DNA_methylation_in_target_region_common_bin.csv")

##save the number of targeted regions covered with more than 3 CpGs in  with no less than six corvarage
data_list_region_meth_n2<-readRDS("/mnt/data/chenwei/huahua/4.methy_result/0.global_compasion/number_of_all_region_mean_data.rds")
str(data_list_region_meth_n2);head(data_list_region_meth_n2)
data_list_region_meth_n2$region<- as.character(unlist(lapply(strsplit(rownames(data_list_region_meth_n2),"-"), function(x) x[1])))
data_list_region_meth_n2$sample<- as.character(unlist(lapply(strsplit(rownames(data_list_region_meth_n2),"-"), function(x) paste(x[2],x[3],x[4],sep="-"))))
data_list_region_meth_n2$sample<- gsub("-NA","",data_list_region_meth_n2$sample) 
unique(data_list_region_meth_n2$sample)
head(data_list_region_meth_n2)
huahua_meta <-read.csv(file="/mnt/data/chenwei/huahua/0.hua_script/AMA_analysis_metadata.csv", header = T,row.names= 1)
huahua_meta$sample<-huahua_meta$library_code
head(huahua_meta)
mean_data<-merge(data_list_region_meth_n2,huahua_meta,by="sample")
head(mean_data)
mean_data2<-mean_data[,c("analysis_name","region","data_list_region_meth_n2" )]
mean_data3<-dcast(mean_data2,analysis_name~region)
head(mean_data3)
write.csv(mean_data3,"/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/The_number_of_element_covered_3more_6X_CpGs_for_all_60sample.csv")

## plot target regions DNA methylation
##for all region
mean_data2<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/AMA_all_60sample_6X_CpG_DNA_methylation_in_target_region_all_bin.csv",row.names = 1)
head(mean_data2,n=7)
unique(as.character(mean_data2$group_region))
mean_data2$group_age<-paste(mean_data2$Sample.types,mean_data2$Age_group,sep="-")
#target_region_list<-c("Low_complexity","LTR","LINE","SINE","Retroposon_SVA","Transposon",
#                      "Satellite","Micro_satellites","CpGislands","ICRs_380","ICRs_67","MEs",
#                      "Intragenic_region_KSM","Intergenic_region_KSM",
#                      "HCP","ICP","LCP","UTR5","UTR3","Exons","Introns","All_Genes","Intergenic_1kb")
target_region_list2<-c("Low_complexity","CpGislands","HCP","ICP","LCP","UTR5","UTR3","All_Genes",
                       "LTR", "LINE","SINE","Retroposon_SVA","Transposon","Satellite","Micro_satellites","MEs","ICRs_67")
plot_mean_data<-mean_data2[which(mean_data2$region %in% target_region_list2),]
dim(plot_mean_data)
unique(as.character(plot_mean_data$group_region))

#for kids
kids_list<-c("10-3-1E","16-3-4D","17-3-2F","18-3","19-3-4E","20Q-4H","6-3","7-3","8-3","E16C",
               "11-3-3F","12-3-4A","13-3-1G","14-3","15-3-2C","1Q","2-3-3C","3-3","4-3","5-3")
#for mother
mother_list<-c("10-1-in-9","16-1","17-1-2D","18-1","19-1-2G","20M-4F","6-1","7-1","8-1","E16M",
               "11-1","12-1-3G","13-1-1F","14-1","15-1","1M","2-1","3M","4-1","5-1")
#for father
father_list<-c("10-2-1D","16-2-4C","17-2-2E","18-2","19-2-2H","20F-4G","6-2","7F","8-2","E16F",
               "11-2","12-2-3H","13-2-1H","14-2","15-2","1F","2-2","3F","4-2","5-2")

#YOUNG_sample<-c("P11","P12","P13","P14","P15","P31","P32","P33","P34","P35","P36","P41",  
#                   "P16","P17","P18","P19","P20","P21","P22","P24","P25","P38")
#AMA_sample<-c("A10","A11","A36","A6","A8","A9","A1","A22","A25","A3","A31","A32","A33","A4")
plot_mean_data$sample<-factor(plot_mean_data$sample,levels=c(kids_list,mother_list,father_list))
plot_mean_data$analysis_name<-factor(plot_mean_data$analysis_name,levels=c(paste0("YOUNG_K_",1:10),paste0("AMA_K_",1:10),paste0("YOUNG_M_",1:10),paste0("AMA_M_",1:10),paste0("YOUNG_F_",1:10),paste0("AMA_F_",1:10)))
plot_mean_data$group_age<-factor(plot_mean_data$group_age,levels=c("kids-YOUNG","kids-AMA","mother-YOUNG","mother-AMA","father-YOUNG","father-AMA"))
#plot_mean_data$group_region<-factor(plot_mean_data$group_region,levels=c("Low_complexity-YOUNG","Low_complexity-AMA","CpGislands-YOUNG","CpGislands-AMA",
#                                                                         "UTR5-YOUNG","UTR5-AMA","HCP-YOUNG","HCP-AMA","ICP-YOUNG","ICP-AMA",
#                                                                         "LCP-YOUNG","LCP-AMA","Intragenic_region_KSM-YOUNG","Intragenic_region_KSM-AMA",
#                                                                         "Intergenic_region_KSM-YOUNG","Intergenic_region_KSM-AMA","Intergenic_1kb-YOUNG","Intergenic_1kb-AMA",
#                                                                         "Exons-YOUNG","Exons-AMA","Introns-YOUNG","Introns-AMA",
#                                                                         "UTR3-YOUNG","UTR3-AMA", "All_Genes-YOUNG","All_Genes-AMA",
#                                                                         "LTR-YOUNG","LTR-AMA", "LINE-YOUNG","LINE-AMA", "SINE-YOUNG","SINE-AMA",
#                                                                         "Retroposon_SVA-YOUNG","Retroposon_SVA-AMA","Transposon-YOUNG","Transposon-AMA",
#                                                                         "Satellite-YOUNG","Satellite-AMA","Micro_satellites-YOUNG","Micro_satellites-AMA",
#                                                                         "MEs-YOUNG","MEs-AMA","ICRs_380-YOUNG","ICRs_380-AMA","ICRs_67-YOUNG","ICRs_67-AMA"))

plot_mean_data$group_region<-factor(plot_mean_data$group_region,levels=c("Low_complexity-YOUNG","Low_complexity-AMA","CpGislands-YOUNG","CpGislands-AMA",
                                                                         "HCP-YOUNG","HCP-AMA","ICP-YOUNG","ICP-AMA","LCP-YOUNG","LCP-AMA", "UTR5-YOUNG","UTR5-AMA",
                                                                         "UTR3-YOUNG","UTR3-AMA", "All_Genes-YOUNG","All_Genes-AMA",
                                                                         "LTR-YOUNG","LTR-AMA", "LINE-YOUNG","LINE-AMA", "SINE-YOUNG","SINE-AMA",
                                                                         "Retroposon_SVA-YOUNG","Retroposon_SVA-AMA","Transposon-YOUNG","Transposon-AMA",
                                                                         "Satellite-YOUNG","Satellite-AMA","Micro_satellites-YOUNG","Micro_satellites-AMA",
                                                                         "MEs-YOUNG","MEs-AMA","ICRs_67-YOUNG","ICRs_67-AMA"))
plot_mean_data$group_region
head(plot_mean_data)    
#respond to reviewer
average_meth<-aggregate(plot_mean_data[,c("Age_group","region","meth_mean")],by=list(plot_mean_data$Age_group,plot_mean_data$region),FUN=mean)
average_meth[which(average_meth$Group.2 == "CpGislands"),]
##   Age_group region meth_mean
#     AMA CpGislands    26.4950
#   YOUNG CpGislands    26.5202

P_mean_bar_1<-ggplot(data=plot_mean_data,aes(x=group_region,y= meth_mean))+
  geom_bar(aes(fill=sample),stat = "identity", position ="dodge", width=0.6)+# position_dodge2(padding = 0.1,preserve = "single")
  scale_fill_manual(values=c(rep(ppCor[7],10),rep(ppCor[8],10),rep(ppCor[2],10),rep(ppCor[1],10),rep(ppCor[3],10),rep(ppCor[5],10)))+
  labs(title = "all region and CpGs served")+xlab("") +ylab("DNA methylation level(mean)")+#coord_flip()+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))
P_mean_bar_2<-ggplot(data=plot_mean_data,aes(x=group_region,y= meth_mean))+
  geom_bar(aes(fill=analysis_name),stat = "identity", position ="dodge", width=0.6)+# position_dodge2(padding = 0.1,preserve = "single")
  scale_fill_manual(values=c(rep(ppCor[7],10),rep(ppCor[8],10),rep(ppCor[2],10),rep(ppCor[1],10),rep(ppCor[3],10),rep(ppCor[5],10)))+
  labs(title = "all region and CpGs served")+xlab("") +ylab("DNA methylation level(mean)")+#coord_flip()+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))
P_mean_bar_1_1<-P_mean_bar_1+ facet_grid(.~Sample.types)
P_mean_bar_2_1<-P_mean_bar_2+ facet_grid(.~Sample.types)
P_mean_bar_1_2<-P_mean_bar_1+ facet_wrap(~Sample.types,scales="free",nrow=3)
P_mean_bar_2_2<-P_mean_bar_2+ facet_wrap(~Sample.types,scales="free",nrow=3)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/bar_meth_means_in_each_element_for_six_groups_all_bin.pdf",P_mean_bar_1_1,width=18, height=6)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/bar_meth_means_in_each_element_for_six_groups_all_bin2.pdf",P_mean_bar_2_1,width=18, height=6)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/bar_meth_means_in_each_element_for_six_groups_all_bin3.pdf",P_mean_bar_1_2,width=10, height=12)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/bar_meth_means_in_each_element_for_six_groups_all_bin4.pdf",P_mean_bar_2_2,width=10, height=12)


P_mean_site <- ggboxplot(plot_mean_data, x = "group_region", y = "meth_mean",color = "Age_group", palette = ppCor[1:2],add = "jitter",size = 0.5)+#ylim(c(0,90))+
  xlab("sample_group") +ylab("DNA_methylation")+#title("DNA methylation in final neonatal DMRs")+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(#plot.title = element_text(size=10,colour = "blue",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),
    axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
P_mean_site1<-P_mean_site+facet_wrap(~Sample.types,scales="free",nrow=3)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_all_bin.pdf",P_mean_site1,width=16, height=30)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_all_bin2.pdf",P_mean_site1,width=8, height=15)

#my_comparisons <- list(c("Low_complexity-YOUNG","Low_complexity-AMA"),c("CpGislands-YOUNG","CpGislands-AMA"),
#                       c("UTR5-YOUNG","UTR5-AMA"),c("HCP-YOUNG","HCP-AMA"),c("ICP-YOUNG","ICP-AMA"),
#                        c("LCP-YOUNG","LCP-AMA"),c("Intragenic_region_KSM-YOUNG","Intragenic_region_KSM-AMA"),
#                        c("Intergenic_region_KSM-YOUNG","Intergenic_region_KSM-AMA"),c("Intergenic_1kb-YOUNG","Intergenic_1kb-AMA"),
#                        c("Exons-YOUNG","Exons-AMA"),c("Introns-YOUNG","Introns-AMA"),
#                        c("UTR3-YOUNG","UTR3-AMA"), c("All_Genes-YOUNG","All_Genes-AMA"),
#                        c("LTR-YOUNG","LTR-AMA"),c( "LINE-YOUNG","LINE-AMA"),c("SINE-YOUNG","SINE-AMA"),
#                        c("Retroposon_SVA-YOUNG","Retroposon_SVA-AMA"),c("Transposon-YOUNG","Transposon-AMA"),
#                        c("Satellite-YOUNG","Satellite-AMA"),c("Micro_satellites-YOUNG","Micro_satellites-AMA"),
#                        c("MEs-YOUNG","MEs-AMA"),c("ICRs_380-YOUNG","ICRs_380-AMA"),c("ICRs_67-YOUNG","ICRs_67-AMA"))
my_comparisons <- list(c("Low_complexity-YOUNG","Low_complexity-AMA"),c("CpGislands-YOUNG","CpGislands-AMA"),
                       c("HCP-YOUNG","HCP-AMA"),c("ICP-YOUNG","ICP-AMA"), c("LCP-YOUNG","LCP-AMA"), 
                       c("UTR5-YOUNG","UTR5-AMA"),c("UTR3-YOUNG","UTR3-AMA"), c("All_Genes-YOUNG","All_Genes-AMA"),
                       c("LTR-YOUNG","LTR-AMA"),c( "LINE-YOUNG","LINE-AMA"),c("SINE-YOUNG","SINE-AMA"),
                       c("Retroposon_SVA-YOUNG","Retroposon_SVA-AMA"),c("Transposon-YOUNG","Transposon-AMA"),
                       c("Satellite-YOUNG","Satellite-AMA"),c("Micro_satellites-YOUNG","Micro_satellites-AMA"),
                       c("MEs-YOUNG","MEs-AMA"),c("ICRs_67-YOUNG","ICRs_67-AMA"))

#P_mean_site2 +#stat_compare_means(method = "kruskal.test",aes(label=paste0(..method..,": "," p = ",..p.format..)), label.x = 1.4, label.y =85)+
plot_region_mean1<-P_mean_site+ylim(0,110) + stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif",label.y = 100) +labs(title = "un_common_region_t.test")
plot_region_mean2<-P_mean_site+ylim(0,110) + stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.format",label.y = 100) +labs(title = "un_common_region_t.test")

plot_region_mean3<-P_mean_site +ylim(0,110) + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif",label.y = 100) +labs(title = "un_common_region_wilcox.test")
plot_region_mean4<-P_mean_site +ylim(0,110) + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format",label.y = 100) +labs(title = "un_common_region_wilcox.test")
plot_region_mean1_1<-plot_region_mean1+facet_wrap(~Sample.types,scales="free",nrow=3)
plot_region_mean2_1<-plot_region_mean2+facet_wrap(~Sample.types,scales="free",nrow=3)
plot_region_mean3_1<-plot_region_mean3+facet_wrap(~Sample.types,scales="free",nrow=3)
plot_region_mean4_1<-plot_region_mean4+facet_wrap(~Sample.types,scales="free",nrow=3)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_all_bin3.pdf",plot_region_mean1_1,width=8, height=15)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_all_bin4.pdf",plot_region_mean2_1,width=8, height=15)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_all_bin5.pdf",plot_region_mean3_1,width=8, height=15)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_all_bin6.pdf",plot_region_mean4_1,width=8, height=15)



##using common region

mean_data_common<-read.csv("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/AMA_all_60sample_6X_CpG_DNA_methylation_in_target_region_common_bin.csv",row.names = 1)
head(mean_data_common,n=7)
dim(mean_data_common)
mean_data_common$group_age<-paste(mean_data_common$Sample.types,mean_data_common$Age_group,sep="-")
plot_mean_data2<-mean_data_common[which(mean_data_common$region %in% target_region_list2),]
dim(plot_mean_data2)
head(plot_mean_data2)
plot_mean_data2$sample<-factor(plot_mean_data2$sample,levels=c(kids_list,mother_list,father_list))
plot_mean_data2$analysis_name<-factor(plot_mean_data2$analysis_name,levels=c(paste0("YOUNG_K_",1:10),paste0("AMA_K_",1:10),paste0("YOUNG_M_",1:10),paste0("AMA_M_",1:10),paste0("YOUNG_F_",1:10),paste0("AMA_F_",1:10)))
plot_mean_data2$group_age<-factor(plot_mean_data2$group_age,levels=c("kids-YOUNG","kids-AMA","mother-YOUNG","mother-AMA","father-YOUNG","father-AMA"))

plot_mean_data2$group_region<-factor(plot_mean_data$group_region,levels=c("Low_complexity-YOUNG","Low_complexity-AMA","CpGislands-YOUNG","CpGislands-AMA",
                                                                         "HCP-YOUNG","HCP-AMA","ICP-YOUNG","ICP-AMA","LCP-YOUNG","LCP-AMA", "UTR5-YOUNG","UTR5-AMA",
                                                                         "UTR3-YOUNG","UTR3-AMA", "All_Genes-YOUNG","All_Genes-AMA",
                                                                         "LTR-YOUNG","LTR-AMA", "LINE-YOUNG","LINE-AMA", "SINE-YOUNG","SINE-AMA",
                                                                         "Retroposon_SVA-YOUNG","Retroposon_SVA-AMA","Transposon-YOUNG","Transposon-AMA",
                                                                         "Satellite-YOUNG","Satellite-AMA","Micro_satellites-YOUNG","Micro_satellites-AMA",
                                                                         "MEs-YOUNG","MEs-AMA","ICRs_67-YOUNG","ICRs_67-AMA"))
plot_mean_data2$group_region
head(plot_mean_data2)    
P_mean_bar_1<-ggplot(data=plot_mean_data2,aes(x=group_region,y= meth_mean))+
  geom_bar(aes(fill=sample),stat = "identity", position ="dodge", width=0.6)+# position_dodge2(padding = 0.1,preserve = "single")
  scale_fill_manual(values=c(rep(ppCor[7],10),rep(ppCor[8],10),rep(ppCor[2],10),rep(ppCor[1],10),rep(ppCor[3],10),rep(ppCor[5],10)))+
  labs(title = "common region and CpGs served")+xlab("") +ylab("DNA methylation level(mean)")+#coord_flip()+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))
P_mean_bar_2<-ggplot(data=plot_mean_data2,aes(x=group_region,y= meth_mean))+
  geom_bar(aes(fill=analysis_name),stat = "identity", position ="dodge", width=0.6)+# position_dodge2(padding = 0.1,preserve = "single")
  scale_fill_manual(values=c(rep(ppCor[7],10),rep(ppCor[8],10),rep(ppCor[2],10),rep(ppCor[1],10),rep(ppCor[3],10),rep(ppCor[5],10)))+
  labs(title = "common region and CpGs served")+xlab("") +ylab("DNA methylation level(mean)")+#coord_flip()+
  theme(plot.title = element_text(hjust=0.5,size=10,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))
P_mean_bar_1_1<-P_mean_bar_1+ facet_grid(.~Sample.types)
P_mean_bar_2_1<-P_mean_bar_2+ facet_grid(.~Sample.types)
P_mean_bar_1_2<-P_mean_bar_1+ facet_wrap(~Sample.types,scales="free",nrow=3)
P_mean_bar_2_2<-P_mean_bar_2+ facet_wrap(~Sample.types,scales="free",nrow=3)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/bar_meth_means_in_each_element_for_six_groups_common_bin.pdf",P_mean_bar_1_1,width=18, height=6)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/bar_meth_means_in_each_element_for_six_groups_common_bin2.pdf",P_mean_bar_2_1,width=18, height=6)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/bar_meth_means_in_each_element_for_six_groups_common_bin3.pdf",P_mean_bar_1_2,width=10, height=12)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/bar_meth_means_in_each_element_for_six_groups_common_bin4.pdf",P_mean_bar_2_2,width=10, height=12)

P_mean_site <- ggboxplot(plot_mean_data2, x = "group_region", y = "meth_mean",color = "Age_group", palette = ppCor[1:2],add = "jitter")+#ylim(c(0,90))+
  xlab("sample_group") +ylab("DNA_methylation")+#title("common region and CpGs served")+
  theme(plot.title = element_text(hjust=0.5,size=5,vjust=0.5),axis.text.x=element_text(angle=90,hjust=1, vjust=0.5),axis.line = element_line(colour="black"))+
  theme(#plot.title = element_text(size=10,colour = "blue",face = "bold"),axis.title.x = element_text(size=10,colour = "black",face = "bold"),
    axis.title.y = element_text(size=10,colour = "black",face = "bold"),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10,colour = "black"))
P_mean_site1<-P_mean_site+facet_wrap(~Sample.types,scales="free",nrow=3)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_common_bin.pdf",P_mean_site1,width=16, height=30)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_common_bin2.pdf",P_mean_site1,width=8, height=15)

my_comparisons <- list(c("Low_complexity-YOUNG","Low_complexity-AMA"),c("CpGislands-YOUNG","CpGislands-AMA"),
                       c("HCP-YOUNG","HCP-AMA"),c("ICP-YOUNG","ICP-AMA"), c("LCP-YOUNG","LCP-AMA"), 
                       c("UTR5-YOUNG","UTR5-AMA"),c("UTR3-YOUNG","UTR3-AMA"), c("All_Genes-YOUNG","All_Genes-AMA"),
                       c("LTR-YOUNG","LTR-AMA"),c( "LINE-YOUNG","LINE-AMA"),c("SINE-YOUNG","SINE-AMA"),
                       c("Retroposon_SVA-YOUNG","Retroposon_SVA-AMA"),c("Transposon-YOUNG","Transposon-AMA"),
                       c("Satellite-YOUNG","Satellite-AMA"),c("Micro_satellites-YOUNG","Micro_satellites-AMA"),
                       c("MEs-YOUNG","MEs-AMA"),c("ICRs_67-YOUNG","ICRs_67-AMA"))

#P_mean_site2 +#stat_compare_means(method = "kruskal.test",aes(label=paste0(..method..,": "," p = ",..p.format..)), label.x = 1.4, label.y =85)+
#P_mean_site2 +#stat_compare_means(method = "kruskal.test",aes(label=paste0(..method..,": "," p = ",..p.format..)), label.x = 1.4, label.y =85)+
plot_region_mean1<-P_mean_site+ylim(0,110) + stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.signif",label.y = 100) +labs(title = "un_common_region_t.test")
plot_region_mean2<-P_mean_site+ylim(0,110) + stat_compare_means(method = "t.test", comparisons = my_comparisons, label = "p.format",label.y = 100) +labs(title = "un_common_region_t.test")

plot_region_mean3<-P_mean_site +ylim(0,110) + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif",label.y = 100) +labs(title = "un_common_region_wilcox.test")
plot_region_mean4<-P_mean_site +ylim(0,110) + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format",label.y = 100) +labs(title = "un_common_region_wilcox.test")
plot_region_mean1_1<-plot_region_mean1+facet_wrap(~Sample.types,scales="free",nrow=3)
plot_region_mean2_1<-plot_region_mean2+facet_wrap(~Sample.types,scales="free",nrow=3)
plot_region_mean3_1<-plot_region_mean3+facet_wrap(~Sample.types,scales="free",nrow=3)
plot_region_mean4_1<-plot_region_mean4+facet_wrap(~Sample.types,scales="free",nrow=3)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_common_bin3.pdf",plot_region_mean1_1,width=8, height=15)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_common_bin4.pdf",plot_region_mean2_1,width=8, height=15)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_common_bin5.pdf",plot_region_mean3_1,width=8, height=15)
ggsave("/mnt/data/chenwei/huahua/4.methy_result/2.meth_target_region/box_meth_means_in_each_element_for_six_groups_common_bin6.pdf",plot_region_mean4_1,width=8, height=15)

